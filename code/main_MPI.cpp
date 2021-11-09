#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>

#include "Denominator.h"
#include "DependencyGraph.h"
#include "FracSimulation.h"
#include "Graph.h"
#include "Initializer.h"
#include "InitializerMPI.h"
#include "MatchMethod.h"

using namespace std;
using namespace chrono;

#define EPS 0.01
#define MAXTURN 20

int main(int argc, char* argv[]) {
    // Step 1: read parameter configuration
    // sample command: mpiexec -n 7 -f host_file ./calsim dataset/ results/ mz_test mz_test identical null sim N 0.4 0.4 0.2 0 32 hdrf 6
    // sample command: mpiexec -n 7 -f host_file ./calsim /import/vldb01/2/scratch/mazhuo/HIN_data results/ yeast yeast identical null sim N 0.4 0.4 0.2 0 32 hdrf 6
    string dataset_path = argv[1];
    string result_path = argv[2];
    string input_graph_1 = argv[3];
    string input_graph_2 = argv[4];
    string initializer_type = argv[5];
    string optimization = argv[6];
    string simulation_type = argv[7];
    string label_constrainted = argv[8];
    float w_o = stof(argv[9]);
    float w_i = stof(argv[10]);
    float w_l = stof(argv[11]);
    float theta = stof(argv[12]);
    uint32_t num_thus = stoi(argv[13]);
    string partition_method = argv[14];
    uint32_t num_partitions = stoi(argv[15]);

    // Step 2: read graph
    GraphNL graph_1(dataset_path, input_graph_1);
    GraphNL graph_2(dataset_path, input_graph_2);

    // Step 3: Initialize MPI communication environment
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);

    // Step 4: master task
    if (task_id == 0) {
        string config_summary = input_graph_1 + "_" + input_graph_2 + "_" + initializer_type + "_" + optimization + "_" + simulation_type + "_" + label_constrainted + "_i" + to_string(w_i) + "_o" + to_string(w_o) + "_l" + to_string(w_l) + "_t" + to_string(theta) + "_thr" + to_string(num_thus);
        cout << "Configuration Summaries: " << config_summary << endl;

        // Step 4-1: chose hyper-parameter based on config
        unique_ptr<MatchMethod> match_method;
        unique_ptr<Denominator> denominator;

        if (simulation_type == "sim") {
            match_method = make_unique<InjectMatch>(label_constrainted);
            denominator = make_unique<QuerySide>();
        } else if (simulation_type == "dpsim") {
            match_method = make_unique<BijectMatch>(label_constrainted);
            denominator = make_unique<QuerySide>();
        } else if (simulation_type == "bisim") {
            match_method = make_unique<BiMatch>(label_constrainted);
            denominator = make_unique<Add>();
        } else if (simulation_type == "bijectsim") {
            match_method = make_unique<BijectMatch>(label_constrainted);
            denominator = make_unique<Root>();
        }

        // Step 4-2: perform calculation
        if (optimization == "NULL" || optimization == "null") {
            cout << "Start computation, initializer type: " << initializer_type << ", simulation type: " << simulation_type << ", optimization: " << optimization << endl;

            double start_total, end_total;
            start_total = MPI_Wtime();

            InitializerMPI initializer(initializer_type, graph_1.get_label_info(), graph_2.get_label_info());

            // ****************************** random partitioning ******************************
            if (partition_method == "random") {
                MPI_Status status;

                int64_t total_pairs = 0;
                for (int i = 1; i <= num_partitions; i++) {
                    int64_t temp;
                    MPI_Recv(&temp, 1, MPI_INT64_T, i, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD, &status);
                    cout << "Worker " << i - 1 << ", task id = " << i << " has received its partition of size = " << temp << "." << endl;
                    total_pairs += temp;
                }
                cout << "Total node-pair assigned = " << total_pairs << endl;

                // values initialization on all workers
                cout << "\nStart initialization..." << endl;
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int temp;
                    MPI_Send(&temp, 1, MPI_INT, i, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD);
                }
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int initialization_complete;
                    MPI_Recv(&initialization_complete, 1, MPI_INT, i, WORKER_INIT_COMPLETE, MPI_COMM_WORLD, &status);
                }
                cout << "Workers initialization completed." << endl;

                cout << "\nStarting computation." << endl;
                int turn = 1;
                do {
                    cout << "***********************round : " << turn << "**********************" << endl;

                    cout << "Workers copy sim values" << endl;
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Send(&temp, 1, MPI_INT, i, FROM_MASTER_COPY, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_COPY_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    cout << "Workers have copied their sim values." << endl;

                    cout << "Workers start computing..." << endl;
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int start_compute = 1;
                        MPI_Send(&start_compute, 1, MPI_INT, i, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {    // thread 1 of workers finished computation
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {    // stop thread 2 of workers
                        int start_compute = 1;
                        MPI_Send(&start_compute, 1, MPI_INT, i, QUERY_VALUE, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {    // thread 2 of workers stoped
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    cout << "Workers have finished all local computations." << endl;
                    turn += 1;
                } while (turn <= 5);
                cout << "All computation finished." << endl;

            // ******************************** hdrf partitioning ********************************
            } else if (partition_method == "hdrf") {
                MPI_Status status;
                DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type);

                // initialize dependency graph
                dependency_graph.init_dependency_graph(graph_1, graph_2);

                // perform HDRF partition
                dependency_graph.HDRF_partition(graph_1, graph_2);
                int64_t replication_count = 0;
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int64_t temp;
                    MPI_Recv(&temp, 1, MPI_INT64_T, i, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD, &status);
                    cout << "Worker " << i - 1 << ", task id = " << i << " has received its partition." << endl;
                    replication_count += temp;
                }
                cout << "Total node-pair assigned = " << replication_count << endl;

                // distribute partition records across workers
                dependency_graph.distribute_partition_record();
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int temp;
                    MPI_Recv(&temp, 1, MPI_INT, i, WORKER_RECEIVED_PART_RECORD, MPI_COMM_WORLD, &status);
                }
                cout << "All workers received their partition records." << endl;

                // values initialization on all workers
                cout << "\nStart initialization..." << endl;
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int temp;
                    MPI_Send(&temp, 1, MPI_INT, i, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD);
                }
                for (uint32_t i = 1; i <= num_partitions; i++) {
                    int initialization_complete;
                    MPI_Recv(&initialization_complete, 1, MPI_INT, i, WORKER_INIT_COMPLETE, MPI_COMM_WORLD, &status);
                }
                cout << "Workers initialization completed." << endl;

                cout << "\nStarting computation." << endl;
                int turn = 1;
                do {
                    cout << "***********************round : " << turn << "**********************" << endl;

                    cout << "Workers copy sim values" << endl;
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Send(&temp, 1, MPI_INT, i, FROM_MASTER_COPY, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_COPY_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    cout << "Workers have copied their sim values." << endl;

                    cout << "Workers start computing..." << endl;
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int start_compute = 1;
                        MPI_Send(&start_compute, 1, MPI_INT, i, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    cout << "Workers have finished all local computations." << endl;

                    cout << "Start synchronzation: " << endl;
                    for (int i = 1; i <= num_partitions; i++) {
                        for (uint32_t j = 1; j <= num_partitions; j++) {
                            cout << "Synchronizing worker " << i << "...\n";
                            int synch_number = i;
                            MPI_Send(&synch_number, 1, MPI_INT, j, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD);
                        }
                        for (uint32_t j = 1; j <= num_partitions; j++) {
                            int synch_number;
                            MPI_Recv(&synch_number, 1, MPI_INT, j, WORKER_SYNCHRONIZE_COMPLETE, MPI_COMM_WORLD, &status);
                        }
                        cout << "Worker " << i << " finished synchronization." << endl;
                    }
                    cout << "All workers have finihsed synchronization." << endl;

                    cout << "Workers start final computation..." << endl;
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int start_aggre = 1;
                        MPI_Send(&start_aggre, 1, MPI_INT, i, FROM_MASTER_AGGRE, MPI_COMM_WORLD);
                    }
                    for (uint32_t i = 1; i <= num_partitions; i++) {
                        int temp;
                        MPI_Recv(&temp, 1, MPI_INT, i, WORKER_AGGRE_COMPLETE, MPI_COMM_WORLD, &status);
                    }
                    cout << "Workers have finished all computations." << endl;

                    // TODO: calculate difference

                    turn += 1;
                } while (turn <= 5);
                // } while (!fracsim.diff_relative(premat, EPS) && turn <= MAXTURN);
                end_total = MPI_Wtime();
                cout << "Total computation time: " << end_total - start_total << " s." << endl;
                //fracsim.save_simmat_to_bin(result_path, config_summary);
                //fracsim.print_simmat();
                cout << "*******************************************************" << endl;
            }
        } else {
            cout << "Wrong Optimization Code!" << endl;
            exit(0);
        }
    }

    // Step 5: worker task
    if (task_id > 0) {
        if (partition_method == "random") {
            MPI_Status status;
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type);
            dependency_graph.worker_random_partition(graph_1, graph_2);
            dependency_graph.worker_initialize(graph_1, graph_2);

            int turn = 0;
            while (turn < 5) {
                dependency_graph.worker_copy_values();
                dependency_graph.worker_random_compute(graph_1, graph_2, w_i, w_o, w_l, label_constrainted);
                turn++;
            }
        }

        if (partition_method == "hdrf") {
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type);
            dependency_graph.worker_receive_partition(graph_1, graph_2);
            dependency_graph.worker_receive_partition_record(graph_2);
            dependency_graph.worker_initialize(graph_1, graph_2);

            int turn = 0;
            while (turn < 5) {
                dependency_graph.worker_copy_values();
                dependency_graph.worker_compute(graph_1, graph_2, w_i, w_o, w_l, label_constrainted);
                turn++;
            }
        }
        /*
        uint32_t u, v;
        for (int i = 0; i < values.size(); i++) {
            u = values_node_pair[i][0];
            v = values_node_pair[i][1];
            vector<uint32_t> u_out_neighbors = graph_1.get_out_neighbors_list(u);
            vector<uint32_t> u_in_neighbors = graph_1.get_in_neighbors_list(u);
            vector<uint32_t> v_out_neighbors = graph_2.get_out_neighbors_list(v);
            vector<uint32_t> v_in_neighbors = graph_2.get_in_neighbors_list(v);

            float in_neighbor_sim = 0.0;
            float out_neighbor_sim = 0.0;

            // if both u,v neighbors are empty, or one of them is empty
            if (u_out_neighbors.size() == 0 && v_out_neighbors.size() == 0) {
                out_neighbor_sim = 1.0;
                break;
            } else if (u_out_neighbors.size() == 0 || v_out_neighbors.size() == 0) {
                out_neighbor_sim = 0.0;
                break;
            }
            if (u_in_neighbors.size() == 0 && v_in_neighbors.size() == 0) {
                in_neighbor_sim = 1.0;
                break;
            } else if (u_in_neighbors.size() == 0 || v_in_neighbors.size() == 0) {
                in_neighbor_sim = 0.0;
                break;
            }

            // simple simulation
            if (simulation_type == "sim") {

                vector<uint32_t>::iterator iItr, jItr;
                float numerator = 0.0;
                for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
                    float max_score = 0.0;
                    for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                        weight = basemat[*iItr][*jItr];
                        if (weight > max_score) {
                            max_score = weight;
                        }
                    }
                    numerator += max_score;
                }

            }
            // bisimulation
            else if (simulation_type == "bisim") {
                // ..
                // ..
            }
            // degree-preserving simulation or bijective simulation
            else if (simulation_type == "dpsim" || simulation_type == "bijectsim") {
                // ..
                // ..
            } else {
                cout << "Incorrect simulation type" << endl;
            }
        }
        */
    }

    MPI_Finalize();
    return 0;
}