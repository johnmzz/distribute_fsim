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

            if (partition_method == "random") {
                MPI_Status status;
                DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type);
                // prepare to initialize sim matrix and allocate values by naively divide into chunk
                //float value;
                //uint32_t num_workers, avg_row, rows, extra, count;
                //vector<uint32_t> count_per_partition(num_tasks);

                //num_workers = num_tasks - 1;
                //avg_row = graph_1.get_num_vertices() / num_workers;
                //extra = graph_1.get_num_vertices() % num_workers;
                //cout << "Rows per partition: " << avg_row << ", extra: " << extra << endl;

                //VertexType g1_vertex_labels = graph_1.get_vertex_type();
                //VertexType g2_vertex_labels = graph_2.get_vertex_type();

                // initialize sim matrix based on label values
                //vector<vector<float>> sim_matrix(graph_1.get_num_vertices(), vector<float>(graph_2.get_num_vertices(), 0));
                //for (int i = 0; i < graph_1.get_num_vertices(); i++) {
                //    for (int j = 0; j < graph_2.get_num_vertices(); j++) {
                //        sim_matrix[i][j] = initializer.get_sim_value(g1_vertex_labels[i], g2_vertex_labels[j]);
                //    }
                //}

                //int offset = 0;
                //for (int worker = 1; worker <= num_workers; worker++) {
                //    rows = (worker <= extra) ? avg_row + 1 : avg_row;
                //    MPI_Send(&offset, 1, MPI_INT, worker, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
                //    MPI_Send(&rows, 1, MPI_INT, worker, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
                //    for (int i = offset; i < sim_matrix.size(); i++) {
                //        MPI_Send(&sim_matrix[i][0], graph_2.get_num_vertices(), MPI_FLOAT, worker, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
                //    }
                //    offset += rows;
                //}
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

                    // out-neighbors
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
            }
            /*
            // create similarity matrix
            SimMatrix premat;
            for (uint32_t i = 0; i < graph_1.get_num_vertices(); ++i) {
                premat.push_back(vector<float>(graph_2.get_num_vertices(), 0.0));
            }

            // initialize similarity matrix
            FracSimulation fracsim(graph_1, graph_2, std::move(initializer), std::move(match_method), std::move(denominator));
            fracsim.simmat_copy_to(premat);
            fracsim.initialize(num_thus);

            // perform iteration
            int turn = 1;
            do {
                cout << "***********************round : " << turn << "**********************" << endl;
                fracsim.simmat_copy_to(premat);
                fracsim.update_simmatrix(premat, w_i, w_o, w_l, num_thus);
                turn += 1;
            } while (!fracsim.diff_relative(premat, EPS) && turn <= MAXTURN);

            end_total = MPI_Wtime();
            cout << "Total computation time : " << end_total - start_total << " s." << endl;

            //fracsim.save_simmat_to_bin(result_path, config_summary);
            fracsim.print_simmat();
            cout << "*******************************************************" << endl;
            */
        } else {
            cout << "Wrong Optimization Code!" << endl;
            exit(0);
        }
    }

    // Step 5: worker task
    if (task_id > 0) {
        if (partition_method == "random") {
            MPI_Status status;
            int offset, rows, extra, num_workers, avg_row;

            num_workers = num_tasks - 1;

            //printf("Worker processor %s, rank %d out of %d processors\n", processor_name, task_id, num_tasks);

            MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);

            vector<float> temp(graph_2.get_num_vertices());
            vector<vector<float>> sim_matrix(rows, vector<float>(graph_2.get_num_vertices()));

            for (int i = 0; i < rows; i++) {
                MPI_Recv(&temp.front(), graph_2.get_num_vertices(), MPI_FLOAT, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
                sim_matrix[i] = temp;
            }
            /*
            if (task_id == 3) {
                cout << "Hello 3: " << endl;
                for (auto i : sim_matrix) {
                    for (auto j : i) {
                        cout << j << " ";
                    }
                    cout << endl;
                }
            }

            cout << endl;
            */
        }

        if (partition_method == "hdrf") {
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type);
            dependency_graph.worker_receive_partition(graph_1, graph_2);
            dependency_graph.worker_receive_partition_record(graph_2);
            dependency_graph.worker_initialize(graph_1, graph_2);

            int copy_values;
            while (true) {
                dependency_graph.worker_copy_values();
                dependency_graph.worker_compute(graph_1, graph_2, w_i, w_o, w_l, label_constrainted);
            }
            /*
            // printf("Worker processor %s, rank %d out of %d processors\n", processor_name, task_id, num_tasks);
            MPI_Status status;
            int num_workers = num_tasks - 1;
            int64_t count = 0;

            vector<unordered_map<uint32_t, float>> sim_values;
            sim_values.resize(graph_1.get_num_vertices());

            // vector<unordered_map> edges;

            int64_t from;
            int64_t to;

            while (true) {
                MPI_Recv(&from, 1, MPI_INT64_T, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
                if (from == -1) {
                    break;
                }

                MPI_Recv(&to, 1, MPI_INT64_T, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);

                uint32_t from_1 = from / graph_2.get_num_vertices();
                uint32_t from_2 = from % graph_2.get_num_vertices();
                uint32_t to_1 = to / graph_2.get_num_vertices();
                uint32_t to_2 = to % graph_2.get_num_vertices();

                if (sim_values[from_1].find(from_2) == sim_values[from_1].end()) {
                    sim_values[from_1].insert({from_2, 0.0});
                    count++;
                }
                if (sim_values[to_1].find(to_2) == sim_values[to_1].end()) {
                    sim_values[to_1].insert({to_2, 0.0});
                    count++;
                }

                if (task_id == 4) {
                    // printf("Worker processor %s, rank %d out of %d processors, received from master (%d,%d)\n", processor_name, task_id, num_tasks, from, to);
                }
            }
            if (task_id == 4) {
                // printf("Worker processor %s, rank %d out of %d processors finished receiving %d node-pairs.\n", processor_name, task_id, num_tasks, count);
            }
            /*
            if (task_id == 4) {
                cout << sim_values.size() << endl;
                for (int i = 0; i < sim_values.size(); i++) {
                    cout << i << ": ";
                    for (auto v : sim_values[i]) {
                        cout << "(" << v.first << "," << v.second << "), ";
                    }
                    cout << endl;
                }
            }


            MPI_Send(&count, 1, MPI_INT64_T, MASTER, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD);

            int start_initialize;
            MPI_Recv(&start_initialize, 1, MPI_INT, MASTER, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD, &status);
            printf("Worker processor %s, rank %d out of %d processors starting to initialize sim values.\n", processor_name, task_id, num_tasks);
            */
        }

        /*
        uint32_t num_workers, , extra;
        float value;
        vector<float> values;
        vector<float> label_sim;

        num_workers = num_tasks - 1;
        v_per_partition = graph_1.get_num_vertices() * graph_2.get_num_vertices() / num_workers;
        extra = graph_1.get_num_vertices() * graph_2.get_num_vertices() % num_workers;

        vector<vector<int>> values_node_pair(v_per_partition + extra, vector<int>(2, -1));

        int x, y;
        int i = 0;
        // poll for values
        while (i < (v_per_partition + extra)) {
            MPI_Recv(&value, 1, MPI_FLOAT, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
            MPI_Recv(&x, 1, MPI_INT, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
            MPI_Recv(&y, 1, MPI_INT, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
            if (value >= 0) {
                values.push_back(value);
                values_node_pair[i][0] = x;
                values_node_pair[i][1] = y;
            }
            i++;
        }
        values_node_pair.resize(values.size());
        label_sim = values;

        // receive partition allocation information
        // int a[13];
        // vector<vector<int>> partition_record(graph_1.get_num_vertices(), vector<int>(graph_2.get_num_vertices()));
        // for (int j = 0; j < graph_1.get_num_vertices(); j++) {
        //    MPI_Recv(&partition_record[j].front(), partition_record[j].size(), MPI_INT, MASTER, FROM_MASTER_PARTITION, MPI_COMM_WORLD, &status);
        //}

        // printf("Worker processor %s, rank %d out of %d processors, has %d values\n", processor_name, task_id, num_tasks, values.size());

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