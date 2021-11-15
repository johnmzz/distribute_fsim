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

#include "DependencyGraph.h"
#include "Graph.h"
#include "InitializerMPI.h"

using namespace std;
using namespace chrono;

#define EPS 0.01
#define MAXTURN 20

int main(int argc, char* argv[]) {
    // Read parameter configuration
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

    // Read graph
    GraphNL graph_1(dataset_path, input_graph_1);
    GraphNL graph_2(dataset_path, input_graph_2);

    // Initialize MPI communication environment
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);

    // Master task
    if (task_id == 0) {
        string config_summary = input_graph_1 + "_" + input_graph_2 + "_" + initializer_type + "_" + optimization + "_" + simulation_type + "_" + label_constrainted + "_i" + to_string(w_i) + "_o" + to_string(w_o) + "_l" + to_string(w_l) + "_t" + to_string(theta) + "_thr" + to_string(num_thus);
        cout << "Configuration Summaries: " << config_summary << endl;
        cout << "Start computation, initializer type: " << initializer_type << ", simulation type: " << simulation_type << ", optimization: " << optimization << endl;
        double start_total, end_total;
        start_total = MPI_Wtime();

        InitializerMPI initializer(initializer_type, graph_1.get_label_info(), graph_2.get_label_info());

        // ****************************** random partitioning ******************************
        if (partition_method == "random") {
            int signal = 1;
            MPI_Status status;

            // step 1: partition sim_matrix
            int64_t total_pairs = 0;
            for (int p = 1; p <= num_partitions; p++) {
                int64_t count;
                MPI_Recv(&count, 1, MPI_INT64_T, p, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD, &status);
                cout << "Worker " << p - 1 << ", task id = " << p << " has received its partition of size = " << count << "." << endl;
                total_pairs += count;
            }
            cout << "Total node-pair assigned = " << total_pairs << endl << endl;

            // step 2: initialize each worker's sim_values
            cout << "Start initialization..." << endl;
            for (int p = 1; p <= num_partitions; p++) {
                MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD);
            }
            for (int p = 1; p <= num_partitions; p++) {
                MPI_Recv(&signal, 1, MPI_INT, p, WORKER_INIT_COMPLETE, MPI_COMM_WORLD, &status);
            }
            cout << "Workers initialization completed." << endl << endl;

            // step 3: computation (iteration)
            cout << "Starting computation." << endl;
            int turn = 1;
            do {
                cout << "***********************round : " << turn << "**********************" << endl;

                // step 3-1: worker copy sim_values before computation
                cout << "Workers copy sim values" << endl;
                for (int p = 1; p <= num_partitions; p++) {
                    MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_COPY, MPI_COMM_WORLD);
                }
                for (int p = 1; p <= num_partitions; p++) {
                    MPI_Recv(&signal, 1, MPI_INT, p, WORKER_COPY_COMPLETE, MPI_COMM_WORLD, &status);
                }
                cout << "Workers have copied their sim values." << endl << endl;

                // step 3-2: workers calculation one by one
                for (int round = 1; round <= num_partitions; round++) {
                    cout << "Worker " << round << " performing computation..." << endl;
                    for (int p = 1; p <= num_partitions; p++) {
                        MPI_Send(&round, 1, MPI_INT, p, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD);
                    }
                    MPI_Recv(&signal, 1, MPI_INT, round, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);  // the computing worker finishes
                    for (int p = 1; p <= num_partitions; p++) {
                        if (p != round) {
                            uint32_t stop_listenner[2] = {0, 0};
                            MPI_Send(&stop_listenner[0], 2, MPI_UINT32_T, p, QUERY_VALUE, MPI_COMM_WORLD);
                            MPI_Recv(&signal, 1, MPI_INT, p, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);
                        }
                    }
                    cout << "Worker " << round << " finished its computation." << endl;
                }
                turn += 1;
            } while (turn <= 7);
            cout << "All computation finished." << endl << endl;

            // step 4: output sim_values
            cout << "Final result:" << endl;
            for (int p = 1; p <= num_partitions; p++) {
                MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_PRINT, MPI_COMM_WORLD);
                MPI_Recv(&signal, 1, MPI_INT, p, WORKER_FINISH_PRINT, MPI_COMM_WORLD, &status);
            }
        // ******************************** hdrf partitioning ********************************
        } else if (partition_method == "hdrf") {
            int signal = 1;
            MPI_Status status;
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type, partition_method);

            // step 1: initialize dependency graph
            cout << "Initializing dependency graph..." << endl;
            dependency_graph.init_dependency_graph(graph_1, graph_2);
            cout << "Dependency graph initialization completed." << endl << endl;

            // step 2: HDRF partition
            cout << "Partitioning using HDRF algorithm..." << endl;
            dependency_graph.HDRF_partition(graph_1, graph_2);
            cout << "Partitioning complete." << endl;
            int64_t replication_count = 0;
            for (int p = 1; p <= num_partitions; p++) {
                int64_t temp;
                MPI_Recv(&temp, 1, MPI_INT64_T, p, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD, &status);
                cout << "Worker task id = " << p << " has received its partition." << endl;
                replication_count += temp;
            }
            cout << "Total node-pair assigned = " << replication_count << endl << endl;

            // step 3: distribute partition records across workers
            cout << "Sending workers partition record..." << endl;
            dependency_graph.distribute_partition_record();
            for (uint32_t p = 1; p <= num_partitions; p++) {
                MPI_Recv(&signal, 1, MPI_INT, p, WORKER_RECEIVED_PART_RECORD, MPI_COMM_WORLD, &status);
            }
            cout << "All workers received their partition records." << endl << endl;

            // step 4: initialization on all workers
            cout << "Start initialization..." << endl;
            for (uint32_t p = 1; p <= num_partitions; p++) {
                MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD);
            }
            for (uint32_t p = 1; p <= num_partitions; p++) {
                MPI_Recv(&signal, 1, MPI_INT, p, WORKER_INIT_COMPLETE, MPI_COMM_WORLD, &status);
            }
            cout << "Workers initialization completed." << endl;

            // step 5: computation (iteration)
            cout << "Starting computation." << endl;
            int turn = 1;
            do {
                cout << "***********************round : " << turn << "**********************" << endl;

                // step 5-1: worker copy sim_values before computation
                cout << "Workers copy sim values" << endl;
                for (int p = 1; p <= num_partitions; p++) {
                    MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_COPY, MPI_COMM_WORLD);
                }
                for (int p = 1; p <= num_partitions; p++) {
                    MPI_Recv(&signal, 1, MPI_INT, p, WORKER_COPY_COMPLETE, MPI_COMM_WORLD, &status);
                }
                cout << "Workers have copied their sim values." << endl << endl;

                // step 5-2: workers calculation one by one
                for (int round = 1; round <= num_partitions; round++) {
                    cout << "Worker " << round << " performing computation..." << endl;
                    for (int p = 1; p <= num_partitions; p++) {
                        MPI_Send(&round, 1, MPI_INT, p, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD);
                    }
                    MPI_Recv(&signal, 1, MPI_INT, round, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);  // the computing worker finishes
                    for (int p = 1; p <= num_partitions; p++) {
                        if (p != round) {
                            uint32_t stop_listenner[2] = {0, 0};
                            MPI_Send(&stop_listenner[0], 2, MPI_UINT32_T, p, QUERY_PARTITION_RECORD, MPI_COMM_WORLD);
                            MPI_Recv(&signal, 1, MPI_INT, p, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD, &status);
                        }
                    }
                    cout << "Worker " << round << " finished its computation." << endl;
                }
                turn += 1;
            } while (turn <= 6);
            cout << "All computation finished." << endl << endl;

            // step 6: output sim_values
            cout << "Final result:" << endl;
            for (int p = 1; p <= num_partitions; p++) {
                MPI_Send(&signal, 1, MPI_INT, p, FROM_MASTER_PRINT, MPI_COMM_WORLD);
                MPI_Recv(&signal, 1, MPI_INT, p, WORKER_FINISH_PRINT, MPI_COMM_WORLD, &status);
            }

        }
    }

    if (task_id > 0) {
        if (partition_method == "random") {
            MPI_Status status;
            int signal = 1;
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type, partition_method);

            // step 1: partition sim_matrix
            uint64_t count = dependency_graph.worker_random_partition(graph_1, graph_2);
            MPI_Send(&count, 1, MPI_INT64_T, MASTER, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD);

            // step 2: worker initialize sim_values
            MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD, &status);
            dependency_graph.worker_initialize(graph_1, graph_2);
            MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_INIT_COMPLETE, MPI_COMM_WORLD);

            // step 3: computation (iteration)
            int turn = 1;
            do {
                // step 3-1: copy sim_values before computation
                MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_COPY, MPI_COMM_WORLD, &status);
                dependency_graph.worker_copy_values();
                MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COPY_COMPLETE, MPI_COMM_WORLD);

                // step 3-2: calculation or listenning
                for (int round = 1; round <= num_partitions; round++) {
                    int p;
                    MPI_Recv(&p, 1, MPI_INT, MASTER, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD, &status);
                    if (p == task_id) {
                        dependency_graph.worker_random_compute(graph_1, graph_2, w_i, w_o, w_l, label_constrainted);
                        MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
                    } else {
                        dependency_graph.worker_listen();
                        MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
                    }
                }
                turn += 1;
            } while (turn <= 7);

            // step 4: output result
            MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_PRINT, MPI_COMM_WORLD, &status);
            dependency_graph.worker_print();
            dependency_graph.write_result_txt();
            MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_FINISH_PRINT, MPI_COMM_WORLD);
        }

        if (partition_method == "hdrf") {
            MPI_Status status;
            int signal = 1;
            DependencyGraph dependency_graph(graph_1, graph_2, num_partitions, initializer_type, simulation_type, partition_method);

            // step 1: initialize dependency graph (done by master)

            // step 2: receive partition
            uint64_t count = dependency_graph.worker_receive_partition(graph_1, graph_2);
            MPI_Send(&count, 1, MPI_INT64_T, MASTER, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD);

            // step 3: receive partition record
            dependency_graph.worker_receive_partition_record(graph_2);
            MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_RECEIVED_PART_RECORD, MPI_COMM_WORLD);

            // step 4: initialize sim value
            MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD, &status);
            dependency_graph.worker_initialize(graph_1, graph_2);
            MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_INIT_COMPLETE, MPI_COMM_WORLD);

            // step 5: computation (iteration)
            int turn = 1;
            do {
                // step 5-1: copy sim_values before computation
                MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_COPY, MPI_COMM_WORLD, &status);
                dependency_graph.worker_copy_values();
                MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COPY_COMPLETE, MPI_COMM_WORLD);

                // step 5-2: calculation or listenning
                for (int round = 1; round <= num_partitions; round++) {
                    int p;
                    MPI_Recv(&p, 1, MPI_INT, MASTER, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD, &status);
                    if (p == task_id) {
                        dependency_graph.worker_random_compute(graph_1, graph_2, w_i, w_o, w_l, label_constrainted);
                        MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
                    } else {
                        //dependency_graph.worker_hdrf_listen();
                        MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
                    }
                }
                turn += 1;
            } while (turn <= 6);

            // step 6: output result
            MPI_Recv(&signal, 1, MPI_INT, MASTER, FROM_MASTER_PRINT, MPI_COMM_WORLD, &status);
            dependency_graph.worker_print();
            MPI_Send(&signal, 1, MPI_INT, MASTER, WORKER_FINISH_PRINT, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}