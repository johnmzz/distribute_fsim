#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <time.h>

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>

#include "Denominator.h"
#include "FracSimulation.h"
#include "Graph.h"
#include "Initializer.h"
#include "MatchMethod.h"

using namespace std;
using namespace chrono;

#define EPS 0.01
#define MAXTURN 20

int main(int argc, char* argv[]) {
    // Step 1: read parameter configuration
    // sample command: mpiexec -n 5 -f host_file ./calsim dataset/ results/ mz_test mz_test identical null sim N 0.4 0.4 0.2 0 32
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

    // Step 2: Initialize MPI communication environment
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);

    // Step 3: master task
    if (task_id == 0) {
        string config_summary = input_graph_1 + "_" + input_graph_2 + "_" + initializer_type + "_" + optimization + "_" + simulation_type + "_" + label_constrainted + "_i" + to_string(w_i) + "_o" + to_string(w_o) + "_l" + to_string(w_l) + "_t" + to_string(theta) + "_thr" + to_string(num_thus);
        cout << "Configuration Summaries: " << config_summary << endl;

        // Step 3-1: read graph
        GraphNL graph_1(dataset_path, input_graph_1);
        GraphNL graph_2(dataset_path, input_graph_2);

        // Step 3-2: chose hyper-parameter based on config
        unique_ptr<MatchMethod> match_method;
        unique_ptr<Denominator> denominator;
        unique_ptr<Initializer> initializer;

        if (initializer_type == "IdenticalInitializer" || initializer_type == "identical") {
            initializer = make_unique<IdenticalInitializer>();
        } else if (initializer_type == "HammingInitializer" || initializer_type == "hamming") {
            initializer = make_unique<HammingInitializer>();
        } else if (initializer_type == "EditDistInitializer" || initializer_type == "edit") {
            initializer = make_unique<EditDistInitializer>();
        } else if (initializer_type == "JarowinklerInitializer" || initializer_type == "jarowinkler") {
            initializer = make_unique<JarowinklerInitializer>();
        } else {
            cout << "Initializer Type Error, " << initializer_type << endl;
            exit(0);
        }

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

        // Step 3-3: perform calculation
        if (optimization == "NULL" || optimization == "null") {
            cout << "Start computation, initializer type: " << initializer_type << ", simulation type: " << simulation_type << ", optimization: " << optimization << endl;

            // create similarity matrix
            double start_n, end_n;
            start_n = omp_get_wtime();
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

            end_n = omp_get_wtime();
            cout << "computation time : " << end_n - start_n << " s." << endl;

            //fracsim.save_simmat_to_bin(result_path, config_summary);
            fracsim.print_simmat();
            cout << "*******************************************************" << endl;
        } else {
            cout << "Wrong Optimization Code!" << endl;
            exit(0);
        }
    }

    // Step 4: worker task
    if (task_id > 0) {
        // printf("Worker processor %s, rank %d out of %d processors\n", processor_name, task_id, num_tasks);
    }

    MPI_Finalize();
    return 0;
}