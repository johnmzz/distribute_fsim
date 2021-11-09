#pragma once

#include <mpi.h>
#include <omp.h>

#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Graph.h"
#include "InitializerMPI.h"

#define MASTER 0
#define FROM_MASTER_PARTITION 1
#define WORKER_RECEIVED_PARTITION 2
#define FROM_MASTER_PART_RECORD 3
#define WORKER_RECEIVED_PART_RECORD 4
#define FROM_MASTER_INITIALIZE 5
#define WORKER_INIT_COMPLETE 6
#define FROM_MASTER_COPY 7
#define WORKER_COPY_COMPLETE 8
#define FROM_MASTER_START_COMPUTE 9
#define WORKER_COMPUTE_COMPLETE 10
#define FROM_MASTER_SYNCHRONIZE 11
#define WORKER_SYNCHRONIZE_COMPLETE 12
#define FROM_MASTER_AGGRE 13
#define WORKER_AGGRE_COMPLETE 14
#define QUERY_PARTITION_RECORD 15
#define ANSWER_PARTITION_RECORDD 16
#define QUERY_VALUE 17
#define ANSWER_VALUE 18

using namespace std;

class DependencyGraph {
  private:
    double lambda;
    const double BALANCE_RATIO_ = 1.0;
    double capacity_;
    int64_t dependency_count_;  // number of edges in dependency graph
    int64_t node_pairs_count_;  // number of total node-pairs across the two input graphs
    uint32_t num_partitions_;
    int64_t max_size_;
    int64_t min_size_;
    string simulation_type_;
    InitializerMPI initializer_;
    vector<int64_t> occupied_;  // number of edges in each partition, eg. [p0_num_edges, p1_num_edges, ...]
    vector<int64_t> count_;     // degree of each node-pair (in dependency graph, not original graph), eg. [d(v0-v0, d(v0-v1), ...)]  TODO: distribute?
    vector<boost::dynamic_bitset<>> is_boundarys_;
    vector<unordered_map<uint32_t, float>> sim_values_;
    vector<unordered_map<uint32_t, float>> pre_values_;
    vector<unordered_map<uint32_t, boost::dynamic_bitset<>>> part_record_;
    vector<uint32_t> node_label_1_;
    vector<uint32_t> node_label_2_;

    int64_t encode(uint32_t i, uint32_t j, uint32_t offset);
    uint32_t HDRF_calculate_partition(uint32_t from, uint32_t to);
    double compute_partition_score(uint32_t u, uint32_t v, uint32_t i);

  public:
    DependencyGraph(GraphNL& graph_1, GraphNL& graph_2, uint32_t num_partitions, string initializer_type, string simulation_type);
    void init_dependency_graph(GraphNL& graph_1, GraphNL& graph_2);
    void HDRF_partition(GraphNL& graph_1, GraphNL& graph_2);
    void get_vertex_partition(string graph_1_name, string graph_2_name);
    void count_replication();
    void print_label_mat();
    void worker_receive_partition(GraphNL& graph_1, GraphNL& graph_2);
    void worker_initialize(GraphNL& graph_1, GraphNL& graph_2);
    void distribute_partition_record();
    void worker_random_partition(GraphNL& graph_1, GraphNL& graph_2);
    void worker_receive_partition_record(GraphNL& graph_2);
    void worker_copy_values();
    void worker_compute(GraphNL& graph_1, GraphNL& graph_2, float w_i, float w_o, float w_l, string label_constrainted);
    void worker_random_compute(GraphNL& graph_1, GraphNL& graph_2, float w_i, float w_o, float w_l, string label_constrainted);
};

// eg.
// is_boundary_ = [
//                  [0000000010010100101001],     // vertices in bucket 1
//                  [1100000000001011110111],     // vertices in bucket 2
//                  ...
//                ]
//
// complexity = p * v^2 / 8