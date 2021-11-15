#include "DependencyGraph.h"
DependencyGraph::DependencyGraph(GraphNL& graph_1, GraphNL& graph_2, uint32_t num_partitions, string initializer_type, string simulation_type, string partition_method)
    : initializer_(initializer_type, graph_1.get_label_info(), graph_2.get_label_info()) {
    lambda = 1.1;  // default value for HDRF
    capacity_ = 0;
    dependency_count_ = 0;
    node_pairs_count_ = graph_1.get_num_vertices() * graph_2.get_num_vertices();
    num_partitions_ = num_partitions;
    simulation_type_ = simulation_type;
    partition_method_ = partition_method;
    occupied_.assign(num_partitions, 0);
    node_label_1_ = graph_1.get_vertex_type();
    node_label_2_ = graph_2.get_vertex_type();

    int name_len;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks_);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id_);
    MPI_Get_processor_name(processor_name_, &name_len);
    num_workers_ = num_tasks_ - 1;
}

int64_t DependencyGraph::encode(uint32_t i, uint32_t j, uint32_t offset) {
    return i * offset + j;
}

void DependencyGraph::init_dependency_graph(GraphNL& graph_1, GraphNL& graph_2) {
    node_pairs_count_ = graph_1.get_num_vertices() * graph_2.get_num_vertices();
    count_.resize(node_pairs_count_, 0);  // to record degree of each node_pair

    // *********************************************************************************************************************
    // TODO: distribute across partitions
    for (int i = 0; i < num_partitions_; i++) {
        boost::dynamic_bitset<> bs(node_pairs_count_);
        is_boundarys_.push_back(bs);
        // cout << "Partition " << i << " bitset with size: " << is_boundarys_[i].size() << " = " << is_boundarys_[i] << endl;
    }
    // *********************************************************************************************************************

    for (uint32_t i = 0; i < graph_1.get_num_vertices(); i++) {
        auto adj_out_i = graph_1.get_out_neighbors_list(i);
        auto adj_in_i = graph_1.get_in_neighbors_list(i);
        for (uint32_t j = 0; j < graph_2.get_num_vertices(); j++) {
            set<pair<uint32_t, uint32_t>> neighbor_pairs;
            // cout << "For node pair " << i << "," << j << " (" << (i * graph_2.get_num_vertices() + j) << "): ";

            auto adj_out_i = graph_2.get_out_neighbors_list(i);
            auto adj_out_j = graph_2.get_out_neighbors_list(j);

            for (auto n_i : adj_out_i) {
                for (auto n_j : adj_out_j) {
                    pair<uint32_t, uint32_t> neighbors(n_i, n_j);
                    neighbor_pairs.insert(neighbors);
                }
            }

            auto adj_in_i = graph_2.get_in_neighbors_list(i);
            auto adj_in_j = graph_2.get_in_neighbors_list(j);

            for (auto n_i : adj_in_i) {
                for (auto n_j : adj_in_j) {
                    pair<uint32_t, uint32_t> neighbors(n_i, n_j);
                    neighbor_pairs.insert(neighbors);
                }
            }

            for (auto v : neighbor_pairs) {
                int64_t from = encode(i, j, graph_2.get_num_vertices());
                int64_t to = encode(v.first, v.second, graph_2.get_num_vertices());
                // cout << to << ", ";
                dependency_count_++;
                count_[from]++;
                count_[to]++;
            }
            // cout << endl;
        }
        // cout << endl;
        if (i % 1000 == 0) {
            cout << "read line " << i << endl;
        }
    }
    capacity_ = (double)dependency_count_ * BALANCE_RATIO_ / num_partitions_ + 1;

    cout << "Node pairs count = " << node_pairs_count_ << endl;
    cout << "Total dependencies = " << dependency_count_ << endl;
    // cout << "Degree of each node-pair = ";
    // for (auto c : count_) {
    //    cout << c << ", ";
    //}
    //cout << endl;
}

void DependencyGraph::HDRF_partition(GraphNL& graph_1, GraphNL& graph_2) {
    int64_t num_edges = dependency_count_;
    uint32_t bucket;
    max_size_ = 0;
    min_size_ = 0;

    // ofstream fout("dataset/mz_test/" + graph_1.get_graph_name() + "_" + graph_2.get_graph_name() + "_dependency_graph_" + to_string(num_partitions_) + ".txt");
    for (uint32_t i = 0; i < graph_1.get_num_vertices(); i++) {
        auto adj_out_i = graph_1.get_out_neighbors_list(i);
        auto adj_in_i = graph_1.get_in_neighbors_list(i);
        for (uint32_t j = 0; j < graph_2.get_num_vertices(); j++) {
            set<pair<uint32_t, uint32_t>> neighbor_pairs;
            // cout << "For node pair " << i << "," << j << " (" << (i * graph_2.get_num_vertices() + j) << "): ";

            auto adj_out_i = graph_2.get_out_neighbors_list(i);
            auto adj_out_j = graph_2.get_out_neighbors_list(j);

            for (auto n_i : adj_out_i) {
                for (auto n_j : adj_out_j) {
                    pair<uint32_t, uint32_t> neighbors(n_i, n_j);
                    neighbor_pairs.insert(neighbors);
                }
            }

            auto adj_in_i = graph_2.get_in_neighbors_list(i);
            auto adj_in_j = graph_2.get_in_neighbors_list(j);

            for (auto n_i : adj_in_i) {
                for (auto n_j : adj_in_j) {
                    pair<uint32_t, uint32_t> neighbors(n_i, n_j);
                    neighbor_pairs.insert(neighbors);
                }
            }

            for (auto v : neighbor_pairs) {
                int64_t from = encode(i, j, graph_2.get_num_vertices());
                int64_t to = encode(v.first, v.second, graph_2.get_num_vertices());

                bucket = HDRF_calculate_partition(from, to);

                occupied_[bucket]++;
                is_boundarys_[bucket].set(from);
                is_boundarys_[bucket].set(to);

                if (occupied_[bucket] > max_size_) {
                    max_size_ = occupied_[bucket];
                }

                int min_size_bucket_count = 0;
                for (int i = 0; i < num_partitions_; i++) {
                    if (occupied_[i] == min_size_) {
                        min_size_bucket_count++;
                    }
                }
                if (min_size_bucket_count == 0) {
                    min_size_++;
                }

                // fout << from << " " << to << " " << bucket << endl;
                // cout << "max_size = " << max_size_ << ", min_size = " << min_size_ << endl;

                // cout << "sending (" << i << "," << j << ") to bucket " << bucket << endl;
                MPI_Send(&from, 1, MPI_INT64_T, bucket + 1, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
                MPI_Send(&to, 1, MPI_INT64_T, bucket + 1, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
            }
            // cout << endl;
        }
        // cout << endl;
        if (i % 100 == 0) {
            cout << "partitioning line " << i << endl;
        }
        // cout << endl;
    }

    for (int i = 1; i <= num_partitions_; i++) {
        int64_t terminate = -1;
        MPI_Send(&terminate, 1, MPI_INT64_T, i, FROM_MASTER_PARTITION, MPI_COMM_WORLD);
        // cout << "sending -1 to bucket" << i - 1 << endl;
    }
    // fout.close();
}

/*
void DependencyGraph::HDRF_partition(GraphNL& graph_1, GraphNL& graph_2) {
    cout << "Partitioning using HDRF algorithm..." << endl;

    ifstream fin("dataset/mz_test/" + graph_1.get_graph_name() + "_" + graph_2.get_graph_name() + "_dependency_graph.txt");
    if (!fin) {
        cout << "Cannot open dependency graph!" << endl;
        exit(0);
    }

    ofstream fout("dataset/mz_test/" + graph_1.get_graph_name() + "_" + graph_2.get_graph_name() + "_dependency_graph_" + to_string(num_partitions_) + ".txt");

    int64_t chunk_size;
    int64_t num_edges = dependency_count_;

    max_size_ = 0;
    min_size_ = 0;

    // TODO: optimization - chunk read
    string line, data;
    uint32_t idx;
    uint32_t from, to;
    uint32_t bucket;
    while (getline(fin, line)) {
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> from;

        line.erase(0, idx + 1);
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> to;

        bucket = HDRF_calculate_partition(from, to);

        occupied_[bucket]++;
        is_boundarys_[bucket].set(from);
        is_boundarys_[bucket].set(to);

        if (occupied_[bucket] > max_size_) {
            max_size_ = occupied_[bucket];
        }

        int min_size_bucket_count = 0;
        for (int i = 0; i < num_partitions_; i++) {
            if (occupied_[i] == min_size_) {
                min_size_bucket_count++;
            }
        }
        if (min_size_bucket_count == 0) {
            min_size_++;
        }
        // TODO: error?
        //if (occupied_[bucket] == min_size_) {
        //    int min_sized_bucket_count = 0;
        //    for (int i = 0; i < num_partitions_; i++) {
        //        if (occupied_[i] == min_size_) {
        //            min_sized_bucket_count++;
        //        }
        //    }
        //    if (min_sized_bucket_count == 1) {
        //        min_size_++;
        //    }
        //}

        fout << from << " " << to << " " << bucket << endl;
    }

    fout.close();
    fin.close();

    cout << "Partitioning complete." << endl;
}
*/

uint32_t DependencyGraph::HDRF_calculate_partition(uint32_t from, uint32_t to) {
    double best_score = -1.0;
    uint32_t best_partition = 0;

    for (uint32_t i = 0; i < num_partitions_; i++) {
        double score = compute_partition_score(from, to, i);
        // cout << "score for partition " << i << " is " << score << endl;
        if (score > best_score) {
            best_score = score;
            best_partition = i;
        }
    }
    return best_partition;
}

double DependencyGraph::compute_partition_score(uint32_t u, uint32_t v, uint32_t bucket_id) {
    if (occupied_[bucket_id] >= capacity_) {
        // cout << "partition " << bucket_id << " is full with " << occupied[bucket_id] << endl;
        return -1.0;  // partition is full, do not choose it
    }
    size_t degree_u = count_[u];
    size_t degree_v = count_[v];
    size_t sum = degree_u + degree_v;
    double gu = 0.0, gv = 0.0;
    if (is_boundarys_[bucket_id][u]) {
        gu = degree_u;
        gu /= sum;
        gu = 1 + (1 - gu);
    }
    if (is_boundarys_[bucket_id][v]) {
        gv = degree_v;
        gv /= sum;
        gv = 1 + (1 - gv);
    }

    double bal = (max_size_ - occupied_[bucket_id]) / (1.0 + max_size_ - min_size_);

    double score = gu + gv + lambda * bal;
    return score;
}

void DependencyGraph::count_replication() {
    uint32_t total = 0;
    int i = 0;
    for (auto bs : is_boundarys_) {
        total += bs.count();
        cout << "Partition " << i << " has " << bs.count() << " vertices." << endl;
        i++;
    }
    cout << "Total vertices = " << total << endl;
}

// read edge-partitioned edge list, trasnform into vertex partition
void DependencyGraph::get_vertex_partition(string graph_1_name, string graph_2_name) {
    ifstream fin("dataset/mz_test/" + graph_1_name + "_" + graph_2_name + "_dependency_graph.txt_convert.txt");
    if (!fin) {
        cout << "Cannot open convert file!" << endl;
        exit(0);
    }

    unordered_map<uint32_t, uint32_t> vertex_conversion;
    string line, data;
    uint32_t before, after;
    uint32_t idx;

    while (getline(fin, line)) {
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> before;

        line.erase(0, idx + 1);
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> after;

        vertex_conversion[before] = after;
    }
    fin.close();

    fin.open("dataset/mz_test/" + graph_1_name + "_" + graph_2_name + "_dependency_graph.txt.edgepart." + to_string(num_partitions_));
    if (!fin) {
        cout << "Cannot open partitioned dependency graph!" << endl;
        exit(0);
    }

    vector<vector<uint32_t>> vertex_partition_count(node_pairs_count_, vector<uint32_t>(num_partitions_, 0));
    uint32_t node_1, node_2, partition;
    while (getline(fin, line)) {
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> node_1;

        line.erase(0, idx + 1);
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> node_2;

        line.erase(0, idx + 1);
        idx = line.find(" ");
        data = line.substr(0, idx);
        istringstream(data) >> partition;

        vertex_partition_count[vertex_conversion[node_1]][partition]++;
        vertex_partition_count[vertex_conversion[node_2]][partition]++;
    }
    fin.close();

    cout << "Vertex partition: " << endl;
    for (int i = 0; i < vertex_partition_count.size(); i++) {
        cout << "vertex " << i << ": ";
        for (int j = 0; j < vertex_partition_count[i].size(); j++) {
            cout << vertex_partition_count[i][j] << " ";
        }
        cout << endl;
    }

    vector<uint32_t> vertex_partition(node_pairs_count_, -1);
    for (int i = 0; i < vertex_partition_count.size(); i++) {
        int max = 0;
        int max_idx = -1;
        for (int j = 0; j < vertex_partition_count[i].size(); j++) {
            if (vertex_partition_count[i][j] > max) {
                max = vertex_partition_count[i][j];
                max_idx = j;
            }
        }
    }
}

void DependencyGraph::print_label_mat() {
    initializer_.print_label_info();
    initializer_.print_label_sim();
}

uint64_t DependencyGraph::worker_receive_partition(GraphNL& graph_1, GraphNL& graph_2) {
    MPI_Status status;
    int64_t count = 0;

    sim_values_.resize(graph_1.get_num_vertices());

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

        if (sim_values_[from_1].find(from_2) == sim_values_[from_1].end()) {
            sim_values_[from_1].insert({from_2, 0.0});
            count++;
        }
        if (sim_values_[to_1].find(to_2) == sim_values_[to_1].end()) {
            sim_values_[to_1].insert({to_2, 0.0});
            count++;
        }
        /*
        if (task_id == 4) {
            printf("Worker processor %s, rank %d out of %d processors, received from master (%d,%d)\n", processor_name, task_id, num_tasks, from, to);
        }
        */
    }

    /*
    if (task_id == 1) {
        cout << sim_values_.size() << endl;
        for (int i = 0; i < sim_values_.size(); i++) {
            cout << i << ": ";
            for (auto v : sim_values_[i]) {
                cout << "(" << v.first << "," << v.second << "), ";
            }
            cout << endl;
        }
    }
    */
    return count;
}

void DependencyGraph::worker_initialize(GraphNL& graph_1, GraphNL& graph_2) {
    for (uint32_t i = 0; i < sim_values_.size(); i++) {
        for (auto v : sim_values_[i]) {
            sim_values_[i][v.first] = initializer_.label_sim_matrix_[node_label_1_[i]][node_label_2_[v.first]];
        }
    }
    /*
    if (task_id_ == 1) {
        cout << sim_values_.size() << endl;
        for (int i = 0; i < sim_values_.size(); i++) {
            cout << i << ": ";
            for (auto v : sim_values_[i]) {
                cout << setprecision(3) << "(" << v.first << "," << v.second << "), ";
            }
            cout << endl;
        }
        //print_label_mat();
    }
    */
}

void DependencyGraph::distribute_partition_record() {
    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers_;
    uint64_t extra = node_pairs_count_ % num_workers_;
    // cout << node_pairs_per_part << endl;
    /*
    cout << "bit set = " << endl;
    for (int i = 0; i < is_boundarys_.size(); i++) {
        cout << is_boundarys_[i] << endl;
    }
    cout << endl;
    */

    // for (boost::dynamic_bitset<>::size_type i = 0; i < is_boundarys_[0].size(); ++i) {
    uint64_t i = 0;
    uint64_t idx = 0;
    for (int p = 1; p <= num_workers_; p++) {
        idx += (p <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part;
        // cout << "Partition " << p << " assigned with " << ((p <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part) << " values, has values: " << endl;
        for (; i < idx; i++) {
            // cout << "node_pair " << i << " = ";
            uint16_t tmp = 0;
            for (int j = 0; j < is_boundarys_.size(); j++) {
                tmp = (tmp << 1) | is_boundarys_[j][i];
                //cout << j << "," << i << "=" << is_boundarys_[j][i] << " ";
            }
            // cout << tmp << endl;
            MPI_Send(&tmp, 1, MPI_UINT16_T, p, FROM_MASTER_PART_RECORD, MPI_COMM_WORLD);
        }
        // cout << endl;
    }
}

void DependencyGraph::worker_receive_partition_record(GraphNL& graph_2) {
    MPI_Status status;
    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers_;
    uint64_t extra = node_pairs_count_ % num_workers_;
    uint64_t num_receive = (task_id_ <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part;

    uint64_t offset;
    if (task_id_ <= extra) {
        offset = (node_pairs_per_part + 1) * (task_id_ - 1);
    } else {
        offset = ((node_pairs_per_part + 1) * extra + node_pairs_per_part * (task_id_ - extra - 1));
    }
    uint32_t u_offset = offset / graph_2.get_num_vertices();
    uint32_t u_finish = (offset + num_receive) / graph_2.get_num_vertices();

    //if (task_id) {
    // printf("Worker processor %s, rank %d out of %d processors withh offset %d, starts with index %d, ends with index %d\n", processor_name, task_id, num_tasks, offset, u_offset, u_finish);
    //}

    part_record_.resize(u_finish - u_offset + 1);
    uint32_t current_u = 0;
    for (uint64_t i = offset; i < offset + num_receive; i++) {
        uint16_t tmp;
        MPI_Recv(&tmp, 1, MPI_UINT16_T, MASTER, FROM_MASTER_PART_RECORD, MPI_COMM_WORLD, &status);
        boost::dynamic_bitset<> bs(16, tmp);
        bs.resize(num_partitions_);

        uint32_t u = i / graph_2.get_num_vertices();
        uint32_t v = i % graph_2.get_num_vertices();

        part_record_[u - u_offset].insert({v, bs});
    }
    /*
    if (task_id == 2)
        printf("Worker processor %s, rank %d out of %d processors, partition record = \n", processor_name, task_id, num_tasks, offset);
    for (int i = 0; i < part_record_.size(); i++) {
        for (auto v : part_record_[i]) {
            //if (task_id == 2)
            //    cout << "(" << i + u_offset << "," << v.first << ") = " << v.second << endl;
        }
    }
    */
}

void DependencyGraph::worker_copy_values() {
    pre_values_.assign(sim_values_.begin(), sim_values_.end());
    /*
    if (task_id_ == 1) {
        cout << pre_values_.size() << endl;
        for (int i = 0; i < pre_values_.size(); i++) {
            cout << i << ": ";
            for (auto v : pre_values_[i]) {
                cout << setprecision(3) << "(" << v.first << "," << v.second << "), ";
            }
            cout << endl;
        }
    }
    */
}

uint64_t DependencyGraph::worker_random_partition(GraphNL& graph_1, GraphNL& graph_2) {
    int64_t count = 0;

    uint32_t avg_row, rows, extra;
    avg_row = graph_1.get_num_vertices() / num_workers_;
    extra = graph_1.get_num_vertices() % num_workers_;
    rows = (task_id_ <= extra) ? avg_row + 1 : avg_row;
    sim_values_.resize(graph_1.get_num_vertices());

    uint32_t offset;
    if (task_id_ <= extra) {
        offset = (rows + 1) * (task_id_ - 1);
    } else {
        offset = ((rows + 1) * extra + rows * (task_id_ - extra - 1));
    }

    for (uint32_t i = offset; i < offset+rows; i++) {
        for (uint32_t j = 0; j < graph_2.get_num_vertices(); j++) {
            if (sim_values_[i].find(j) == sim_values_[i].end()) {
                sim_values_[i].insert({j, 0.0});
                count++;
            }
        }
    }
    /*
    if (task_id_ == 4) {
        cout << sim_values_.size() << endl;
        for (int i = 0; i < sim_values_.size(); i++) {
            cout << i << ": ";
            for (auto v : sim_values_[i]) {
                cout << "(" << v.first << "," << v.second << "), ";
            }
            cout << endl;
        }
    }
    */
    return count;
}

void DependencyGraph::worker_random_compute(GraphNL& graph_1, GraphNL& graph_2, float w_i, float w_o, float w_l, string label_constrainted){
    MPI_Status status;
    uint32_t avg_row, rows, extra;
    avg_row = graph_1.get_num_vertices() / num_workers_;
    extra = graph_1.get_num_vertices() % num_workers_;
    rows = (task_id_ <= extra) ? avg_row + 1 : avg_row;

    uint32_t u_offset;
    if (task_id_ <= extra) {
        u_offset = (avg_row + 1) * (task_id_ - 1);
    } else {
        u_offset = ((avg_row + 1) * extra + avg_row * (task_id_ - extra - 1));
    }

    vector<uint32_t> in_degree_vec1 = graph_1.get_in_degree_vec();
    vector<uint32_t> out_degree_vec1 = graph_1.get_out_degree_vec();
    vector<uint32_t> in_degree_vec2 = graph_2.get_in_degree_vec();
    vector<uint32_t> out_degree_vec2 = graph_2.get_out_degree_vec();

    for (uint32_t i = 0; i < sim_values_.size(); i++) {
        NeighborsMap map_i_out = graph_1.get_out_neighbors_map(i);
        NeighborsMap map_i_in = graph_1.get_in_neighbors_map(i);
        /*
        if (task_id_ == 2) {
            cout << "node " << i << ": {";
            for (auto v : map_i_in) {
                cout << v.first << ": [";
                for (auto u : v.second) {
                    cout << u << ",";
                }
                cout << "],";
            }
            cout << "}" << endl;
        }
        */
        if (sim_values_[i].empty()) continue;

        for (auto j : sim_values_[i]) {
            float label_sim = 0.0, in_neighbor_sim = 0.0, out_neighbor_sim = 0.0;
            NeighborsMap map_j_out = graph_2.get_out_neighbors_map(j.first);
            NeighborsMap map_j_in = graph_2.get_in_neighbors_map(j.first);

            vector<uint32_t> neighbors_out_u, neighbors_out_v, neighbors_in_u, neighbors_in_v;
            if (w_i > 0 && map_i_in.size() > 0 && map_j_in.size() > 0) {
                printf("i = %d, j = %d: \n", i, j);
                in_neighbor_sim = match_score(map_i_in, map_j_in, avg_row, extra, label_constrainted) / denom_func(in_degree_vec1[i], in_degree_vec2[j.first]);
                //if (task_id_ == 1) printf("numerator = %.3f\n", match_score(map_i_in, map_j_in, avg_row, extra, label_constrainted));
                //if (task_id_ == 1) printf("in_neighbor_sim = %.3f\n", in_neighbor_sim);
            }

            if (map_i_out.size() > 0 && map_j_out.size() > 0) {
                out_neighbor_sim = match_score(map_i_out, map_j_out, avg_row, extra, label_constrainted) / denom_func(out_degree_vec1[i], out_degree_vec2[j.first]);
                //if (task_id_ == 1) printf("out_neighbor_sim = %.3f\n", out_neighbor_sim);
            }
            label_sim = initializer_.label_sim_matrix_[node_label_1_[i]][node_label_2_[j.first]];
            //if (task_id_ == 1) printf("label_sim = %.3f\n", label_sim);

            if (map_i_in.size() == 0 && map_j_in.size() == 0) {
                in_neighbor_sim = 1.0;
            }
            if (map_i_out.size() == 0 && map_j_out.size() == 0) {
                out_neighbor_sim = 1.0;
            }
            sim_values_[i][j.first] = out_neighbor_sim * w_o + in_neighbor_sim * w_i + label_sim * w_l;
            //if (task_id_ == 1) printf("sim_value = %.3f\n\n", sim_values_[i][j.first]);

        }
    }
}

float DependencyGraph::match_score(NeighborsMap& map_i, NeighborsMap& map_j, uint32_t avg_row, uint32_t extra, string label_constrainted) {
    vector<uint32_t> sub_i, sub_j, neighbors_i, neighbors_j;
    NeighborsMap::iterator iter, jter;
    if (map_i.size() == 0 && map_j.size() == 0) {
        return 1.0;
    } else if (map_i.size() == 0 || map_j.size() == 0) {
        return 0.0;
    }

    float match_score = 0.0;
    if (label_constrainted == "Y") {
        for (iter = map_i.begin(); iter != map_i.end(); iter++) {
            sub_i = iter->second;
            if (map_j.find(iter->first) != map_j.end()) {
                sub_j = map_j[iter->first];
                match_score += match_score_without_label_constraint(sub_i, sub_j, avg_row, extra);
            }
        }
    }
    else {
        for (iter = map_i.begin(); iter != map_i.end(); iter++) {
            neighbors_i.insert(neighbors_i.end(), iter->second.begin(), iter->second.end());
        }
        for (jter = map_j.begin(); jter != map_j.end(); jter++) {
            neighbors_j.insert(neighbors_j.end(), jter->second.begin(), jter->second.end());
        }/*
        if (task_id_ == 1){
            cout << "neighbors_i = [";
            for (auto v : neighbors_i) {
                cout << v << ", ";
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
            cout << "]\n";

            cout << "neighbors_j = [";
            for (auto v : neighbors_j) {
                cout << v << ", ";
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
            cout << "]\n";
        }*/
        if (partition_method_ == "random") {
            match_score = match_score_without_label_constraint(neighbors_i, neighbors_j, avg_row, extra);
        }
        else {
            match_score = match_score_without_label_constraint_hdrf(neighbors_i, neighbors_j, avg_row, extra);
        }
        //if (task_id_ == 1) printf("match_score = %.3f\n", match_score);
    }
    return match_score;
}

float DependencyGraph::match_score_without_label_constraint(Neighbors& N_i, Neighbors& N_j, uint32_t avg_row, uint32_t extra) {
    MPI_Status t0_status;
    // Inject match
    if (simulation_type_ == "sim") {
        float weight;
        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        }
        else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        //if (task_id_ == 2) cout << "hello\n";
        vector<uint32_t>::iterator iItr, jItr;
        float numerator = 0.0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            float max_score = 0.0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                //if (task_id_ == 1) cout << "neighbor pair = (" << *iItr << "," << *jItr << ") ";
                if (sim_values_[*iItr].empty()){
                    //if (task_id_ == 1) cout << "NEED TO FETCH!  ";
                    // if (task_id_ == 2) cout << "Yes, it's empty" << endl;
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    //if (task_id_ == 2) cout << "bucket = " << bucket << endl;
                    uint32_t i_value = *iItr;
                    uint32_t np[2] = {*iItr, *jItr};
                    //if (task_id_ == 2) cout << "np = {" << np[0] << "," << np[1] << "}\n";
                    //MPI_Send(&i_value, 1, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                    //if (task_id_ == 1) cout << "weight = " << weight << endl;
                    //weight = 0.0;
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                    //if (task_id_ == 1) cout << "weight = " << weight << endl;
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }
        //exit(0);
        //if (task_id_ == 2) cout << "hello" << endl;
        //if (task_id_ == 2) printf("numerator = %.3f\n", numerator);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1));
        return numerator;
    }
    // BiMatch
    else if (simulation_type_ == "bisim") {
        float weight;
        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        } else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        vector<uint32_t>::iterator iItr, jItr;
        float numerator = 0.0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            float max_score = 0.0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }

        for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
            float max_score = 0.0;
            for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }
        return numerator;
    }
    // BijectMatch, if simulatin = "dpsim" or "bijectsim"
    else {
        vector<float> Weight;
        float weight;
        priority_queue<MatchItem, vector<MatchItem>, WeightLessThan<MatchItem>> matchQueue;

        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        } else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        vector<uint32_t>::iterator iItr, jItr;
        MatchItem matchitem;
        uint32_t x = 0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            uint32_t y = 0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > 0.0) {
                    matchitem.x = x;
                    matchitem.y = y;
                    matchitem.weight = weight;
                    matchQueue.push(matchitem);
                }
                y += 1;
            }
            x += 1;
        }
        vector<bool> Flag1 = vector<bool>(size_i, false);
        vector<bool> Flag2 = vector<bool>(size_j, false);
        float numerator = 0.0;
        uint32_t matchSize = min(size_i, size_j);
        uint32_t matchCount = 0;
        uint32_t x_pos, y_pos;
        while (!matchQueue.empty() && matchCount < matchSize) {
            matchitem = matchQueue.top();
            x_pos = matchitem.x;
            y_pos = matchitem.y;
            if (Flag1[x_pos] == false && Flag2[y_pos] == false) {
                Flag1[x_pos] = true;
                Flag2[y_pos] = true;
                numerator += matchitem.weight;
                matchCount++;
            }
            matchQueue.pop();
        }
        return numerator;
    }
}

float DependencyGraph::denom_func(uint32_t i, uint32_t j) {
    if (simulation_type_ == "sim" || simulation_type_ == "dpsim") { // query side
        return float(i);
    }
    else if (simulation_type_ == "bisim") {   // add
        return float(i+j);
    }
    else {  // root, simulation_type_ == "bijectsim"
        return float((sqrt(i) * sqrt(j)));
    }
    // max
    // return float(i > j ? i : j)
    //
    // halfAdd
    // return float((i+j) / 2.0)
}

void DependencyGraph::worker_print() {
    for (int i = 0; i < sim_values_.size(); i++) {
        if (sim_values_[i].empty()) continue;
        for (auto v : sim_values_[i]) {
            cout << setw(2) << setprecision(3) << "(" << i << "," << v.first << ")" << "=" << v.second;
        }
        cout << endl << endl;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void DependencyGraph::worker_listen() {
    MPI_Status status;
    while (true) {
        float value = -1;
        // uint32_t i_value = 0;
        uint32_t node_pair[2] = {0,0};

        // printf("Listenning...\n");
        MPI_Recv(&node_pair[0], 2, MPI_UINT32_T, MPI_ANY_SOURCE, QUERY_VALUE, MPI_COMM_WORLD, &status);
        if (status.MPI_SOURCE == MASTER) break;

        //if (task_id == 1) printf("Received %d", i_value);//printf("Received {%d,%d}\n",np[0],np[1]);
        value = pre_values_[node_pair[0]][node_pair[1]];
        MPI_Send(&value, 1, MPI_FLOAT, status.MPI_SOURCE, ANSWER_VALUE, MPI_COMM_WORLD);
    }
}

float DependencyGraph::match_score_without_label_constraint_hdrf(Neighbors& N_i, Neighbors& N_j, uint32_t avg_row, uint32_t extra) {
    MPI_Status t0_status;
    // Inject match
    if (simulation_type_ == "sim") {
        float weight;
        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        }
        else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        //if (task_id_ == 2) cout << "hello\n";
        vector<uint32_t>::iterator iItr, jItr;
        float numerator = 0.0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            float max_score = 0.0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                // if don't have value
                if (sim_values_[*iItr].empty()){
                    //if (task_id_ == 1) cout << "NEED TO FETCH!  ";
                    // if (task_id_ == 2) cout << "Yes, it's empty" << endl;
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    //if (task_id_ == 2) cout << "bucket = " << bucket << endl;
                    uint32_t i_value = *iItr;
                    uint32_t np[2] = {*iItr, *jItr};
                    //if (task_id_ == 2) cout << "np = {" << np[0] << "," << np[1] << "}\n";
                    //MPI_Send(&i_value, 1, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                    //if (task_id_ == 1) cout << "weight = " << weight << endl;
                    //weight = 0.0;
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                    //if (task_id_ == 1) cout << "weight = " << weight << endl;
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }
        //exit(0);
        //if (task_id_ == 2) cout << "hello" << endl;
        //if (task_id_ == 2) printf("numerator = %.3f\n", numerator);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1));
        return numerator;
    }
    // BiMatch
    else if (simulation_type_ == "bisim") {
        float weight;
        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        } else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        vector<uint32_t>::iterator iItr, jItr;
        float numerator = 0.0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            float max_score = 0.0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }

        for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
            float max_score = 0.0;
            for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > max_score) {
                    max_score = weight;
                }
            }
            numerator += max_score;
        }
        return numerator;
    }
    // BijectMatch, if simulatin = "dpsim" or "bijectsim"
    else {
        vector<float> Weight;
        float weight;
        priority_queue<MatchItem, vector<MatchItem>, WeightLessThan<MatchItem>> matchQueue;

        uint32_t size_i = N_i.size();
        uint32_t size_j = N_j.size();
        if (size_i == 0 && size_j == 0) {
            return 1.0;
        } else if (size_i == 0 || size_j == 0) {
            return 0.0;
        }
        vector<uint32_t>::iterator iItr, jItr;
        MatchItem matchitem;
        uint32_t x = 0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            uint32_t y = 0;
            for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
                if (sim_values_[*iItr].empty()){
                    int bucket = 0;
                    if (*iItr < extra * (avg_row+1)) {
                        bucket = *iItr / (avg_row+1) + 1;
                    } else {
                        bucket = extra + (*iItr - extra*(avg_row+1)) /avg_row + 1;
                    }
                    uint32_t np[2] = {*iItr, *jItr};
                    MPI_Send(&np[0], 2, MPI_UINT32_T, bucket, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&weight, 1, MPI_FLOAT, bucket, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                }
                else {
                    weight = pre_values_[*iItr][*jItr];
                }
                if (weight > 0.0) {
                    matchitem.x = x;
                    matchitem.y = y;
                    matchitem.weight = weight;
                    matchQueue.push(matchitem);
                }
                y += 1;
            }
            x += 1;
        }
        vector<bool> Flag1 = vector<bool>(size_i, false);
        vector<bool> Flag2 = vector<bool>(size_j, false);
        float numerator = 0.0;
        uint32_t matchSize = min(size_i, size_j);
        uint32_t matchCount = 0;
        uint32_t x_pos, y_pos;
        while (!matchQueue.empty() && matchCount < matchSize) {
            matchitem = matchQueue.top();
            x_pos = matchitem.x;
            y_pos = matchitem.y;
            if (Flag1[x_pos] == false && Flag2[y_pos] == false) {
                Flag1[x_pos] = true;
                Flag2[y_pos] = true;
                numerator += matchitem.weight;
                matchCount++;
            }
            matchQueue.pop();
        }
        return numerator;
    }
}

void DependencyGraph::write_result_txt() {
    ofstream fout("result_partition" + to_string(task_id_) + ".txt");
    for (int i = 0; i < sim_values_.size(); i++) {
        for (auto v : sim_values_[i]) {
            fout << setprecision(3) << i << " " << v.first << " " << v.second << endl;
        }
    }
    fout.close();
}