#include "DependencyGraph.h"
DependencyGraph::DependencyGraph(GraphNL& graph_1, GraphNL& graph_2, uint32_t num_partitions, string initializer_type, string simulation_type)
    : initializer_(initializer_type, graph_1.get_label_info(), graph_2.get_label_info()) {
    lambda = 1.1;  // default value for HDRF
    capacity_ = 0;
    dependency_count_ = 0;
    node_pairs_count_ = graph_1.get_num_vertices() * graph_2.get_num_vertices();
    num_partitions_ = num_partitions;
    simulation_type_ = simulation_type;
    occupied_.assign(num_partitions, 0);
    node_label_1_ = graph_1.get_vertex_type();
    node_label_2_ = graph_2.get_vertex_type();
}

int64_t DependencyGraph::encode(uint32_t i, uint32_t j, uint32_t offset) {
    return i * offset + j;
}

void DependencyGraph::init_dependency_graph(GraphNL& graph_1, GraphNL& graph_2) {
    cout << endl
         << "Initializing dependency graph..." << endl;
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
    cout << "Capacity of each partition = " << capacity_ << endl;
    // cout << "Degree of each node-pair = ";
    // for (auto c : count_) {
    //    cout << c << ", ";
    //}
    //cout << endl;
    cout << "Dependency graph initialization completed." << endl;
}

void DependencyGraph::HDRF_partition(GraphNL& graph_1, GraphNL& graph_2) {
    cout << endl
         << "Partitioning using HDRF algorithm..." << endl;

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
    cout << "Partitioning complete." << endl;
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

void DependencyGraph::worker_receive_partition(GraphNL& graph_1, GraphNL& graph_2) {
    // printf("Worker processor %s, rank %d out of %d processors\n", processor_name, task_id, num_tasks);
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;
    int num_workers = num_tasks - 1;
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

        if (task_id == 4) {
            // printf("Worker processor %s, rank %d out of %d processors, received from master (%d,%d)\n", processor_name, task_id, num_tasks, from, to);
        }
    }
    if (task_id == 4) {
        // printf("Worker processor %s, rank %d out of %d processors finished receiving %d node-pairs.\n", processor_name, task_id, num_tasks, count);
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
    MPI_Send(&count, 1, MPI_INT64_T, MASTER, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD);
}

void DependencyGraph::worker_initialize(GraphNL& graph_1, GraphNL& graph_2) {
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;

    int start_initialize;
    MPI_Recv(&start_initialize, 1, MPI_INT, MASTER, FROM_MASTER_INITIALIZE, MPI_COMM_WORLD, &status);
    // printf("Worker processor %s, rank %d out of %d processors starting to initialize sim values.\n", processor_name, task_id, num_tasks);

    for (uint32_t i = 0; i < sim_values_.size(); i++) {
        for (auto v : sim_values_[i]) {
            sim_values_[i][v.first] = initializer_.label_sim_matrix_[node_label_1_[i]][node_label_2_[v.first]];
        }
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
        print_label_mat();
    }
    */

    int intialization_complete = 1;
    // printf("Worker processor %s, rank %d out of %d processors initialization completed.\n", processor_name, task_id, num_tasks);
    MPI_Send(&intialization_complete, 1, MPI_INT, MASTER, WORKER_INIT_COMPLETE, MPI_COMM_WORLD);
}

void DependencyGraph::distribute_partition_record() {
    cout << endl
         << "Sending workers partition record..." << endl;

    MPI_Status status;
    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    int num_workers = num_tasks - 1;

    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers;
    uint64_t extra = node_pairs_count_ % num_workers;
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
    for (int p = 1; p <= num_workers; p++) {
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
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;
    int num_workers = num_tasks - 1;
    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers;
    uint64_t extra = node_pairs_count_ % num_workers;
    uint64_t num_receive = (task_id <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part;

    uint64_t offset;
    if (task_id <= extra) {
        offset = (node_pairs_per_part + 1) * (task_id - 1);
    } else {
        offset = ((node_pairs_per_part + 1) * extra + node_pairs_per_part * (task_id - extra - 1));
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

    int received_part_record = 1;
    MPI_Send(&received_part_record, 1, MPI_INT, MASTER, WORKER_RECEIVED_PART_RECORD, MPI_COMM_WORLD);
}

void DependencyGraph::worker_copy_values() {
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;

    int start_copy;
    MPI_Recv(&start_copy, 1, MPI_INT, MASTER, FROM_MASTER_COPY, MPI_COMM_WORLD, &status);
    // printf("Worker processor %s, rank %d out of %d processors start copying...\n", processor_name, task_id);

    pre_values_.assign(sim_values_.begin(), sim_values_.end());

    /*
    if (task_id == 1) {
        cout << pre_values_.size() << endl;
        for (int i = 0; i < pre_values_.size(); i++) {
            cout << i << ": ";
            for (auto v : pre_values_[i]) {
                cout << "(" << v.first << "," << v.second << "), ";
            }
            cout << endl;
        }
        print_label_mat();
    }
    */

    int copy_completed = 1;
    MPI_Send(&copy_completed, 1, MPI_INT, MASTER, WORKER_COPY_COMPLETE, MPI_COMM_WORLD);
}

void DependencyGraph::worker_compute(GraphNL& graph_1, GraphNL& graph_2, float w_i, float w_o, float w_l, string label_constrainted) {
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;
    int num_workers = num_tasks - 1;
    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers;
    uint64_t extra = node_pairs_count_ % num_workers;
    uint64_t num_receive = (task_id <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part;

    uint64_t offset;
    if (task_id <= extra) {
        offset = (node_pairs_per_part + 1) * (task_id - 1);
    } else {
        offset = ((node_pairs_per_part + 1) * extra + node_pairs_per_part * (task_id - extra - 1));
    }
    uint32_t u_offset = offset / graph_2.get_num_vertices();

    // *********************************** Local computation round ***********************************
    int start_compute;
    MPI_Recv(&start_compute, 1, MPI_INT, MASTER, FROM_MASTER_COPY, MPI_COMM_WORLD, &status);
    // printf("Worker processor %s, rank %d out of %d processors start computing...\n", processor_name, task_id);

    vector<uint32_t> in_degree_vec1 = graph_1.get_in_degree_vec();
    vector<uint32_t> out_degree_vec1 = graph_1.get_out_degree_vec();
    vector<uint32_t> in_degree_vec2 = graph_2.get_in_degree_vec();
    vector<uint32_t> out_degree_vec2 = graph_2.get_out_degree_vec();

    for (uint32_t i = 0; i < sim_values_.size(); i++) {
        uint32_t u = i + u_offset;
        NeighborsMap map_i_out = graph_1.get_out_neighbors_map(u);
        NeighborsMap map_i_in = graph_1.get_in_neighbors_map(u);
        for (auto v : sim_values_[i]) {
            float label_sim = 0.0, in_neighbor_sim = 0.0, out_neighbor_sim = 0.0;
            NeighborsMap map_j_out = graph_2.get_out_neighbors_map(v.first);
            NeighborsMap map_j_in = graph_2.get_in_neighbors_map(v.first);

            vector<uint32_t> neighbors_out_u, neighbors_out_v, neighbors_in_u, neighbors_in_v;
            if (label_constrainted != "Y") {
                neighbors_out_u = graph_1.get_out_neighbors_list(u);
                neighbors_out_v = graph_2.get_out_neighbors_list(v.first);
                neighbors_in_u = graph_1.get_in_neighbors_list(u);
                neighbors_in_v = graph_2.get_in_neighbors_list(v.first);
            }

            if (neighbors_out_u.size() == 0 && neighbors_out_v.size() == 0) {
                out_neighbor_sim = 1;
            } else if (neighbors_out_u.size() == 0 || neighbors_out_v.size() == 0) {
                out_neighbor_sim = 0;
            } else {
                float numerator = 0.0;
                if (simulation_type_ == "sim") {
                    float weight;
                    vector<uint32_t>::iterator uItr, vItr;
                    for (uItr = neighbors_out_u.begin(); uItr != neighbors_out_u.end(); uItr++) {
                        if (*uItr - u_offset < 0 || *uItr - u_offset >= num_receive) {
                            // node-pair (u,v) not in this partition
                            continue;
                        }
                        float max_score = 0.0;
                        for (vItr = neighbors_out_v.begin(); vItr != neighbors_out_v.end(); vItr++) {
                            if (sim_values_[*uItr - u_offset].find(*vItr) == sim_values_[*uItr - u_offset].end()) {
                                // if node-pair (u,v) not in this partition (v not in u's value list)
                                continue;
                            }
                            weight = pre_values_[*uItr - u_offset][*vItr];
                            if (weight > max_score) {
                                max_score = weight;
                            }
                        }
                        numerator += max_score;
                    }
                } else if (simulation_type_ == "bisim") {
                } else if (simulation_type_ == "dpsim") {
                } else if (simulation_type_ == "bjsim") {
                } else {
                    cout << "Incorrect simulation type!" << endl;
                    exit(0);
                }
                out_neighbor_sim = numerator;
            }
            // sim_values_[u-u_offset][v].first = out_neighbor_sim;

            if (neighbors_in_u.size() == 0 && neighbors_in_v.size() == 0) {
                in_neighbor_sim = 1;
            } else if (neighbors_in_u.size() == 0 || neighbors_out_u.size() == 0) {
                in_neighbor_sim = 0;
            } else {
                float numerator = 0.0;
                if (simulation_type_ == "sim") {
                    float weight;
                    vector<uint32_t>::iterator uItr, vItr;
                    for (uItr = neighbors_in_u.begin(); uItr != neighbors_in_u.end(); uItr++) {
                        if (*uItr - u_offset < 0 || *uItr - u_offset >= num_receive) {
                            // node-pair (u,v) not in this partition
                            continue;
                        }
                        float max_score = 0.0;
                        for (vItr = neighbors_in_v.begin(); vItr != neighbors_in_v.end(); vItr++) {
                            if (sim_values_[*uItr - u_offset].find(*vItr) == sim_values_[*uItr - u_offset].end()) {
                                // if node-pair (u,v) not in this partition (v not in u's value list)
                                continue;
                            }
                            weight = pre_values_[*uItr - u_offset][*vItr];
                            if (weight > max_score) {
                                max_score = weight;
                            }
                        }
                        numerator += max_score;
                    }
                } else if (simulation_type_ == "bisim") {
                } else if (simulation_type_ == "dpsim") {
                } else if (simulation_type_ == "bjsim") {
                } else {
                    cout << "Incorrect simulation type!" << endl;
                    exit(0);
                }
                in_neighbor_sim = numerator;
            }
            // sim_values_[i][v].second = in_neighbor_sim;
        }
    }

    int compute_complete = 1;
    MPI_Send(&compute_complete, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
    // *********************************** Synchronization round ***********************************
    for (int p = 1; p <= num_partitions_; p++) {
        int synch_id;
        MPI_Recv(&synch_id, 1, MPI_INT, MASTER, FROM_MASTER_SYNCHRONIZE, MPI_COMM_WORLD, &status);

        if (synch_id == task_id) {
            // do computation, ask others
        } else {
            // receive queries from others
        }
        MPI_Send(&synch_id, 1, MPI_INT, MASTER, WORKER_SYNCHRONIZE_COMPLETE, MPI_COMM_WORLD);
    }

    // *********************************** Local aggregation round ***********************************
    for (uint32_t i = 0; i < sim_values_.size(); i++) {
        uint32_t u = i + u_offset;
        for (auto v : sim_values_[i]) {
            // sim_values_[i][v.first].first = sim_values_[i][v.first].first / demon(out_deg_vec_1[u], out_deg_vec_2[v.first]);
            // sim_values_[i][v.first].second = sim_values[i][v.first].second / demon(in_deg_vec_1[u], in_deg_vec_2[v.first]);
            // label_sim = labelmat_[node_label_1_[u]][node_label_2_[v.first]];

            // final value
            // sim_values_[i][v.first].first = w_o * sim_values_[i][v.first].first + w_i * sim_values_[i][v.first].second + w_l * label_sim;
        }
    }

    int aggre_complete = 1;
    MPI_Send(&aggre_complete, 1, MPI_INT, MASTER, WORKER_AGGRE_COMPLETE, MPI_COMM_WORLD);
}

void DependencyGraph::worker_random_partition(GraphNL& graph_1, GraphNL& graph_2) {

    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;
    int num_workers = num_tasks - 1;
    int64_t count = 0;

    uint32_t avg_row, rows, extra;
    avg_row = graph_1.get_num_vertices() / num_workers;
    extra = graph_1.get_num_vertices() % num_workers;
    rows = (task_id <= extra) ? avg_row + 1 : avg_row;
    sim_values_.resize(graph_1.get_num_vertices());

    uint32_t offset;
    if (task_id <= extra) {
        offset = (rows + 1) * (task_id - 1);
    } else {
        offset = ((rows + 1) * extra + rows * (task_id - extra - 1));
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
    if (task_id == 4) {
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
    MPI_Send(&count, 1, MPI_INT64_T, MASTER, WORKER_RECEIVED_PARTITION, MPI_COMM_WORLD);
}

void DependencyGraph::worker_random_compute(GraphNL& graph_1, GraphNL& graph_2, float w_i, float w_o, float w_l, string label_constrainted){
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);
    MPI_Status status;
    int num_workers = num_tasks - 1;
    uint64_t node_pairs_per_part = node_pairs_count_ / num_workers;
    uint64_t extra = node_pairs_count_ % num_workers;
    uint64_t num_receive = (task_id <= extra) ? node_pairs_per_part + 1 : node_pairs_per_part;

    uint64_t offset;
    if (task_id <= extra) {
        offset = (node_pairs_per_part + 1) * (task_id - 1);
    } else {
        offset = ((node_pairs_per_part + 1) * extra + node_pairs_per_part * (task_id - extra - 1));
    }
    uint32_t u_offset = offset / graph_2.get_num_vertices();

    int start_compute;
    MPI_Recv(&start_compute, 1, MPI_INT, MASTER, FROM_MASTER_START_COMPUTE, MPI_COMM_WORLD, &status);
    // printf("Worker processor %s, rank %d out of %d processors start computing...\n", processor_name, task_id);

    int num_threads = 0;
    int thread_id = 0;
    #pragma omp parallel default(shared) private(thread_id, num_threads) num_threads(2)
    {
        num_threads = omp_get_num_threads();
        thread_id = omp_get_thread_num();
        int tmp = 0;

        // Calculation thread
        if (thread_id == 0) {
            MPI_Status t0_status;
            vector<uint32_t> in_degree_vec1 = graph_1.get_in_degree_vec();
            vector<uint32_t> out_degree_vec1 = graph_1.get_out_degree_vec();
            vector<uint32_t> in_degree_vec2 = graph_2.get_in_degree_vec();
            vector<uint32_t> out_degree_vec2 = graph_2.get_out_degree_vec();

            for (int worker = 1; worker < num_tasks; worker++) {
                if (worker == task_id) continue;
                for (int j = 0; j < 2; j++) {
                    MPI_Send(&task_id, 1, MPI_INT, worker, QUERY_VALUE, MPI_COMM_WORLD);
                    MPI_Recv(&tmp, 1, MPI_INT, worker, ANSWER_VALUE, MPI_COMM_WORLD, &t0_status);
                    a[worker-1][j] = tmp;
                }
            }

            //printf("Thread %d on worker processor %s, rank %d out of %d processors\n", thread_id, processor_name, task_id, num_tasks);
            int finish = 1;
            MPI_Send(&finish, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
        }

        // Listenner thread
        else {
            MPI_Status t1_status;
            while (true) {
                MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE, QUERY_VALUE, MPI_COMM_WORLD, &t1_status);
                if (t1_status.MPI_SOURCE == MASTER) break;
                MPI_Send(&task_id, 1, MPI_INT, t1_status.MPI_SOURCE, ANSWER_VALUE, MPI_COMM_WORLD);
            }
            //printf("Thread %d on worker processor %s, rank %d out of %d processors, finished listenning.\n", thread_id, processor_name, task_id, num_tasks);
            int finish = 1;
            MPI_Send(&finish, 1, MPI_INT, MASTER, WORKER_COMPUTE_COMPLETE, MPI_COMM_WORLD);
        }
    }
}