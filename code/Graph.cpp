#include "Graph.h"

//*************************************************************************************
//*********************************** Base Class **************************************
// Initialization methods
Graph::Graph() {
    graph_out_ = AdjMap();
    graph_in_ = AdjMap();
    list_in_ = AdjList();
    list_out_ = AdjList();

    vertex_type_ = VertexType();
    label_info_ = vector<string>();
    in_degree_vec_ = vector<uint32_t>();
    out_degree_vec_ = vector<uint32_t>();
    label_degree_ = vector<uint32_t>();
    graph_name_ = string("***");
    vsize_ = 0;
    esize_ = 0;
    lsize_ = 0;
}

Graph::Graph(string dataset_path, string graph_name) {
    graph_name_ = graph_name;
    read_vertex_type(dataset_path, graph_name);
    read_label_info(dataset_path, graph_name);
}

void Graph::read_vertex_type(string dataset_path, string graph_name) {
    ifstream input_file(dataset_path + "/" + graph_name + "/node_label.txt");
    if (!input_file) {
        cout << "Cannot open vertex_label!" << endl;
        exit(0);
    }

    string line, node, label;
    uint32_t node_id, label_id;
    uint32_t idx;
    while (getline(input_file, line)) {
        idx = line.find(" ");
        node = line.substr(0, idx);
        istringstream(node) >> node_id;
        line.erase(0, idx + 1);
        idx = line.find(" ");
        label = line.substr(0, idx);
        istringstream(label) >> label_id;

        if (vertex_type_.size() <= node_id) {
            vertex_type_.resize(node_id + 1, 99999);
        }

        vertex_type_[node_id] = label_id;

        if (label_degree_.size() <= label_id) {
            label_degree_.resize(label_id + 1, 0);
        }
        label_degree_[label_id]++;
    }
}

void Graph::read_label_info(string dataset_path, string graph_name) {
    ifstream input_file(dataset_path + "/" + graph_name + "/label_info.txt");
    if (!input_file) {
        cout << "Cannot open label_info!" << endl;
        exit(0);
    }

    string line, label_id, info;
    uint32_t label;
    uint32_t idx;
    while (getline(input_file, line)) {
        idx = line.find(" ");
        label_id = line.substr(0, idx);
        istringstream(label_id) >> label;
        line.erase(0, idx + 1);
        idx = line.find(" ");
        info = line.substr(0, idx);
        label_info_.push_back(info);
    }
    lsize_ = label_info_.size();
}

// Accessor method
uint32_t Graph::get_num_vertices() {
    assert(vsize_ == graph_out_.size() && vsize_ == graph_in_.size());
    return vsize_;
}

uint32_t Graph::get_num_labels() {
    cout << "in func: " << lsize_ << ", " << label_info_.size() << endl;
    assert(lsize_ == label_info_.size());
    return lsize_;
}

string Graph::get_graph_name() {
    return graph_name_;
}

VertexType& Graph::get_vertex_type() {
    return vertex_type_;
}

vector<string>& Graph::get_label_info() {
    return label_info_;
}

vector<uint32_t>& Graph::get_in_degree_vec() {
    return in_degree_vec_;
}

vector<uint32_t>& Graph::get_out_degree_vec() {
    return out_degree_vec_;
}

vector<uint32_t>& Graph::get_label_degree() {
    return label_degree_;
}

NeighborsMap& Graph::get_out_neighbors_map(uint32_t node) {
    return graph_out_[node];
}

NeighborsMap& Graph::get_in_neighbors_map(uint32_t node) {
    return graph_in_[node];
}

vector<uint32_t>& Graph::get_out_neighbors_list(uint32_t node) {
    return list_out_[node];
}

vector<uint32_t>& Graph::get_in_neighbors_list(uint32_t node) {
    return list_in_[node];
}

// Printing method
void Graph::print_label_info() {
    for (uint32_t i = 0; i < label_info_.size(); i++) {
        cout << i << " : " << label_info_[i] << endl;
    }
}

void Graph::print_label_degree() {
    for (uint32_t i = 0; i < label_degree_.size(); i++) {
        cout << label_degree_[i] << " ";
    }
    cout << endl;
}

void Graph::print_vertex_type() {
    for (uint32_t i = 0; i < vertex_type_.size(); i++) {
        cout << vertex_type_[i] << " ";
    }
    cout << endl;
}

void Graph::print_in_degree_vector() {
    for (uint32_t i = 0; i < in_degree_vec_.size(); i++) {
        cout << in_degree_vec_[i] << " ";
    }
    cout << endl;
}

void Graph::print_out_degree_vector() {
    for (uint32_t i = 0; i < out_degree_vec_.size(); i++) {
        cout << out_degree_vec_[i] << " ";
    }
    cout << endl;
}

void Graph::print_list_in() {
    cout << "Adjacency List (in): " << endl;
    for (uint32_t i = 0; i < list_in_.size(); i++) {
        for (uint32_t j = 0; j < list_in_[i].size(); j++) {
            cout << list_in_[i][j] << " ";
        }
        cout << endl;
    }
}

void Graph::print_list_out() {
    cout << "Adjacency List (out): " << endl;
    for (uint32_t i = 0; i < list_out_.size(); i++) {
        for (uint32_t j = 0; j < list_out_[i].size(); j++) {
            cout << list_out_[i][j] << " ";
        }
        cout << endl;
    }
}

void Graph::print_graph_out() {
    cout << "Adjacency Map (out): " << endl;
    for (int i = 0; i < graph_out_.size(); i++) {
        cout << "Node: " << i << " {";
        for (auto it : graph_out_[i]) {
            cout << "label: " << it.first << ": [";
            for (int j = 0; j < it.second.size(); j++) {
                cout << it.second[j] << " ";
            }
            cout << "], ";
        }
        cout << "}, " << endl;
    }
}

void Graph::print_graph_in() {
    cout << "Adjacency Map (in): " << endl;
    for (int i = 0; i < graph_in_.size(); i++) {
        cout << "Node: " << i << " {";
        for (auto it : graph_in_[i]) {
            cout << "label: " << it.first << ": [";
            for (int j = 0; j < it.second.size(); j++) {
                cout << it.second[j] << " ";
            }
            cout << "], ";
        }
        cout << "}, " << endl;
    }
}

void Graph::get_graph_statistics() {
    cout << "----------------Graph Statistics----------------" << endl;
    cout << "vsize: " << vsize_ << ", esize: " << esize_ << ", lsize: " << lsize_ << endl;

    uint32_t max_degree = 0;
    for (auto it : out_degree_vec_) {
        if (it > max_degree) {
            max_degree = it;
        }
    }
    cout << "Maximum out-degree is " << max_degree << endl;

    uint32_t max_in_degree = 0;
    for (auto it : in_degree_vec_) {
        if (it > max_in_degree) {
            max_in_degree = it;
        }
    }
    cout << "Maximum in-degree is " << max_in_degree << endl;

    unordered_map<uint32_t, vector<uint32_t>> label_node;

    uint32_t count = 0;
    for (auto it : vertex_type_) {
        if (label_node.find(it) == label_node.end()) {
            label_node[it] = {count};
        } else {
            label_node[it].push_back(count);
        }
        count++;
    }

    for (auto it : label_node) {
        cout << it.first << ": " << it.second.size() << endl;
    }
}

void Graph::print_graph() {
    cout << "print graph out" << endl;
    for (uint32_t i = 0; i < graph_out_.size(); i++) {
        cout << i << ": ";
        for (auto it : graph_out_[i]) {
            cout << " <" << it.first << ": ";
            for (auto a : it.second) {
                cout << a << " ";
            }
            cout << ">, ";
        }
        cout << endl;
    }
    cout << "print graph in" << endl;
    for (uint32_t i = 0; i < graph_in_.size(); i++) {
        cout << i << ": ";
        for (auto it : graph_in_[i]) {
            cout << " <" << it.first << ": ";
            for (auto a : it.second) {
                cout << a << " ";
            }
            cout << ">, ";
        }
        cout << endl;
    }
}

//*************************************************************************************
//************************************ Graph NL ***************************************
GraphNL::GraphNL(string dataset_path, string graph_name) {
    graph_name_ = graph_name;
    read_vertex_type(dataset_path, graph_name);
    read_label_info(dataset_path, graph_name);

    ifstream input_file(dataset_path + "/" + graph_name + "/graph.txt");
    if (!input_file) {
        cout << "Cannot open graph edgelist!" << endl;
        cout << dataset_path + "/" + graph_name + "/graph.txt" << endl;
        exit(0);
    }
    read_edgelist(input_file);
}

void GraphNL::read_edgelist(ifstream& input_file) {
    vsize_ = 0;
    esize_ = 0;
    string line, from, to;
    uint32_t from_id, to_id;
    uint32_t idx;

    while (getline(input_file, line)) {
        if (line.empty()) {
            continue;
        }
        idx = line.find(" ");
        from = line.substr(0, idx);
        istringstream(from) >> from_id;
        line.erase(0, idx + 1);
        idx = line.find(" ");
        to = line.substr(0, idx);
        istringstream(to) >> to_id;

        if (from_id >= to_id && from_id >= graph_out_.size()) {
            graph_out_.resize(from_id + 1);
            graph_in_.resize(from_id + 1);
            list_out_.resize(from_id + 1);
            list_in_.resize(from_id + 1);
            in_degree_vec_.resize(from_id + 1, 0);
            out_degree_vec_.resize(from_id + 1, 0);
        } else if (to_id > from_id && to_id >= graph_out_.size()) {
            graph_out_.resize(to_id + 1);
            graph_in_.resize(to_id + 1);
            list_out_.resize(to_id + 1);
            list_in_.resize(to_id + 1);
            in_degree_vec_.resize(to_id + 1, 0);
            out_degree_vec_.resize(to_id + 1, 0);
        }

        if (graph_out_[from_id].find(vertex_type_[to_id]) == graph_out_[from_id].end()) {
            graph_out_[from_id][vertex_type_[to_id]] = {to_id};
        } else {
            graph_out_[from_id][vertex_type_[to_id]].push_back(to_id);
        }

        if (graph_in_[to_id].find(vertex_type_[from_id]) == graph_in_[to_id].end()) {
            graph_in_[to_id][vertex_type_[from_id]] = {from_id};
        } else {
            graph_in_[to_id][vertex_type_[from_id]].push_back(from_id);
        }

        list_out_[from_id].push_back(to_id);
        list_in_[to_id].push_back(from_id);
        in_degree_vec_[to_id]++;
        out_degree_vec_[from_id]++;
        esize_++;
    }
    vsize_ = vertex_type_.size();
}