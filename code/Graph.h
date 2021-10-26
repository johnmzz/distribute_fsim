#ifndef _GRAPH_H
#define _GRAPH_H

#include <assert.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

typedef uint32_t Vertex;
typedef vector<Vertex> Vertices;
typedef vector<Vertex> Neighbors;
typedef unordered_map<uint32_t, Neighbors> NeighborsMap;
typedef vector<NeighborsMap> AdjMap;
typedef vector<vector<uint32_t>> AdjList;
typedef vector<uint32_t> VertexType;

class Graph {
  protected:
    string graph_name_;
    uint32_t vsize_;                   // number of vertices
    uint32_t esize_;                   // number of edges
    uint32_t lsize_;                   // number of labels
    AdjMap graph_out_;                 // eg. see below
    AdjMap graph_in_;                  // ^
    AdjList list_in_;                  // eg. [[v0_n0, v0_n1, ...], [v1_n0, v1_n1, ...], ...]
    AdjList list_out_;                 // ^
    VertexType vertex_type_;           // eg. [v0_label, v1_label, v2_label, ...]
    vector<string> label_info_;        // eg. [label_0, label_1, label_2, ...]
    vector<uint32_t> in_degree_vec_;   // eg. [v0_in_deg, v1_in_deg, v2_in_deg, ...]
    vector<uint32_t> out_degree_vec_;  // eg. [v0_out_deg, v1_out_deg, v2_out_deg, ...]
    vector<uint32_t> label_degree_;    // eg. [label_0_num_v, label_1_num_v, label_2_num_v, ...]

  public:
    Graph();
    Graph(string path, string graph);
    void read_vertex_type(string dataset_path, string graph_name);
    void read_label_info(string dataset_path, string graph_name);
    virtual void read_edgelist(ifstream& input_file) = 0;

    // accessor methods
    uint32_t get_num_vertices();
    uint32_t get_num_labels();
    string get_graph_name();
    VertexType& get_vertex_type();
    vector<string>& get_label_info();
    vector<uint32_t>& get_in_degree_vec();
    vector<uint32_t>& get_out_degree_vec();
    vector<uint32_t>& get_label_degree();
    NeighborsMap& get_out_neighbors_map(uint32_t node);
    NeighborsMap& get_in_neighbors_map(uint32_t node);
    vector<uint32_t>& get_out_neighbors_list(uint32_t node);
    vector<uint32_t>& get_in_neighbors_list(uint32_t node);

    // printing methods
    void print_label_info();
    void print_label_degree();
    void print_vertex_type();
    void print_in_degree_vector();
    void print_out_degree_vector();
    void print_list_in();
    void print_list_out();
    void print_graph_out();
    void print_graph_in();
    void get_graph_statistics();
    void print_graph();
};

// Node Labeled Graph
class GraphNL : public Graph {
  public:
    GraphNL() : Graph(){};
    GraphNL(string, string);
    void read_edgelist(ifstream& input_file);
};
#endif

// eg. AdjMap
// graph_out_ = [
//     {1:[3], 0:[1]},                          // node 0 has neighbor 3 with label 1, neighbor 1 with label 0
//     {1:[8], 0:[2]},
//     {1:[3,9], 0:[0]},
//     {0:[1]},
//     {7:[7], 4:[6], 3:[5]},
//     {1:[8], 7:[7]},
//     {6:[11], 1:[8], 3:[5]},
//     {1:[8 ], 4:[6]},
//     {5:[10], 2:[4], 1:[3], 0:[0,2]},
//     {9:[12], 5:[10], 1:[8]},
//     {},
//     {5:[10]},
//     {}
// ]