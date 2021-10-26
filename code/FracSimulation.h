#ifndef _FRACSIMULATION_H
#define _FRACSIMULATION_H

#include <assert.h>
#include <omp.h>
#include <time.h>

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "Denominator.h"
#include "Graph.h"
#include "Initializer.h"
#include "MatchMethod.h"

using namespace std;
using namespace chrono;

typedef vector<vector<float>> SimMatrix;
typedef vector<float> VectorFloatType;
typedef vector<uint32_t> VertexType;

class FracSimulation {
  private:
    GraphNL graph_1_, graph_2_;
    unique_ptr<Initializer> ini_ptr_;
    unique_ptr<MatchMethod> mm_ptr_;
    unique_ptr<Denominator> denom_ptr_;
    SimMatrix simmat_;
    SimMatrix labelmat_;
    uint32_t n1_, n2_, l1_, l2_;
    VertexType node_label_1_, node_label_2_;  // eg. [v0_type, v1_type, v2_type, ...]

  public:
    FracSimulation(GraphNL& graph_1, GraphNL& graph_2, unique_ptr<Initializer> ini_ptr, unique_ptr<MatchMethod> mm_ptr, unique_ptr<Denominator> denom_ptr);

    // basic operations
    SimMatrix& get_curr_simmat();
    void simmat_copy_to(SimMatrix& premat);
    float get_value(uint32_t i, uint32_t j);
    void set_value(uint32_t i, uint32_t j, float sim);
    void simmat_clear();
    void print_simmat();
    void print_labelmat();

    // initialization
    void initialize(uint32_t num_thus);

    // iteration
    bool diff(SimMatrix& premat, float threshold, string mode);
    bool diff_relative(SimMatrix& premat, float threshold);
    bool diff_maxdiff(SimMatrix& premat, float threshold);
    bool diff_totaldiff(SimMatrix& premat, float threshold);

    void update_simmatrix(SimMatrix& premat, float w_i, float w_o, float w_l, uint32_t num_thus);
};

#endif
