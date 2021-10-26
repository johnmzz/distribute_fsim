#include "FracSimulation.h"

FracSimulation::FracSimulation(GraphNL& graph_1, GraphNL& graph_2, unique_ptr<Initializer> ini_ptr, unique_ptr<MatchMethod> mm_ptr, unique_ptr<Denominator> denom_ptr) {
    graph_1_ = graph_1;
    graph_2_ = graph_2;
    ini_ptr_ = std::move(ini_ptr);
    mm_ptr_ = std::move(mm_ptr);
    denom_ptr_ = std::move(denom_ptr);
    n1_ = graph_1_.get_num_vertices();
    n2_ = graph_2_.get_num_vertices();
    l1_ = graph_1_.get_num_labels();
    l2_ = graph_2_.get_num_labels();
    node_label_1_ = graph_1_.get_vertex_type();
    node_label_2_ = graph_2_.get_vertex_type();
    VectorFloatType temp;
    for (uint32_t i = 0; i < n1_; i++) {
        temp = VectorFloatType(n2_, 0.0);
        simmat_.push_back(temp);
    }

    for (uint32_t i = 0; i < l1_; i++) {
        temp = VectorFloatType(l2_, 0.0);
        labelmat_.push_back(temp);
    }
}

//*************************************************************************************
//******************************** Basic operations ***********************************
SimMatrix& FracSimulation::get_curr_simmat() {
    return simmat_;
}

void FracSimulation::simmat_copy_to(SimMatrix& premat) {
    cout << "begin simmat copy" << endl;
    double start, end;
    start = omp_get_wtime();
    premat.assign(simmat_.begin(), simmat_.end());
    end = omp_get_wtime();
    cout << "fracsimulation, simmatrix copy time : " << end - start << " s." << endl;
}

// return value of similarity matrix
float FracSimulation::get_value(uint32_t i, uint32_t j) {
    assert(i < n1_ && j < n2_);
    return simmat_[i][j];
}

// set value for similarity matrix
void FracSimulation::set_value(uint32_t i, uint32_t j, float sim) {
    assert(i < n1_ && j < n2_);
    simmat_[i][j] = sim;
}

// clear content of the similarity matrix and reset all the values to 0.0
void FracSimulation::simmat_clear() {
    VectorFloatType::iterator cell;
    for (uint32_t i = 0; i < n1_; i++) {
        for (cell = simmat_[i].begin(); cell != simmat_[i].end(); cell++) {
            *cell = 0.0;
        }
    }
}

void FracSimulation::print_simmat() {
    for (uint32_t i = 0; i < n1_; i++) {
        for (uint32_t j = 0; j < n2_; j++) {
            cout << setw(7) << setprecision(3) << simmat_[i][j] << " ";
        }
        cout << endl;
    }
}

void FracSimulation::print_labelmat() {
    for (uint32_t i = 0; i < l1_; i++) {
        for (uint32_t j = 0; j < l2_; j++) {
            cout << setw(7) << setprecision(3) << labelmat_[i][j] << " ";
        }
        cout << endl;
    }
}

//*************************************************************************************
//********************************* Initialization ************************************
// initialize the similarity matrix by identity  --- for simrank
void FracSimulation::initialize(uint32_t num_thus) {
    double start, end;
    start = omp_get_wtime();
    ini_ptr_->initialize(labelmat_, graph_1_.get_label_info(), graph_2_.get_label_info(), num_thus);
#pragma omp parallel num_threads(num_thus)
    {
        uint32_t pid = omp_get_thread_num(), np = omp_get_num_threads();
        for (uint32_t i = pid; i < n1_; i += np) {
            for (uint32_t j = 0; j < n2_; j++) {
                simmat_[i][j] = labelmat_[node_label_1_[i]][node_label_2_[j]];
            }
        }
    }
    // cout << "labelmat_: ";
    // print_labelmat();
    cout << "initialize simmat: " << endl;
    print_simmat();
    end = omp_get_wtime();
    cout << "fracsimulation, initialize time : " << end - start << " s." << endl;
}

//*************************************************************************************
//************************************ Iteration **************************************
// calculate the error
bool FracSimulation::diff(SimMatrix& premat, float threshold, string mode) {
    float diff;
    float max_diff = 0.0;
    double total_diff = 0.0;
    double sum = 0;
    for (uint32_t i = 0; i < n1_; i++) {
        for (uint32_t j = 0; j < n2_; j++) {
            diff = abs(simmat_[i][j] - premat[i][j]);
            total_diff += diff;
            if (diff > max_diff) {
                max_diff = diff;
            }
            sum += premat[i][j];
        }
    }
    if (mode == "RELATIVE") {
        cout << "max_diff = " << max_diff << ", total_diff = " << total_diff << ", sum * threshold = " << sum * threshold << endl;
        return (total_diff < (sum * threshold));
    } else if (mode == "MAXDIFF") {
        cout << "max_diff = " << max_diff << ", total_diff = " << total_diff << ", threshold = " << threshold << endl;
        return (max_diff < threshold);
    } else {
        cout << "max_diff = " << max_diff << ", total_diff = " << total_diff << ", threshold = " << threshold << endl;
        return (total_diff < threshold);
    }
}

bool FracSimulation::diff_relative(SimMatrix& premat, float threshold) {
    return diff(premat, threshold, "RELATIVE");
}

bool FracSimulation::diff_totaldiff(SimMatrix& premat, float threshold) {
    return diff(premat, threshold, "TOTALDIFF");
}

bool FracSimulation::diff_maxdiff(SimMatrix& premat, float threshold) {
    return diff(premat, threshold, "MAXDIFF");
}

void FracSimulation::update_simmatrix(SimMatrix& premat, float w_i, float w_o, float w_l, uint32_t num_thus) {
    cout << "update simmatrix" << endl;
    double start, end;
    start = omp_get_wtime();
    vector<uint32_t> in_degree_vec1 = graph_1_.get_in_degree_vec();
    vector<uint32_t> out_degree_vec1 = graph_1_.get_out_degree_vec();
    vector<uint32_t> in_degree_vec2 = graph_2_.get_in_degree_vec();
    vector<uint32_t> out_degree_vec2 = graph_2_.get_out_degree_vec();
#pragma omp parallel num_threads(num_thus)
    {
        uint32_t pid = omp_get_thread_num(), np = omp_get_num_threads();
        for (uint32_t i = pid; i < n1_; i += np) {
            NeighborsMap map_i_out = graph_1_.get_out_neighbors_map(i);
            NeighborsMap map_i_in = graph_1_.get_in_neighbors_map(i);
            for (uint32_t j = 0; j < n2_; j++) {
                float label_sim = 0.0, inneighbor_sim = 0.0, outneighbor_sim = 0.0;
                NeighborsMap map_j_out = graph_2_.get_out_neighbors_map(j);
                NeighborsMap map_j_in = graph_2_.get_in_neighbors_map(j);
                if (w_i > 0 && map_i_in.size() > 0 && map_j_in.size() > 0) {
                    inneighbor_sim = mm_ptr_->match_score(premat, map_i_in, map_j_in) / denom_ptr_->denom_func(in_degree_vec1[i], in_degree_vec2[j]);
                }
                if (map_i_out.size() > 0 && map_j_out.size() > 0) {
                    outneighbor_sim = mm_ptr_->match_score(premat, map_i_out, map_j_out) / denom_ptr_->denom_func(out_degree_vec1[i], out_degree_vec2[j]);
                }
                label_sim = labelmat_[node_label_1_[i]][node_label_2_[j]];

                if (map_i_in.size() == 0 && map_j_in.size() == 0) {
                    inneighbor_sim = 1.0;
                }

                if (map_i_out.size() == 0 && map_j_out.size() == 0) {
                    outneighbor_sim = 1.0;
                }
                simmat_[i][j] = outneighbor_sim * w_o + inneighbor_sim * w_i + label_sim * w_l;
            }
        }
    }
    end = omp_get_wtime();
    cout << "fracsimulation, update simmatrix time : " << end - start << " s." << endl;
}
