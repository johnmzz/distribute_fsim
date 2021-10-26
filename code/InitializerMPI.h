#ifndef _INITIALIZERMPI_H
#define _INITIALIZERMPI_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

typedef vector<vector<float>> SimMatrix;
typedef vector<float> VectorFloatType;

class InitializerMPI {
  private:
    string method_;
    vector<string> label_info_1_;  // eg. [label_0, label_1, label_2, ...]
    vector<string> label_info_2_;  // ^

  public:
    SimMatrix label_sim_matrix_;  // label_size * label_size
    InitializerMPI(string method, const vector<string>& label_info_1, const vector<string>& label_info_2);
    int my_min(int, int, int);
    float cal_string_sim(const string& s1, const string& s2);
    float get_sim_value(uint32_t i, uint32_t j);
    void print_label_info();
    void print_label_sim();
};
#endif
