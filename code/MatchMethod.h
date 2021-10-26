#ifndef _MATCHMETHOD_H
#define _MATCHMETHOD_H

#include <assert.h>

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Graph.h"

using namespace std;

typedef vector<vector<float>> SimMatrix;
typedef vector<uint32_t> Neighbors;
typedef unordered_map<uint32_t, Neighbors> NeighborsMap;
typedef vector<NeighborsMap> AdjMap;
typedef vector<unordered_map<uint32_t, float>> SimMap;

struct MatchItem {
    uint32_t x;
    uint32_t y;
    float weight;
};

template <typename T>
struct WeightLessThan {  // new for SimMatrix.cpp
    bool operator()(const T& a, const T& b) {
        return (a.weight < b.weight);
    }
};

class MatchMethod {
  protected:
    string label_constrainted_;

  public:
    MatchMethod();
    MatchMethod(string flag);

    float match_score(SimMatrix& basemat, NeighborsMap& map_i, NeighborsMap& map_j);
    virtual float match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j) = 0;
};

class InjectMatch : public MatchMethod {
  public:
    InjectMatch() : MatchMethod(){};
    InjectMatch(string flag) : MatchMethod(flag){};
    float match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j);
};

class BiMatch : public MatchMethod {
  public:
    BiMatch() : MatchMethod(){};
    BiMatch(string flag) : MatchMethod(flag){};
    float match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j);
};

class BijectMatch : public MatchMethod {
  public:
    BijectMatch() : MatchMethod(){};
    BijectMatch(string flag) : MatchMethod(flag){};
    float match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j);
};
#endif
