#ifndef _DENOMINATOR_H
#define _DENOMINATOR_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Graph.h"
#include "MatchMethod.h"

using namespace std;

typedef vector<vector<float>> SimMatrix;
typedef vector<float> VectorFloatType;

class Denominator {
  public:
    // virtual ~Denominator() = 0;
    virtual float denom_func(uint32_t, uint32_t) = 0;
};

class QuerySide : public Denominator {
  public:
    float denom_func(uint32_t, uint32_t);
};

class Add : public Denominator {
  public:
    float denom_func(uint32_t, uint32_t);
};

class Root : public Denominator {
  public:
    float denom_func(uint32_t, uint32_t);
};

class Max : public Denominator {
  public:
    float denom_func(uint32_t, uint32_t);
};

class HalfAdd : public Denominator {
  public:
    float denom_func(uint32_t, uint32_t);
};
#endif
