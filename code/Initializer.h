#ifndef _INITIALIZER_H
#define _INITIALIZER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

class Initializer {
  public:
    // virtual void initialize(vector<vector<float>>&, vector<string>&, vector<string>&) = 0;
    virtual float cal_string_sim(const string&, const string&) = 0;
    void initialize(vector<vector<float>>&, const vector<string>&, const vector<string>&, uint32_t);
    int my_min(int, int, int);
};

class IdenticalInitializer : public Initializer {
  public:
    float cal_string_sim(const string&, const string&);
};

class JaccardInitializer : public Initializer {
  public:
    float cal_string_sim(const string&, const string&);
};

class HammingInitializer : public Initializer {
  public:
    float cal_string_sim(const string&, const string&);
};

class EditDistInitializer : public Initializer {
  public:
    float cal_string_sim(const string&, const string&);
};

class JarowinklerInitializer : public Initializer {
  public:
    float cal_string_sim(const string&, const string&);
};

#endif
