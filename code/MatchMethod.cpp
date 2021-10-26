#include "MatchMethod.h"

//*************************************************************************************
//*********************************** Base Class **************************************
MatchMethod::MatchMethod() {
    label_constrainted_ = "N";
}

MatchMethod::MatchMethod(string flag) {
    label_constrainted_ = flag;
}

float MatchMethod::match_score(SimMatrix& basemat, NeighborsMap& map_i, NeighborsMap& map_j) {
    vector<uint32_t> sub_i, sub_j, neighbors_i, neighbors_j;
    NeighborsMap::iterator iter, jter;
    if (map_i.size() == 0 && map_j.size() == 0) {
        return 1.0;
    } else if (map_i.size() == 0 || map_j.size() == 0) {
        return 0.0;
    }
    float match_score = 0.0;
    if (label_constrainted_ == "Y") {
        for (iter = map_i.begin(); iter != map_i.end(); iter++) {
            sub_i = iter->second;
            if (map_j.find(iter->first) != map_j.end()) {
                sub_j = map_j[iter->first];
                match_score += match_score_without_label_constraint(basemat, sub_i, sub_j);
            }
        }
    } else {
        for (iter = map_i.begin(); iter != map_i.end(); iter++) {
            neighbors_i.insert(neighbors_i.end(), iter->second.begin(), iter->second.end());
        }

        for (jter = map_j.begin(); jter != map_j.end(); jter++) {
            neighbors_j.insert(neighbors_j.end(), jter->second.begin(), jter->second.end());
        }
        match_score = match_score_without_label_constraint(basemat, neighbors_i, neighbors_j);
    }
    return match_score;
}

//**************************************************************************************
//********************************* Injective Match ************************************
float InjectMatch::match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j) {
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
            weight = basemat[*iItr][*jItr];
            if (weight > max_score) {
                max_score = weight;
            }
        }
        numerator += max_score;
    }
    return numerator;
}

//**************************************************************************************
//************************************ Bi-Match ****************************************
float BiMatch::match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j) {
    float weight;
    uint32_t size_i = N_i.size();
    uint32_t size_j = N_j.size();
    if (size_i == 0 && size_j == 0) {
        return 1.0;
    } else if (size_i == 0 || size_j == 0) {
        return 0.0;
    }
    vector<uint32_t>::iterator iItr, jItr;
    MatchItem matchitem;
    float numerator = 0.0;
    for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
        float max_score = 0.0;
        for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
            weight = basemat[*iItr][*jItr];
            if (weight > max_score) {
                max_score = weight;
            }
        }
        numerator += max_score;
    }

    for (jItr = N_j.begin(); jItr != N_j.end(); jItr++) {
        float max_score = 0.0;
        for (iItr = N_i.begin(); iItr != N_i.end(); iItr++) {
            weight = basemat[*iItr][*jItr];
            if (weight > max_score) {
                max_score = weight;
            }
        }
        numerator += max_score;
    }
    return numerator;
}

//**************************************************************************************
//********************************** Biject Match **************************************
float BijectMatch::match_score_without_label_constraint(SimMatrix& basemat, Neighbors& N_i, Neighbors& N_j) {
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
            weight = basemat[*iItr][*jItr];
            if (weight > 0.0) {
                matchitem.x = x;
                matchitem.y = y;
                matchitem.weight = weight;
                matchQueue.push(matchitem);
                // cout << "( " << *iItr << ", " << *jItr << "), (" << x << ", " << y << "), " << weight << endl;
            }
            y += 1;
        }
        x += 1;
    }
    // Find the top weighted match-pairs
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
            // cout << x_pos << ", " << y_pos << matchitem.weight << " is selected "<< endl;
        }
        matchQueue.pop();
    }
    return numerator;
}