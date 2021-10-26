#include "Denominator.h"

float QuerySide::denom_func(uint32_t i, uint32_t j) {
    return float(i);
}

float Add::denom_func(uint32_t i, uint32_t j) {
    return float(i + j);
}

float Root::denom_func(uint32_t i, uint32_t j) {
    return float(sqrt(i) * sqrt(j));
}

float Max::denom_func(uint32_t i, uint32_t j) {
    return float(i > j ? i : j);
}

float HalfAdd::denom_func(uint32_t i, uint32_t j) {
    return float((i + j) / 2.0);
}
