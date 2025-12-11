#include "initial_distributions.hpp"

// ---- w0_uniform ----
w0_uniform::w0_uniform() {}

double w0_uniform::operator()(double x, double y) {
    return 1.0;
}

void w0_uniform::print() {
    std::cout << "w0_uniform distribution: constant 1" << std::endl;
}


// ---- j0_uniform ----
j0_uniform::j0_uniform() {}

double j0_uniform::operator()(double x, double y) {
    return 1.0;
}

void j0_uniform::print() {
    std::cout << "j0_uniform distribution: constant 1" << std::endl;
}