#include "initial_distributions.hpp"

// ---- w0_uniform ----
w0_uniform::w0_uniform() {}

double w0_uniform::operator()(double x, double y) {
    return 1.0;
}

void w0_uniform::print() {
    std::cout << "w0_uniform distribution: 1s" << std::endl;
}


// ---- w0_zero ----
w0_zero::w0_zero() {}

double w0_zero::operator()(double x, double y) {
    return 0.0;
}

void w0_zero::print() {
    std::cout << "w0_zerodistribution: 0s" << std::endl;
}




// ---- j0_uniform ----
j0_uniform::j0_uniform() {}

double j0_uniform::operator()(double x, double y) {
    return 1.0;
}

void j0_uniform::print() {
    std::cout << "j0_uniform distribution: 1s" << std::endl;
}

// ---- current sheet ----
j0_current_sheet::j0_current_sheet() {}

double j0_current_sheet::operator()(double x, double y) {
    return 1.0 / (std::cosh(y) * std::cosh(y));
}


void j0_current_sheet::print() {
    std::cout << "j0_current_sheet distribution: 1.0 / (cosh(y))^2" << std::endl;
}

