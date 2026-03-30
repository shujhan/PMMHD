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
    std::cout << "w0_zero distribution: 0s" << std::endl;
}

// ---- w0_alfven----
w0_alfven::w0_alfven(double kx_w, double amp_w):
    kx(kx_w), amp(amp_w) {}

double w0_alfven::operator()(double x, double y) {
    return -1.0 * amp * kx * sin(kx * x);
}

void w0_alfven::print() {
    std::cout << "w0_alfven distribution: -1.0 * amp * kx * sin(kx * x)" << std::endl;
}






// ---- j0_uniform ----
j0_uniform::j0_uniform() {}

double j0_uniform::operator()(double x, double y) {
    return 1.0;
}

void j0_uniform::print() {
    std::cout << "j0_uniform distribution: 1s" << std::endl;
}

// ---- j0_uniform ----
j0_uniform_y::j0_uniform_y() {}

double j0_uniform_y::operator()(double x, double y) {
    return y;
}

void j0_uniform_y::print() {
    std::cout << "j0_uniform_y distribution: y" << std::endl;
}

// ---- current sheet ----
j0_current_sheet::j0_current_sheet(double kx_j, double amp_j):
    kx(kx_j), amp(amp_j) {}

double j0_current_sheet::operator()(double x, double y) {
    return 1.0 / (std::cosh(y) * std::cosh(y)) * ( 1 + amp * cos(kx * x));
}


void j0_current_sheet::print() {
    std::cout << "j0_current_sheet distribution: 1.0 / (cosh(y))^2 * ( 1 + amp * cos(kx * x))" << std::endl;
}


// ---- j0_alfven----
j0_alfven::j0_alfven(double kx_j, double amp_j):
    kx(kx_j), amp(amp_j) {}

double j0_alfven::operator()(double x, double y) {
    return -1.0 * amp * kx * sin(kx * x);
}

void j0_alfven::print() {
    std::cout << "j0_alfven distribution: -1.0 * amp * kx * sin(kx * x)" << std::endl;
}
