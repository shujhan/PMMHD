#include "initial_distributions.hpp"

// ---- w0_uniform ----
w0_uniform::w0_uniform() {}

double w0_uniform::operator()(double x, double y) {
    return 1.0;
}

void w0_uniform::print() {
    std::cout << "w0_uniform distribution: 1s" << std::endl;
}


// ---- w0_current_sheet ----
w0_current_sheet::w0_current_sheet(double kx, double amp, double thickness):
    kx(kx), amp(amp), a(thickness){}

double w0_current_sheet::operator()(double x, double y) {
    // return amp * cos(kx * x);
    double xi = y / a;
    double T  = std::tanh(xi);
    double s  = 1.0 / (std::cosh(xi) * std::cosh(xi));
    double E  = std::exp(-xi * xi);
    double F   = T * E;
    double Fpp = E * (-2.0*s*T - 4.0*xi*s + (4.0*xi*xi - 2.0)*T);
    return amp * std::sin(kx * x) * (Fpp / (a * a * kx) - kx * F);
}

void w0_current_sheet::print() {
    std::cout << "w0_current_sheet distribution: amp * cos(kx * x) " << std::endl;
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


// ---- w0_polarized_alfven----
w0_polarized_alfven::w0_polarized_alfven(double kx_w, double ky_w, double amp_w):
    kx(kx_w), ky(ky_w), amp(amp_w) {
        k_norm = sqrt(kx*kx + ky*ky);
    }

double w0_polarized_alfven::operator()(double x, double y) {
    return amp * k_norm * cos(kx * x + ky* y);
}

void w0_polarized_alfven::print() {
    std::cout << "w0_polarized_alfven: amp * k_norm * cos(kx * x + ky * y)" << std::endl;
}



// ---- w0_orszag_tang ----
w0_orszag_tang::w0_orszag_tang(double kx_w, double ky_w, double amp_w):
    kx(kx_w), ky(ky_w), amp(amp_w) {}

double w0_orszag_tang::operator()(double x, double y) {
    return amp * (cos(kx * x) + cos(ky * y));
}

void w0_orszag_tang::print() {
    std::cout << "w0_orszag_tang: amp * (cos(kx * x) + cos(ky * y))" << std::endl;
}



// // ---- w0_kraus_maj ----
// // omega = laplacian(phi), phi = 2 sin(ky*y) - 2 cos(kx*x)
// w0_kraus_maj::w0_kraus_maj(double kx_w, double ky_w, double amp_w):
//     kx(kx_w), ky(ky_w), amp(amp_w) {}

// double w0_kraus_maj::operator()(double x, double y) {
//     return amp * (sin(ky * y) - cos(kx * x));
// }

// void w0_kraus_maj::print() {
//     std::cout << "w0_kraus_maj: amp * (-2 ky^2 sin(ky*y) + 2 kx^2 cos(kx*x))" << std::endl;
// }

// ---- w0_kraus_maj ----
// omega = laplacian(phi), phi = 2 sin(ky*y) - 2 cos(kx*x)
w0_kraus_maj::w0_kraus_maj(double kx_w, double ky_w, double amp_w):
    kx(kx_w), ky(ky_w), amp(amp_w) {}

double w0_kraus_maj::operator()(double x, double y) {
    return amp * (-2.0 * ky * ky * sin(ky * y) + 2.0 * kx * kx * cos(kx * x));
}

void w0_kraus_maj::print() {
    std::cout << "w0_kraus_maj: amp * (-2 ky^2 sin(ky*y) + 2 kx^2 cos(kx*x))" << std::endl;
}



// ---- w0_coalescence ----
// omega = laplacian(phi), phi = amp * (cos(kx*x) - cos(ky*y))
w0_coalescence::w0_coalescence(double kx_w, double ky_w, double amp_w):
    kx(kx_w), ky(ky_w), amp(amp_w) {}

double w0_coalescence::operator()(double x, double y) {
    return amp * (ky * ky * cos(ky * y) - kx * kx * cos(kx * x));
}

void w0_coalescence::print() {
    std::cout << "w0_coalescence: amp * (ky^2 cos(ky*y) - kx^2 cos(kx*x))" << std::endl;
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
j0_current_sheet::j0_current_sheet(double kx_j, double amp_j, double thickness):
    kx(kx_j), amp(amp_j) , a(thickness){}

double j0_current_sheet::operator()(double x, double y) {
    return 1.0 / a / (std::cosh(y/a) * std::cosh(y/a)) * ( 1 + amp * cos(kx * x));
}


void j0_current_sheet::print() {
    std::cout << "j0_current_sheet distribution: 1.0 / a / (cosh(y/a))^2 * ( 1 + amp * cos(kx * x))" << std::endl;
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


// ---- j0_polarized_alfven----
j0_polarized_alfven::j0_polarized_alfven(double kx_j, double ky_j, double amp_j):
    kx(kx_j), ky(ky_j), amp(amp_j) {
        k_norm = sqrt(kx*kx + ky*ky);
    }

double j0_polarized_alfven::operator()(double x, double y) {
    return amp * k_norm * cos(kx* x + ky* y);
}

void j0_polarized_alfven::print() {
    std::cout << "j0_polarized_alfven: amp * k_norm * cos(kx * x + ky * y)" << std::endl;
}


// ---- j0_orszag_tang ----
j0_orszag_tang::j0_orszag_tang(double kx_j, double ky_j, double amp_j):
    kx(kx_j), ky(ky_j), amp(amp_j) {}

double j0_orszag_tang::operator()(double x, double y) {
    return amp * (cos(kx * x) + 2 * cos(ky * y));
}

void j0_orszag_tang::print() {
    std::cout << "j0_orszag_tang: amp * (cos(kx * x) + cos(ky * y))" << std::endl;
}


// // ---- j0_kraus_maj ----
// // j = -laplacian(psi), psi = cos(2 ky*y) - 2 cos(kx*x)
// j0_kraus_maj::j0_kraus_maj(double kx_j, double ky_j, double amp_j):
//     kx(kx_j), ky(ky_j), amp(amp_j) {}

// double j0_kraus_maj::operator()(double x, double y) {
//     return amp * (2 * cos(ky * y) - cos(kx * x));
// }

// void j0_kraus_maj::print() {
//     std::cout << "j0_kraus_maj: amp * (4 ky^2 cos(2 ky*y) - 2 kx^2 cos(kx*x))" << std::endl;
// }

// ---- j0_kraus_maj ----
// j = -laplacian(psi), psi = cos(2 ky*y) - 2 cos(kx*x)
j0_kraus_maj::j0_kraus_maj(double kx_j, double ky_j, double amp_j):
    kx(kx_j), ky(ky_j), amp(amp_j) {}

double j0_kraus_maj::operator()(double x, double y) {
    return amp * (4.0 * ky * ky * cos(2.0 * ky * y) - 2.0 * kx * kx * cos(kx * x));
}

void j0_kraus_maj::print() {
    std::cout << "j0_kraus_maj: amp * (4 ky^2 cos(2 ky*y) - 2 kx^2 cos(kx*x))" << std::endl;
}



// ---- j0_coalescence ----
// j = -laplacian(psi), psi = A0 sin(kx*x) sin(ky*y), amp = A0 (kx^2 + ky^2)
j0_coalescence::j0_coalescence(double kx_j, double ky_j, double amp_j):
    kx(kx_j), ky(ky_j), amp(amp_j) {}

double j0_coalescence::operator()(double x, double y) {
    return amp * sin(kx * x) * sin(ky * y);
}

void j0_coalescence::print() {
    std::cout << "j0_coalescence: amp * sin(kx*x) * sin(ky*y)" << std::endl;
}