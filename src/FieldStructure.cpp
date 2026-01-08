#include "FieldStructure.hpp"
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cfloat> // dbl_min
#include <cstddef>
using namespace std;
#if OPENACC_ENABLED
#include <accelmath.h>
#endif

Field::~Field() = default;


U_DirectSum::U_DirectSum() {}
U_DirectSum::U_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
U_DirectSum::~U_DirectSum() = default;

void U_DirectSum::operator() (double* u1s, double* u2s, double* x_vals, int nx, 
                        double* y_vals, double* q_ws, int ny)
{    
    const double pi = std::atan(1.0) * 4.0;
    for (int i = 0; i < nx; i++) {
        for(int k = 0; k < ny; k++) {
            double denom = cosh(2* pi / L * (y_vals[i] - y_vals[k])) - cos(2* pi / L * (x_vals[i] - x_vals[k])) + epsilon * epsilon;
            u1s[i] -= 0.5/L * sinh(2 * pi / L * (y_vals[i] - y_vals[k])) / denom * q_ws[k];
            u2s[i] += 0.5/L * sin(2 * pi / L * (x_vals[i] - x_vals[k])) / denom * q_ws[k];
        }
    }
}


void U_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
}
