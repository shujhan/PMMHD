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

ElectricField::~ElectricField() = default;


E_MQ_DirectSum::E_MQ_DirectSum() {}
E_MQ_DirectSum::E_MQ_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
E_MQ_DirectSum::~E_MQ_DirectSum() = default;

void E_MQ_DirectSum::operator() (double* e1s, double* e2s, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
    // double epsLsq = epsilon * epsilon;
    for (int i = 0; i < nt; i++) {
        // e1s[i] = 0.1 * q_ws[i];
        // e2s[i] = 0.2 * q_ws[i];
        e1s[i] = 0.1 * q_ws[i];
        e2s[i] = 0.2 * q_ws[i];
    }


}


void E_MQ_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
}
