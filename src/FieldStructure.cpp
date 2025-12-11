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

void E_MQ_DirectSum::operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
    double epsLsq = epsilon * epsilon;

}


void E_MQ_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
}
