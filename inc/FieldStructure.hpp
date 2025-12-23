#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP
#include <math.h>

using namespace std;

class Field {
    public: 
        virtual void operator()     (double* e1s, double* e2s, double* x_vals, int nx, 
                                    double* y_vals, double* q_ws, int ny) = 0;
        virtual void print_field_obj() = 0;
        virtual ~Field();
};

class U_DirectSum : public Field {
    public:
        U_DirectSum();
        U_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double* e1s, double* e2s, double* x_vals, int nx, 
                        double* y_vals, double* q_ws, int ny);
        void print_field_obj();
        ~U_DirectSum();
};


#endif /* FIELD_STRUCTURE_HPP */
