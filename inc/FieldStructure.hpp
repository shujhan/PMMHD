#ifndef FIELD_STRUCTURE_HPP
#define FIELD_STRUCTURE_HPP
#include <math.h>

using namespace std;

class ElectricField {
    public: 
        virtual void operator()     (double* es, double* targets, int nt, 
                                    double* sources, double* q_ws, int ns) = 0;
        virtual void print_field_obj() = 0;
        virtual ~ElectricField();
};

class E_MQ_DirectSum : public ElectricField {
    public:
        E_MQ_DirectSum();
        E_MQ_DirectSum(double L, double epsilon);
        double epsilon;
        double L;
        void operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns);
        void print_field_obj();
        ~E_MQ_DirectSum();
};


#endif /* FIELD_STRUCTURE_HPP */