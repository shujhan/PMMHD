#ifndef INITIAL_DISTRIBUTIONS_HPP
#define INITIAL_DISTRIBUTIONS_HPP

#include <vector> 
#include <iostream> 
#include <cmath>

class distribution {
    public:
        virtual double operator() (double x, double y)=0;
        virtual void print()=0;
};

class w0_uniform : public distribution {
    public:
        w0_uniform();
        double operator() (double x, double y);
        void print();
};

class w0_zero : public distribution {
    public:
        w0_zero();
        double operator() (double x, double y);
        void print();
};

class j0_uniform : public distribution {
    public:
        j0_uniform();

        double operator() (double x, double y);
        void print();
};

class j0_current_sheet: public distribution {
    public:
        j0_current_sheet();

        double operator() (double x, double y);
        void print();
};


#endif
