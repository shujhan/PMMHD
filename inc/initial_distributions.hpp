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

class w0_alfven : public distribution {
    public:
        w0_alfven(double kx_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
};


class w0_polarized_alfven : public distribution {
    public:
        w0_polarized_alfven(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double k_norm;
    double amp;
};

class w0_orszag_tang : public distribution {
    public:
        w0_orszag_tang(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};



// Ideal/inviscid Orszag-Tang (Kraus & Maj 2018):
//   phi = 2 sin(ky*y) - 2 cos(kx*x),  u0 = curl(phi)
//   => omega = laplacian(phi) = 2 ky^2 sin(ky*y) - 2 kx^2 cos(kx*x)
class w0_kraus_maj : public distribution {
    public:
        w0_kraus_maj(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};



// Coalescence seed (Bhattacharjee et al. 2009 setup):
//   phi = amp * (cos(kx*x) - cos(ky*y)),  u0 = curl(phi)
//   => omega = laplacian(phi) = amp * (ky^2 cos(ky*y) - kx^2 cos(kx*x))
// Small amp (~1e-3) breaks the symmetry of the four-flux-tube equilibrium
// and pushes same-sign current islands toward each other along the diagonal.
class w0_coalescence : public distribution {
    public:
        w0_coalescence(double kx_w, double ky_w, double amp_w);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};





class j0_uniform : public distribution {
    public:
        j0_uniform();

        double operator() (double x, double y);
        void print();
};


class j0_current_sheet: public distribution {
    public:
        j0_current_sheet(double kx_j, double amp_j, double thickness);

        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
    double a;
};


class j0_alfven: public distribution {
    public:
        j0_alfven(double kx_j, double amp_j);

        double operator() (double x, double y);
        void print();
    double kx;
    double amp;
};


class j0_polarized_alfven : public distribution {
    public:
        j0_polarized_alfven(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double k_norm;
    double amp;
};


class j0_orszag_tang : public distribution {
    public:
        j0_orszag_tang(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};



// Ideal/inviscid Orszag-Tang (Kraus & Maj 2018):
//   psi = cos(2 ky*y) - 2 cos(kx*x)  (out-of-plane vector potential A_z)
//   j = -laplacian(psi) = 4 ky^2 cos(2 ky*y) - 2 kx^2 cos(kx*x)
class j0_kraus_maj : public distribution {
    public:
        j0_kraus_maj(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};



// Coalescence / four-flux-tube equilibrium (Bhattacharjee et al. 2009):
//   psi = A0 sin(kx*x) sin(ky*y)
//   => j = -laplacian(psi) = A0 (kx^2 + ky^2) sin(kx*x) sin(ky*y)
// Pass amp = A0 (kx^2 + ky^2); for kx = ky = 1 and A0 = 1, amp = 2.
class j0_coalescence : public distribution {
    public:
        j0_coalescence(double kx_j, double ky_j, double amp_j);
        double operator() (double x, double y);
        void print();
    double kx;
    double ky;
    double amp;
};

#endif