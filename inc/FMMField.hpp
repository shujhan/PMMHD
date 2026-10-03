#ifndef FMM_FIELD_HPP
#define FMM_FIELD_HPP

#include <cstddef>
#include "FieldStructure.hpp"

// BarytreeK (Kokkos) field solver behind the Field interface.
// This header has no Kokkos types, so the OpenACC/nvc++ translation units can
// include it; the implementation (fmm/FMMField.cpp) is built as a separate
// Kokkos library.
//
// Same call convention as U_Treecode: the ny entries of x_vals/y_vals/q_ws are
// the sources, the first nx of them are the targets.
// Supported modes: periodic_xy (free-space Biot-Savart, u1/u2) and
// periodic_xy_potentials (free-space log potential, e1s only).

void fmm_initialize(int& argc, char** argv);   // wraps Kokkos::initialize
void fmm_finalize();                           // wraps Kokkos::finalize

class U_FMM : public Field {
    public:
        U_FMM();
        U_FMM(double epsilon, double mac, int degree, int max_source);
        void operator() (double* e1s, double* e2s, double* x_vals, int nx,
                        double* y_vals, double* q_ws, int ny);
        void print_field_obj();
        void set_mode(KernelMode m) override;
        ~U_FMM();

        // allocate and free U_FMM inside the FMM library, so `new U_FMM` in nvc++
        // code (-gpu=managed) and the delete from the nvcc-compiled destructor use
        // the same allocator
        static void* operator new(size_t size);
        static void operator delete(void* p);

    private:
        double epsilon;
        double mac;       // MAC theta, (r_t + r_s) / R < mac
        int degree;       // Chebyshev interpolation degree (<= 10)
        int max_source;   // leaf size (cluster_size)
        KernelMode mode;
};

#endif /* FMM_FIELD_HPP */