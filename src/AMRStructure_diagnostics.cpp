#include "AMRStructure.hpp"
#include <fstream>
#include <iomanip>

MHDDiagnostics AMRStructure::compute_diagnostics() {
    MHDDiagnostics d{};
    d.iter = iter_num;
    d.t    = t;

    const size_t N = weights.size();
    double E_kin = 0.0, E_mag = 0.0, H_C = 0.0;
    double I_j = 0.0, I_w = 0.0;

    #pragma omp parallel for reduction(+:E_kin, E_mag, H_C, I_j, I_w)
    for (size_t i = 0; i < N; ++i) {
        const double wi = weights[i];
        E_kin += 0.5 * wi * (u1s[i]*u1s[i] + u2s[i]*u2s[i]);
        E_mag += 0.5 * wi * (b1s[i]*b1s[i] + b2s[i]*b2s[i]);
        H_C   +=       wi * (u1s[i]*b1s[i] + u2s[i]*b2s[i]);
        I_j   +=       wi * j0s[i];
        I_w   +=       wi * w0s[i];
    }

    d.E_kin = E_kin;
    d.E_mag = E_mag;
    d.E_tot = E_kin + E_mag;
    d.H_C   = H_C;
    d.I_j   = I_j;
    d.I_w   = I_w;
    return d;
}

int AMRStructure::write_diagnostics(const MHDDiagnostics& d) {
    static bool header_written = false;
    const std::string path = sim_dir + "simulation_output/diagnostics.csv";

    std::ofstream f;
    if (!header_written) {
        f.open(path, std::ios::out | std::ios::trunc);
        f << "iter,t,E_kin,E_mag,E_tot,H_C,I_j,I_w\n";
        header_written = true;
    } else {
        f.open(path, std::ios::out | std::ios::app);
    }

    if (!f.is_open()) {
        std::cerr << "[diag] failed to open " << path << std::endl;
        return 1;
    }

    f << std::setprecision(16);
    f << d.iter << "," << d.t << ","
      << d.E_kin << "," << d.E_mag << "," << d.E_tot << ","
      << d.H_C << "," << d.I_j << "," << d.I_w << "\n";
    return 0;
}