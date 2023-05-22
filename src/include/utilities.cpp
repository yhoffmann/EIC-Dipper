#pragma once

#include"constants.cpp"

typedef long unsigned int luint;

inline double get_b_range_factor() {
    return 15.0;
}

inline double get_r_range_factor() {
    return 20.0;
}


struct CubaConfig {
    char integrator = 'C';
    int num_of_dims;
    double epsrel = 1.0e-6;
    double epsabs = 1.0e-8;
    luint mineval = 2e5;
    luint maxeval = 1e6;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;
    double bessel_tolerance = 1e-5;
    int n_new = 1000;
    int n_min = 2;
};


struct IntegrandParams {
    double (*func)(TEN_D);

    double Q = 0.3;
    const double z = 0.5;

    double min = -999;
    double max = std::sqrt(get_b_range_factor()*2.0*BG);
};