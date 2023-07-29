#pragma once


#include <gsl/gsl_sf.h>


typedef long unsigned int luint;

inline double get_b_range_factor() {
    return 15.0;
}

inline double get_r_range_factor() {
    return 20.0;
}


inline double sqr (double x) {
    return x*x;
}


double bessel_K_safe (int n, double x);