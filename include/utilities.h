#pragma once


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


inline double x (double b, double r) {
    return b+r/2.0;
}

inline double y (double b, double r) {
    return b-r/2.0;
}