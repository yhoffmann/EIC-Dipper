#pragma once


namespace NRPhoton {
    extern double ONE_DIV_EPSILON_MIN;

    double epsilon (double Q, double z);

    double wave_function (double r1, double r2, double Q, double z, TransverseOrLongitudinal transverse_or_longitudinal);
}