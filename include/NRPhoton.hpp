#pragma once


#include "IntegrationRoutines.hpp"


namespace NRPhoton
{
    double epsilon(double Q);

    extern double wave_function_factor_T;

    double wave_function_factor_L(double Q);

    double wave_function_factor(double Q);

    double wave_function(double r1, double r2, double Q);
}