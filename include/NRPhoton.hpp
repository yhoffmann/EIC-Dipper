#pragma once

#include "IntegrationRoutines.hpp"


namespace NRPhoton
{
    double epsilon(double Q);

    double get_wave_function_factor_T(double m_Q, double e_Q);
    void set_wave_function_factor_T(double m_Q, double e_Q);
    inline double wave_function_factor_T = get_wave_function_factor_T(m_Q, e_Q);

    double wave_function_factor_L(double Q);

    double wave_function_factor(double Q);

    double wave_function(double r1, double r2, double Q);
}