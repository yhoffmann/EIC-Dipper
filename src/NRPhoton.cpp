#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "../include/NRPhoton.hpp"
#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include "../include/IntegrationRoutines.hpp"

namespace NRPhoton
{
    double epsilon (double Q)
    {
        return std::sqrt( sqr(Q)*0.25 + sqr(m_Q) );
    }


    double get_wave_function_factor_T (double m_Q, double e_Q)
    {
        return -A_Q * std::sqrt(2.0*m_Q*Nc) * e * e_Q;
    }

    void set_wave_function_factor_T (double m_Q, double e_Q)
    {
        wave_function_factor_T = get_wave_function_factor_T(m_Q, e_Q);
    }


    double wave_function_factor_L (double Q)
    {
        return 2.0 * wave_function_factor_T * Q / m_Q;
    }


    double wave_function_factor (double Q)
    {
        return std::sqrt( sqr(wave_function_factor_T) + sqr( wave_function_factor_L(Q) ) );
    }


    double wave_function (double r1, double r2, double Q)
    {
        double rsqr = sqr(r1) + sqr(r2);

        if (rsqr == 0.0)
            rsqr = 2.0e-15;

        return gsl_sf_bessel_K0( epsilon(Q) * std::sqrt(rsqr) );
    }
}