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
        return std::sqrt( sqr(Q)*0.25 + sqr(m_Q_c) );
    }


    double wave_function_factor_T = -A_Q * sqrt_2m_c_Nc * e * e_Q;


    double wave_function_factor_L (double Q)
    {
        return wave_function_factor_T * Q / (2.0*m_Q_c);
    }


    double wave_function_factor (double Q)
    {
        return std::sqrt(sqr(wave_function_factor_T) + sqr(wave_function_factor_L(Q)));
    }


    double wave_function (double r1, double r2, double Q)
    {
        if (sqr(r1)+sqr(r2) < 2.0e-40)
        {
            r1 = 1.0e-20;
            r2 = 1.0e-20;
        }
        return gsl_sf_bessel_K0( epsilon(Q) * std::sqrt( sqr(r1)+sqr(r2) ) );
    }
}