#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include "../include/IntegrationRoutines.hpp"

namespace NRPhoton
{
    double epsilon (double Q, double z)
    {
        return std::sqrt( sqr(Q)*z*(1.0-z) + sqr(m_Q_c) );
    }

    double wave_function (double r1, double r2, double Q, double z, TransverseOrLongitudinal transverse_or_longitudinal)
    {
        if (sqr(r1)+sqr(r2) < 2.0e-40) { r1 = 1.0e-20; r2 = 1.0e-20;}
        
        if (transverse_or_longitudinal==T)
        {
            return -A_Q * sqrt_2m_c_Nc * e * e_Q * gsl_sf_bessel_K0( epsilon(Q,z) * std::sqrt( sqr(r1)+sqr(r2) ) );
        }
        else
        {
            return -2.0 / m_Q_c * Q * z * (1.0-z) * A_Q * sqrt_2m_c_Nc * e * e_Q * gsl_sf_bessel_K0( epsilon(Q,z) * std::sqrt( sqr(r1)+sqr(r2) ) );
        }
    }
}