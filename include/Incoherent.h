#pragma once


#include <cuba.h>
#include <vector>
#include "IntegrationRoutines.h"


namespace Incoherent {
    double A_integrand_function_simple (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal);

    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal);

    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);

    std::vector<double> dsigma_dt (CubaConfig* cuba_config, IntegrationConfig* integration_config);
}