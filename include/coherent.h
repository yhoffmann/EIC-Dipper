#pragma once


#include <cuba.h>
#include <vector>
#include "integration_routines.h"
#include "models_and_functions.h"


namespace Coherent {
    double dsigma_d2b (double b1, double b2, double r1, double r2);

    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double z, double Delta, bool t_not_l);

    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);

    std::vector<double> calculate_dsigma_dt (CubaConfig,IntegrandParams);
}