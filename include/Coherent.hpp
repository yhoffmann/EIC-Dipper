#pragma once


#include <cuba.h>
#include <vector>
#include "IntegrationRoutines.h"


namespace Coherent {
    double dsigma_d2b (double x1, double x2, double y1, double y2);

    double A_integrand_function (double x1, double x2, double y1, double y2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal);

    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);

    std::vector<double> dsigma_dt (CubaConfig* cuba_config, IntegrationConfig* integration_config);
}