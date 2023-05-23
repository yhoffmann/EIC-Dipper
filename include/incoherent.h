#pragma once


#include <cuba.h>
#include <vector>
#include "integration_routines.h"


namespace Incoherent {
    double DD (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

    double dsigma_d2b_sqr (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

    double A_integrand_function_simple (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, bool t_not_l);

    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, bool t_not_l);

    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);

    std::vector<double> calculate_dsigma_dt (CubaConfig c_config, IntegrandParams i_params);
}