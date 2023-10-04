#pragma once


#include <cuba.h>
#include <tuple>
#include "IntegrationRoutines.hpp"


namespace Incoherent {
    double A_integrand_function_simple(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    double A_integrand_function(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    int integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);
    double dsigma_dt(double Q, double Delta);
    double dsigma_dt(CubaConfig* cuba_config, IntegrationConfig* integration_config);

    int integrand_cubature(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigma_dt_cubature(double Q, double Delta);
    double dsigma_dt_cubature(CubatureConfig* cuba_config, IntegrationConfig* integration_config);
}


namespace Incoherent { namespace DiluteApprox
{
    double A_integrand_function(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    int integrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
    double dsigma_dt(double Q, double Delta);
    double dsigma_dt(CubaConfig* c_config, IntegrationConfig* integration_config);
} }