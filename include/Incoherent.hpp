#pragma once


#include <cuba.h>
#include <tuple>
#include "IntegrationRoutines.hpp"


namespace Incoherent {
    double A_integrand_function_factor(double Q);
    double A_integrand_function(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    int integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);
    double dsigmadt(double Q, double Delta);
    double dsigmadt(CubaConfig* cuba_config, IntegrationConfig* integration_config);

    int integrand_cubature(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigmadt_cubature(double Q, double Delta);
    double dsigmadt_cubature(CubatureConfig* cuba_config, IntegrationConfig* integration_config);
}


namespace Incoherent { namespace GeometryAverage
{
    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigmadt_single_event (double Q, double Delta, const HotspotNucleus& nucleus);
    double dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config);
} }