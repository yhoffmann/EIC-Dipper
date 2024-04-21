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


namespace Incoherent { namespace Demirci
{
    double one_connected_factor (double Q);
    double two_connected_factor (double Q);

    double K_integrand_function (double k, double phik, double kb, double phikb, double A, double Delta, double m2, double epsilon2);

    int one_connected_integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    int two_connected_integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);

    double one_disconnected (double Q, double Delta);
    double two_disconnected (double Q, double Delta);

    double color_fluctuations (double Q, double Delta);
    double color_fluctuations (CubatureConfig* c_config, IntegrationConfig* i_config);

    double hotspot_fluctuations (double Q, double Delta);

    double dsigmadt (double Q, double Delta);
} }


namespace Incoherent { namespace Sampled
{
    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigmadt_single_event (double Q, double Delta, const HotspotNucleus& nucleus);
    double dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config);
} }