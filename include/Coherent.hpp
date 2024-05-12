#pragma once


#include <cuba.h>
#include <tuple>
#include "IntegrationRoutines.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace Coherent {
    double A_integrand_function_factor(double Q);
    double A_integrand_function(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    int integrand(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigmadt(double Q, double Delta);
    double dsigmadt(CubatureConfig* cuba_config, IntegrationConfig* integration_config);

    double dsigmadt_test(double Q, double Delta);
}


namespace Coherent { namespace Demirci
{
    double Psi0 (double Delta, double m_sqr);
    double Psi (double r, double phi, double lambda, double Delta, double m_sqr);
    double Z (double r, double phi, double lambda, double Delta, double m_sqr, double epsilon);

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigmadt (double Q, double Delta);
    double dsigmadt (CubatureConfig* c_config, IntegrationConfig* i_config);
} }


namespace Coherent { namespace Sampled
{
    double A_integrand_function(double b1, double b2, double r1, double r2, double Q, double Delta, const HotspotNucleus* nucleus);

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    std::tuple<double,double> sqrt_dsigmadt_single_event (double Q, double Delta, const HotspotNucleus& nucleus);
    std::tuple<double,double> sqrt_dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config);
} }