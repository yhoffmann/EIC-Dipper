#pragma once


#include <cuba.h>
#include <tuple>
#include "IntegrationRoutines.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace Coherent {
    double A_integrand_function(double b1, double b2, double r1, double r2, double Q, double Delta);

    int integrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
    double dsigma_dt(double Q, double Delta);
    double dsigma_dt(CubaConfig* cuba_config, IntegrationConfig* integration_config);

    int integrand_cubature(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigma_dt_cubature(double Q, double Delta);
    double dsigma_dt_cubature(CubatureConfig* cuba_config, IntegrationConfig* integration_config);
}


namespace Coherent { namespace Demirci
{
    double Psi0 (double Delta, double m_sqr);
    double Psi (double r, double phi, double lambda, double Delta, double m_sqr);
    double Z (double r, double phi, double lambda, double Delta, double m_sqr, double epsilon);

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    double dsigma_dt (double Q, double Delta);
    double dsigma_dt (CubatureConfig* c_config, IntegrationConfig* i_config);
} }


namespace Coherent { namespace GeometryAverage
{
    double A_integrand_function(double b1, double b2, double r1, double r2, double Q, double Delta, const HotspotNucleus* nucleus);

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    std::tuple<double,double> A (double Q, double Delta, const HotspotNucleus& nucleus);
    std::tuple<double,double> A (CubatureConfig* c_config, IntegrationConfig* i_config);
} }