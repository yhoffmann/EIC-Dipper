#pragma once


#include <cuba.h>
#include <tuple>
#include "IntegrationRoutines.hpp"


typedef integrand integrand_cubature;


namespace Coherent {
    double A_integrand_function(double x1, double x2, double y1, double y2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal);

    int integrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
    std::tuple<double,double> dsigma_dt(double Q, double Delta);
    std::tuple<double,double> dsigma_dt(CubaConfig* cuba_config, IntegrationConfig* integration_config);

    int integrand_cubature(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    std::tuple<double,double> dsigma_dt_cubature(double Q, double Delta);
    std::tuple<double,double> dsigma_dt_cubature(CubatureConfig* cuba_config, IntegrationConfig* integration_config);
}


namespace Coherent { namespace DiluteApprox
{
    double A_integrand_function(double b1, double b2, double r1, double r2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal);

    int integrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
    std::tuple<double,double> dsigma_dt(double Q, double Delta);
    std::tuple<double,double> dsigma_dt(CubaConfig* c_config, IntegrationConfig* integration_config);

    int integrand_cubature(unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff);
    std::tuple<double,double> dsigma_dt_cubature(double Q, double Delta);
    std::tuple<double,double> dsigma_dt_cubature(CubatureConfig cubature_config);
} }