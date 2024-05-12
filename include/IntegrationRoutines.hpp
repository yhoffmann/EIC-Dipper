#pragma once


#include <cuba.h>
#include "../include/utilities.hpp"
#include "../external/cubature/cubature.h"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


enum class CubaIntegrator : unsigned char
{
    Cuhre, Suave
};


enum class CubatureIntegrator : unsigned char
{
    H, P
};


struct CubaConfig
{
    CubaIntegrator integrator = CubaIntegrator::Cuhre;
    
    uint num_dims;
    
    double epsrel = 1.0e-6;
    double epsabs = 1.0e-8;
    luint mineval = 2e5;
    luint maxeval = 1e7;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;
    int n_new = 1000;
    int n_min = 2;

    double bessel_tolerance = 1.0e-5;
    uint max_oscillations = 40;
    uint ocillations_per_partial_sum = 3;

    bool progress_monitor = false;
};


struct CubatureConfig
{
    CubatureIntegrator integrator = CubatureIntegrator::H;

    uint num_dims = 2;
    uint num_f_dims = 1;
    size_t max_eval = 1e8;
    double abs_err = 1.0e-14;
    double rel_err = 1.0e-6;
    error_norm err_norm = ERROR_INDIVIDUAL;

    double bessel_tolerance = 1.0e-4;
    uint min_oscillations = 5;
    uint max_oscillations = 40;
    uint ocillations_per_partial_sum = 3;

    bool progress_monitor = false;
};


struct AIntegrandParams
{
    double Delta = 0.001;
    double Q = std::sqrt(0.1);

    bool is_incoherent = false;

    const HotspotNucleus* h_nucleus = nullptr;
};


struct GIntegrandParams
{
    double x1, x2, y1, y2;
};


struct IntegrationConfig
{
    double* min;
    double* max;

    void* integrand_params; // integrand specific parameters, see above
};


namespace IntegrationRoutines
{
    double cuba_integrate(integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config);
    double cuba_integrate_one_bessel(integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config);

    double cubature_integrate(integrand integrand, CubatureConfig* cubature_config, IntegrationConfig* integration_config);
    double cubature_integrate_one_bessel(integrand integrand, CubatureConfig* cubature_config, IntegrationConfig* integration_config);
}