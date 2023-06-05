#pragma once


#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cuba.h>
#include "utilities.h"
#include <gsl/gsl_sf.h>
#include "constants.h"


struct CubaConfig {
    char integrator = 'c';
    
    int num_of_dims;
    
    double epsrel = 1.0e-6;
    double epsabs = 1.0e-8;
    luint mineval = 2e5;
    luint maxeval = 1e7;
    int flags1 = 0;
    int flags2 = 4;
    int seed = 0;
    int n_new = 1000;
    int n_min = 2;

    bool using_bessel_integration = false;
    double bessel_tolerance = 1e-5;
    int max_oscillations = 1e4;
    int ocillations_per_partial_sum = 3;

    bool progress_monitor = false;
};


struct IntegrandParams {
    bool t_not_l = true;

    double Delta = 0.001;
    double Q = 0.3;
    const double z = 0.5;

    double min = -999;
    double max = std::sqrt(get_b_range_factor()*2.0*BG);
};


namespace IntegrationRoutines {
    double cuba_integrate (integrand_t integrand, CubaConfig c_config, IntegrandParams i_params);

    double cuba_integrate_one_bessel (integrand_t integrand, CubaConfig c_config, IntegrandParams i_params);
}