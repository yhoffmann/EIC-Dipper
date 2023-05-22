#pragma once


#define SIX_D double,double,double,double,double,double
#define TEN_D double,double,double,double,double,double,double,double,double,double

#include<vector>
#include<cuba.h>
#include<iostream>
#include<gsl/gsl_math.h>
#include"utilities.cpp"
#include"coherent.cpp"
#include"incoherent.cpp"
#include"constants.cpp"
#include"integrand.cpp"
#include"utilities.cpp"


double cuba_integrate (CubaConfig c_config, IntegrandParams i_params) {
    double ret;

    int num_of_regions(0), num_of_evals(0), error_status(0);
    int num_of_integrals = 1;
    int num_of_points = 1;
    double epsrel = c_config.epsrel;
    double epsabs = c_config.epsabs;
    int flags1 = c_config.flags1;
    int flags2 = c_config.flags2;
    luint mineval = c_config.mineval;
    luint maxeval = c_config.maxeval;
    int seed = c_config.seed;
std::cout << "in integrate" << std::endl;
    cubareal value[num_of_integrals]{0.0}, error[num_of_integrals]{0.0}, probability[num_of_integrals]{0};

    int num_of_dim = c_config.num_of_dims;

    if (i_params.min == -999) {
        i_params.min = 0.0;
        i_params.max = std::sqrt(get_b_range_factor()*2.0*BG);
    }

    switch (c_config.integrator) {
        case 'c':
            Cuhre(num_of_dim, num_of_integrals,
                integrand, &i_params, num_of_points,
                epsrel, epsabs,
                flags1 | flags2,
                mineval, maxeval,
                NULL, NULL, NULL,
                &num_of_regions,&num_of_evals, &error_status,
                value, error, probability
            );
            break;
        case 's':
            Suave(num_of_dim, num_of_integrals,
                integrand, &i_params, num_of_points,
                epsrel, epsabs,
                flags1 | flags2, seed,
                mineval, maxeval,
                c_config.n_new, c_config.n_min,
                0, NULL, NULL,
                &num_of_regions,&num_of_evals, &error_status,
                value, error, probability
            );
            break;
    }

    ret = value[0];

    return ret;
}


double cuba_bessel_integrate (CubaConfig c_config, IntegrandParams i_params) {
    double ret;
        

    return ret;
}