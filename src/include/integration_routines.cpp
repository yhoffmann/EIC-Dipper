#pragma once

#define MINEVAL 100000
#define MAXEVAL 1000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 100.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0


namespace Routines {

    double cuba_integrate (integrand_t integrand, CubaConfig c_config, IntegrandParams i_params) {
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

        cubareal value[num_of_integrals]{0.0}, error[num_of_integrals]{0.0}, probability[num_of_integrals]{0};

        int num_of_dim = c_config.num_of_dims;

        if (i_params.min == -999) {
            i_params.min = 0.0;
            i_params.max = std::sqrt(get_b_range_factor()*2.0*BG);
        }

        if (c_config.integrator == 'c') {
            Cuhre(num_of_dim, num_of_integrals,
                integrand, &i_params, num_of_points,
                epsrel, epsabs,
                flags1 | flags2,
                mineval, maxeval,
                NULL, NULL, NULL,
                &num_of_regions,&num_of_evals, &error_status,
                value, error, probability
            );
        } else if (c_config.integrator == 's') {
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
        }

        ret = value[0];

        return ret;
    }


    double cuba_bessel_integrate (integrand_t integrand, CubaConfig c_config, IntegrandParams i_params) {
        double ret;
        
        

        return ret;
    }
}