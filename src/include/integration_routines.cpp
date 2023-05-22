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

        if (!c_config.using_bessel_integration) {
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


    double cuba_integrate_one_bessel (integrand_t integrand, CubaConfig c_config, IntegrandParams i_params) {
        // for telling cuba_integrate whether to set integration range or not (here, no because range is being set in this function)
        c_config.using_bessel_integration = true;

        bool use_bessel_zeros = (gsl_sf_bessel_zero_J0(1)/i_params.Delta < 4*BMAX);
        
        double partial_sum = 0.0;
        double total_sum = 0.0;

        double current_sum;

        for (int n=0; n<c_config.max_oscillations; n++) {
            // Check if partial_sum is small compared to overall total_sum, if yes then stop integration
            if (n%c_config.ocillations_per_partial_sum == 0) {
                if (std::abs(partial_sum) < c_config.bessel_tolerance*std::abs(total_sum)) {
                    if (c_config.progress_monitor) std::cout << i_params.Q << " " << i_params.Delta << " " << n << " Converged " << total_sum << std::endl;
                    break;
                } else if ((n>0) && (std::abs(partial_sum)<1e-10) && (std::abs(total_sum)<1e-10)) {
                    if (c_config.progress_monitor) std::cout << i_params.Q << " " << i_params.Delta << " " << n << " Small Value " << total_sum << std::endl;
                    break;
                } else {
                    partial_sum = 0.0;
                }
            }
            
            if (n == c_config.max_oscillations-1) std::cout << "### WARNING SLOW CONVERGENCE ###" << std::endl;

            // Assigning integration range based on Delta to hit the zeros of the J0
            if (n!=0) {
                i_params.bmin = i_params.bmax;
                i_params.bmax = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(n+1)/i_params.Delta : (n+1)*brange;
            } else {
                i_params.bmin = 0;
                i_params.bmax = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(1)/i_params.Delta : brange;
            }

            current_oscillation_value = cuba_integrate(integrand, c_config, i_params);

            // Adding current integration results to overall result
            total_sum += current_oscillation_value;
            partial_sum += current_oscillation_value;
            if (c_config.progress_monitor) std::cout << i_params.Q << " " << i_params.Delta << " " << i_params.WhichIntegrand << "\t" << n << "\t(" << i_params.bmin << "," << i_params.bmax << ")\t" << total_sum << "(+" << current_oscillation_value << ")" << std::endl;
        }

        c_config.using_bessel_integration = false;

        return total_sum;
    }
}