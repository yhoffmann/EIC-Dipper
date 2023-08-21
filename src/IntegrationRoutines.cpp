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


#include "../include/IntegrationRoutines.h"
#include "../include/constants.h"
#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_sf.h>


namespace IntegrationRoutines {
    double cuba_integrate (integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config) {
        double ret;

        int num_of_dim = cuba_config->num_of_dims;
        int num_of_integrals = 1;
        int num_of_points = 1;
        double epsrel = cuba_config->epsrel;
        double epsabs = cuba_config->epsabs;
        int flags1 = cuba_config->flags1;
        int flags2 = cuba_config->flags2;
        luint mineval = cuba_config->mineval;
        luint maxeval = cuba_config->maxeval;
        int seed = cuba_config->seed;

        cubareal value[num_of_integrals]{0.0}, error[num_of_integrals]{0.0}, probability[num_of_integrals]{0.0};

        int num_of_regions(0), num_of_evals(0), error_status(0);

        if (cuba_config->integrator==Integrator::Cuhre) {
            Cuhre(num_of_dim, num_of_integrals,
                integrand, integration_config, num_of_points,
                epsrel, epsabs,
                flags1 | flags2,
                mineval, maxeval,
                KEY, NULL, NULL,
                &num_of_regions,&num_of_evals, &error_status,
                value, error, probability
            );
        } else if (cuba_config->integrator==Integrator::Suave) {
            Suave(num_of_dim, num_of_integrals,
                integrand, integration_config, num_of_points,
                epsrel, epsabs,
                flags1 | flags2, seed,
                mineval, maxeval,
                cuba_config->n_new, cuba_config->n_min,
                KEY, NULL, NULL,
                &num_of_regions, &num_of_evals, &error_status,
                value, error, probability
            );
        }

        ret = value[0];

        return ret;
    }


    double cuba_integrate_one_bessel (integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config) {
        // for telling cuba_integrate whether to set integration range or not (here, no because range is being set in this function)
        cuba_config->using_bessel_integration = true;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)integration_config->integrand_params;

        bool use_bessel_zeros = (gsl_sf_bessel_zero_J0(1)/A_integrand_params->Delta < 4*B_MAX); // TODO check if this is neccesary or if it even causes inaccuracies
        
        double partial_sum = 0.0;
        double total_sum = 0.0;

        double current_oscillation_value;

        for (int n=0; n<cuba_config->max_oscillations; n++) {
            // Check if partial_sum is small compared to overall total_sum, if yes then stop integration
            if (n%cuba_config->ocillations_per_partial_sum == 0) {
                if (std::abs(partial_sum) < cuba_config->bessel_tolerance*std::abs(total_sum)) {
                    if (cuba_config->progress_monitor) std::cout << A_integrand_params->Q << " " << A_integrand_params->Delta << " " << n << " Converged " << total_sum << std::endl;
                    break;
                } else if ((n>0) && (std::abs(partial_sum)<1e-10) && (std::abs(total_sum)<1e-10)) {
                    if (cuba_config->progress_monitor) std::cout << A_integrand_params->Q << " " << A_integrand_params->Delta << " " << n << " Small Value " << total_sum << std::endl;
                    break;
                } else {
                    partial_sum = 0.0;
                }
            }
            
            if (n == cuba_config->max_oscillations-1) std::cout << "### WARNING SLOW CONVERGENCE ###" << std::endl;

            // Assigning integration range based on Delta to hit the zeros of the J0
            if (n!=0) {
                integration_config->min = integration_config->max;
                integration_config->max = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(n+1)/A_integrand_params->Delta : (n+1)*B_MAX;
            } else {
                integration_config->min = 0;
                integration_config->max = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(1)/A_integrand_params->Delta : B_MAX;
            }

            current_oscillation_value = cuba_integrate(integrand, cuba_config, integration_config);

            // Adding current integration results to overall result
            total_sum += current_oscillation_value;
            partial_sum += current_oscillation_value;

            if (cuba_config->progress_monitor) std::cout << A_integrand_params->Q << " " << A_integrand_params->Delta << " " << A_integrand_params->transverse_or_longitudinal << "\t" << n << "\t(" << integration_config->min << "," << integration_config->max << ")\t" << total_sum << "(+" << current_oscillation_value << ")" << std::endl;
        }

        cuba_config->using_bessel_integration = false;

        return total_sum;
    }
}