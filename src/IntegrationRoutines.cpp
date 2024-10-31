// #define MINEVAL 100000
// #define MAXEVAL 1000000

// #define NSTART 1000
// #define NINCREASE 500
// #define NBATCH 1000
// #define GRIDNO 0
// #define STATEFILE NULL
// #define SPIN NULL

// #define NNEW 1000
// #define NMIN 2
// #define FLATNESS 100.

// #define KEY1 47
// #define KEY2 1
// #define KEY3 1
// #define MAXPASS 5
// #define BORDER 0.
// #define MAXCHISQ 10.
// #define MINDEVIATION .25
// #define NGIVEN 0
// #define LDXGIVEN NDIM
// #define NEXTRA 0

// #define KEY 0


#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_sf.h>
#include "../include/IntegrationRoutines.hpp"
#include "../include/constants.hpp"
#include "../external/cubature/cubature.h"


namespace IntegrationRoutines
{
    // double cuba_integrate (integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config)
    // {
    //     double ret;

    //     int num_dim = cuba_config->num_dims;
    //     int num_integrals = 1;
    //     int num_points = 1;
    //     double epsrel = cuba_config->epsrel;
    //     double epsabs = cuba_config->epsabs;
    //     int flags1 = cuba_config->flags1;
    //     int flags2 = cuba_config->flags2;
    //     luint mineval = cuba_config->mineval;
    //     luint maxeval = cuba_config->maxeval;
    //     int seed = cuba_config->seed;

    //     cubareal value[1]{0.0}, error[1]{0.0}, probability[1]{0.0};

    //     int num_regions(0), num_evals(0), error_status(0);

    //     if (cuba_config->integrator==CubaIntegrator::Cuhre) {
    //         Cuhre(num_dim, num_integrals,
    //             integrand, integration_config, num_points,
    //             epsrel, epsabs,
    //             flags1 | flags2,
    //             mineval, maxeval,
    //             KEY, NULL, NULL,
    //             &num_regions, &num_evals, &error_status,
    //             value, error, probability
    //         );
    //     } else if (cuba_config->integrator==CubaIntegrator::Suave) {
    //         Suave(num_dim, num_integrals,
    //             integrand, integration_config, num_points,
    //             epsrel, epsabs,
    //             flags1 | flags2, seed,
    //             mineval, maxeval,
    //             cuba_config->n_new, cuba_config->n_min,
    //             FLATNESS, NULL, NULL,
    //             &num_regions, &num_evals, &error_status,
    //             value, error, probability
    //         );
    //     }

    //     ret = value[0];

    //     return ret;
    // }


    // double cuba_integrate_one_bessel (integrand_t integrand, CubaConfig* cuba_config, IntegrationConfig* integration_config)
    // {
    //     AIntegrandParams* p = (AIntegrandParams*)integration_config->integrand_params;

    //     double Delta = std::sqrt(sqr(p->Delta1) + sqr(p->Delta2));
    //     bool use_bessel_zeros = (gsl_sf_bessel_zero_J0(1)/Delta < 2.0*B_MAX); // TODO check if this is neccesary or if it even causes inaccuracies
        
    //     double partial_sum = 0.0;
    //     double total_sum = 0.0;

    //     double current_oscillation_value;

    //     for (uint n=0; n<cuba_config->max_oscillations; n++)
    //     {
    //         // Check if partial_sum is small compared to overall total_sum, if yes then stop integration
    //         if (n%cuba_config->ocillations_per_partial_sum == 0)
    //         {
    //             if (std::abs(partial_sum) < cuba_config->bessel_tolerance*std::abs(total_sum))
    //             {
    //         #ifndef _QUIET
    //                 if (cuba_config->progress_monitor)
    //                     std::cout << p->Q << " " << Delta << " " << n << " Converged " << total_sum << std::endl;
    //         #endif
    //                 break;
    //             }
    //             else
    //             {
    //                 partial_sum = 0.0;
    //             }
    //         }
    // #ifndef _QUIET
    //         if (n == cuba_config->max_oscillations-1)
    //             std::cout << "### WARNING SLOW CONVERGENCE ###" << std::endl;
    // #endif
    //         // Assigning integration range based on Delta to hit the zeros of the J0
    //         if (n!=0)
    //         {
    //             integration_config->min[0] = integration_config->max[0];
    //             integration_config->max[0] = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(n+1)/Delta : (n+1)*B_MAX/2.0;
    //         }
    //         else
    //         {
    //             integration_config->min[0] = 0.0;
    //             integration_config->max[0] = (use_bessel_zeros) ? gsl_sf_bessel_zero_J0(1)/Delta : B_MAX/2.0;
    //         }

    //         current_oscillation_value = cuba_integrate(integrand, cuba_config, integration_config);

    //         // Adding current integration results to overall result
    //         total_sum += current_oscillation_value;
    //         partial_sum += current_oscillation_value;
    // #ifndef _QUIET
    //         if (cuba_config->progress_monitor)
    //         {
    //             std::cout << (((AIntegrandParams*)integration_config->integrand_params)->is_incoherent ? "inco " : "co ") << p->Q << " " << Delta << " " << p->phi << "\t" << n << "\t(" << integration_config->min[0] << "," << integration_config->max[0] << ")\t" << total_sum << "(+" << current_oscillation_value << ")" << std::endl;
    //         }
    // #endif
    //     }

    //     return total_sum;
    // }


    double cubature_integrate (integrand integrand, CubatureConfig* cubature_config, IntegrationConfig* integration_config)
    {
        double result = 0.0;
        double result_err = 0.0;

        if (cubature_config->integrator==CubatureIntegrator::H)
        {
            hcubature(cubature_config->num_f_dims,
                integrand,
                integration_config->integrand_params,
                cubature_config->num_dims,
                integration_config->min, integration_config->max,
                cubature_config->max_eval, cubature_config->abs_err, cubature_config->rel_err,
                cubature_config->err_norm,
                &result, &result_err
            );
        }
        else if (cubature_config->integrator==CubatureIntegrator::P)
        {
            pcubature(cubature_config->num_f_dims,
                integrand,
                integration_config->integrand_params,
                cubature_config->num_dims,
                integration_config->min, integration_config->max,
                cubature_config->max_eval, cubature_config->abs_err, cubature_config->rel_err,
                cubature_config->err_norm,
                &result, &result_err
            );
        }

        return result;
    }


    double cubature_integrate_zeros (integrand integrand, CubatureConfig* cubature_config, IntegrationConfig* integration_config, double (*zeros)(uint n))
    {
        AIntegrandParams* p = (AIntegrandParams*)integration_config->integrand_params;

        double Delta = p->Delta;//std::sqrt(sqr(p->Delta1) + sqr(p->Delta2));
        bool use_zeros = (zeros(1)/Delta < B_MAX/8.0);
        
        double partial_sum = 0.0;
        double total_sum = 0.0;

        double current_oscillation_value;

        for (uint n=0; n<cubature_config->max_oscillations; n++)
        {
            // Check if partial_sum is small compared to overall total_sum, if yes then stop integration
            if (n%cubature_config->ocillations_per_partial_sum == 0)
            {
                if (std::abs(partial_sum) < cubature_config->bessel_tolerance*std::abs(total_sum) && (n+1) > cubature_config->min_oscillations)
                {
                    if (cubature_config->progress_monitor)
                        std::cout << (((AIntegrandParams*)integration_config->integrand_params)->is_incoherent ? "inco " : "co ") << p->Q << " " << Delta << " " << n << " Converged " << total_sum << std::endl;
                    
                    break;
                }
                else
                {
                    partial_sum = 0.0;
                }
            }

            if (n!=0)
            {
                integration_config->min[0] = (use_zeros) ? zeros(n)/Delta : (n)*B_MAX/8.0;
                integration_config->max[0] = (use_zeros) ? zeros(n+1)/Delta : (n+1)*B_MAX/8.0;
            }
            else // if n == 0
            {
                integration_config->min[0] = 0.0;
                integration_config->max[0] = (use_zeros) ? zeros(1)/Delta : B_MAX/8.0;
            }
            
            integration_config->max[2] = std::max(integration_config->max[2], integration_config->max[0]);
            if (p->is_incoherent)
            {
                integration_config->max[4] = std::max(integration_config->max[4], integration_config->max[0]);
                integration_config->max[6] = std::max(integration_config->max[6], integration_config->max[0]);
            }

            current_oscillation_value = cubature_integrate(integrand, cubature_config, integration_config);

            total_sum += current_oscillation_value;
            partial_sum += current_oscillation_value;
            
            if (cubature_config->progress_monitor)
            {
                std::cout << (((AIntegrandParams*)integration_config->integrand_params)->is_incoherent ? "inco " : "co ") << p->Q << " " << Delta << " " << p->phi << "\t" << n << "\t"<< (use_zeros ? "root" : "") <<"(" << integration_config->min[0] << "," << integration_config->max[0] << ")\t" << total_sum << ((current_oscillation_value > 0.0) ? "(+" : "(") << current_oscillation_value << ")" << std::endl;

                if (n == cubature_config->max_oscillations-1)
                    std::cout << "### WARNING SLOW CONVERGENCE ###" << std::endl;

            }
        }

        return total_sum;
    }
}
