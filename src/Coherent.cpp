#include "../include/Coherent.hpp"
#include <math.h>
#include "../include/constants.hpp"
#include "../include/IntegrationRoutines.hpp"
#include "../include/GBWModel.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/SaturationModel.hpp"
#include <gsl/gsl_math.h>
#include <iostream>


namespace Coherent
{
    double A_integrand_function_factor (double Q)
    {
        return 1.0/(4.0*PI) * NRPhoton::wave_function_factor(Q);
    }


    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double Delta)
    {
        return NRPhoton::wave_function(r1, r2, Q) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta ) * SaturationModel::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0);
    }


    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
    {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(integration_config->integrand_params);

        double bmin = integration_config->min[0];
        double bmax = integration_config->max[0];

        double b = bmin + (bmax-bmin)*xx[0];
        double phib = 2.0*PI*xx[1];
        double r = R_MAX*xx[2];
        double phir = 2.0*PI*xx[3];

        double b1 = b*cos(phib);
        double b2 = b*sin(phib);
        double r1 = r*cos(phir);
        double r2 = r*sin(phir);

        double jacobian = b*r*(bmax-bmin)*R_MAX*4.0*PI*PI;

        ff[0] = jacobian * Coherent::A_integrand_function(b1, b2, r1, r2, A_integrand_params->Q, A_integrand_params->Delta);

        return 0;
    }


    double dsigma_dt (double Q, double Delta)
    {
        CubaConfig c_config;
        c_config.progress_monitor = false;

        IntegrationConfig i_config;
        
        AIntegrandParams params;
        params.Q = Q;
        params.Delta = Delta;

        i_config.integrand_params = &params;

        return dsigma_dt(&c_config, &i_config);
    }


    double dsigma_dt (CubaConfig* c_config, IntegrationConfig* integration_config)
    {
        c_config->num_dims = 4;

        integration_config->min = (double*)alloca(sizeof(double));
        integration_config->max = (double*)alloca(sizeof(double));

        double ret;

        ret = A_integrand_function_factor(((AIntegrandParams*)integration_config->integrand_params)->Q) * IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand, c_config, integration_config);
        ret = sqr(ret*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }


    int integrand_cubature (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* params = ((AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params));

        ff[0] = xx[0]*xx[2]*Coherent::A_integrand_function(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), params->Q, params->Delta);

        return 0;
    }


    double dsigma_dt_cubature (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = true;

        IntegrationConfig i_config;
        
        AIntegrandParams params;
        params.Q = Q;
        params.Delta = Delta;

        i_config.integrand_params = &params;

        return dsigma_dt_cubature(&c_config, &i_config);
    }


    double dsigma_dt_cubature (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 4;

        i_config->min = (double*)alloca(4*sizeof(double));
        i_config->max = (double*)alloca(4*sizeof(double));

        i_config->min[1] = 0.0;
        i_config->max[1] = 2.0*M_PI;

        i_config->min[2] = 0.0;
        i_config->max[2] = R_MAX;

        i_config->min[3] = 0.0;
        i_config->max[3] = 2.0*M_PI;

        double ret = A_integrand_function_factor(((AIntegrandParams*)(i_config->integrand_params))->Q) * IntegrationRoutines::cubature_integrate_one_bessel(Coherent::integrand_cubature, c_config, i_config);
        ret = sqr(ret*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }
}


namespace Coherent { namespace Demirci
{
    double Psi0 (double Delta, double m_sqr)
    {
        return 2.0 * atanh( Delta/std::sqrt( sqr(Delta)+4.0*m_sqr ) ) / (Delta * std::sqrt( sqr(Delta)+4.0*m_sqr ));
    }


    double Psi (double r, double phi, double lambda, double Delta, double m_sqr)
    {
        double sqrt = std::sqrt(-sqr(Delta*lambda) + 0.25*sqr(Delta) + m_sqr);

        return cos(lambda*Delta*r*cos(phi)) * r * bessel_K_safe(1, r*sqrt)/sqrt;
    }


    double Z (double r, double phi, double lambda, double Delta, double m_sqr, double epsilon)
    {
        return r * bessel_K_safe(0, epsilon*r) * ( 2.0*cos(0.5*Delta*r*cos(phi))*Psi0(Delta, m_sqr) - Psi(r, phi, lambda, Delta, m_sqr));
    }


    double Z_integrand_function_factor (double Q, double Delta)
    {
        return NRPhoton::wave_function_factor(Q) / (4.0*PI) * g2mu02_demirci * Nq * CF / PI * exp(-RC_sqr*sqr(Delta)/2.0);
    }


    double Z_integrand_function(double r, double phi, double lambda, double Q, double Delta)
    {
        return  Z(r, phi, lambda, Delta, sqr(m), NRPhoton::epsilon(Q));
    }


    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* params = (AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);

        ff[0] = Z_integrand_function(xx[0], xx[1], xx[2], params->Q, params->Delta);

        return 0;
    }


    double dsigma_dt (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = false;
        c_config.rel_err = 1.0e-14;
        c_config.max_eval = 2e7;

        IntegrationConfig i_config;
        
        AIntegrandParams params;
        params.Q = Q;
        params.Delta = Delta;

        i_config.integrand_params = &params;

        return dsigma_dt(&c_config, &i_config);
    }


    double dsigma_dt (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 3;
        
        i_config->min = (double*)alloca(3*sizeof(double));
        i_config->max = (double*)alloca(3*sizeof(double));

        i_config->min[0] = 0.0;
        i_config->max[0] = 40.0;

        i_config->min[1] = 0.0;
        i_config->max[1] = 2.0*PI;

        i_config->min[2] = 0.0;
        i_config->max[2] = 0.5;

        c_config->max_eval = 1e7;

        AIntegrandParams* params = (AIntegrandParams*)(i_config->integrand_params);

        double ret = Z_integrand_function_factor(params->Q, params->Delta) * IntegrationRoutines::cubature_integrate(Demirci::integrand, c_config, i_config);
        ret = sqr(ret*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }
} }


namespace Coherent { namespace GeometryAverage
{
    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double Delta, const Nucleus* nucleus)
    {
        return NRPhoton::wave_function(r1, r2, Q) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta ) * SaturationModel::GeometryAverage::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, nucleus);
    }


    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* params = ((AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params));

        ff[0] = xx[0]*xx[2]*Coherent::GeometryAverage::A_integrand_function(xx[0]*cos(xx[1])-params->b01, xx[0]*sin(xx[1])-params->b02, xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), params->Q, params->Delta, params->nucleus);

        return 0;
    }


    std::tuple<double,double> A (double Q, double Delta, const Nucleus& nucleus)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = false;

        IntegrationConfig i_config;
        AIntegrandParams params;
        
        i_config.integrand_params = &params;
        
        params.Q = Q;
        params.Delta = Delta;
        params.nucleus = &nucleus;

        return A(&c_config, &i_config);
    }


    std::tuple<double,double> A (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 4;

        i_config->min = (double*)alloca(4*sizeof(double));
        i_config->max = (double*)alloca(4*sizeof(double));

        i_config->min[1] = 0.0;
        i_config->max[1] = 2.0*M_PI;

        i_config->min[2] = 0.0;
        i_config->max[2] = R_MAX;

        i_config->min[3] = 0.0;
        i_config->max[3] = 2.0*M_PI;

        return {0.0, A_integrand_function_factor(((AIntegrandParams*)(i_config->integrand_params))->Q) * IntegrationRoutines::cubature_integrate_one_bessel(Coherent::GeometryAverage::integrand, c_config, i_config)};
    }
} }