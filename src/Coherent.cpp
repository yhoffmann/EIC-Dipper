#include "../include/Coherent.hpp"
#include <math.h>
#include "../include/constants.hpp"
#include "../include/IntegrationRoutines.hpp"
#include "../include/GBWModel.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/SaturationModel.hpp"
#include <gsl/gsl_math.h>
#include <iostream>


namespace Coherent {
    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double z, double Delta)
    {
        return 1.0/(8.0*PI*PI) * NRPhoton::wave_function(r1, r2, Q, z) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta ) * SaturationModel::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0);
    }


    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
    {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(integration_config->integrand_params);

        double bmin = integration_config->min[0];
        double bmax = integration_config->max[0];
        //double phibmin = 0.0;
        //double phibmax = 2.0*M_PI;

        //double rmin = 0.0;
        //double rmax = R_MAX;
        //double phirmin = 0.0;
        //double phirmax = 2.0*M_PI;

        double b = bmin + (bmax-bmin)*xx[0];
        double phib = 2.0*PI*xx[1];
        double r = R_MAX*xx[2];
        double phir = 2.0*PI*xx[3];

        double b1 = b*cos(phib);
        double b2 = b*sin(phib);
        double r1 = r*cos(phir);
        double r2 = r*sin(phir);

        double jacobian = b*r*(bmax-bmin)*R_MAX*4.0*PI*PI;

        ff[0] = jacobian * Coherent::A_integrand_function(b1, b2, r1, r2, A_integrand_params->Q, A_integrand_params->z, A_integrand_params->Delta);

        return 0;
    }


    double dsigma_dt (double Q, double Delta)
    {
        CubaConfig c_config;
        c_config.progress_monitor = true;

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

        ret = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand, c_config, integration_config);
        ret = sqr(ret*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }


    int integrand_cubature (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* params = ((AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params));

        ff[0] = xx[0]*xx[2]*Coherent::A_integrand_function(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), params->Q, params->z, params->Delta);

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

        double ret;

        ret = IntegrationRoutines::cubature_integrate_one_bessel(Coherent::integrand_cubature, c_config, i_config);
        ret = sqr(ret*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }
}