#include "../include/Incoherent.hpp"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "../include/NRPhoton.hpp"
#include "../include/constants.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/utilities.hpp"


namespace Incoherent
{
    double A_integrand_function_simple (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta)
    {
        return 1.0/(16.0*M_PI*M_PI) * gsl_sf_bessel_J0( std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta )/(4.0*PI*PI) * NRPhoton::wave_function( r1, r2, Q, z) * NRPhoton::wave_function(rb1, rb2, Q, z) * SaturationModel::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0) * SaturationModel::dsigma_d2b(bb1+rb1/2.0, bb2+rb2/2.0, bb1-rb1/2.0, bb2-rb2/2.0);
    }

    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta)
    {
        return 1.0/(64.0*M_PI*M_PI*M_PI*M_PI) * gsl_sf_bessel_J0( std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta ) * NRPhoton::wave_function(r1,r2,Q,z) * NRPhoton::wave_function(rb1,rb2,Q,z) * SaturationModel::dsigma_d2b_sqr(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, bb1+rb1/2.0, bb2+rb2/2.0, bb1-rb1/2.0, bb2-rb2/2.0);
    }

    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
    {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(integration_config->integrand_params);

        double dbmin = integration_config->min[0];
        double dbmax = integration_config->max[0];
        //double phidbmin = 0.0;
        //double phidbmax = 2*PI;

        //double rmin = 0.0;
        //double rmax = R_MAX;
        //double phirmin = 0.0;
        //double phirmax = 2.0*PI;

        //double Bmin = 0.0;
        //double Bmax = B_MAX;
        //double phiBmin = 0.0;
        //double phiBmax = 2*PI;

        //double rbmin = 0.0;
        //double rbmax = rmax;
        //double phirbmin = 0.0;
        //double phirbmax = 2.0*PI;

        double db = dbmin + (dbmax-dbmin)*xx[0];
        double phib = 2.0*M_PI*xx[1];
        double r = R_MAX*xx[2];
        double phir = 2.0*M_PI*xx[3];

        double B = B_MAX*xx[4];
        double phiB = 2.0*M_PI*xx[5];
        double rb = R_MAX*xx[6];
        double phirb = 2.0*M_PI*xx[7];

        double db1 = db*cos(phib);
        double db2 = db*sin(phib);
        double r1 = r*cos(phir);
        double r2 = r*sin(phir);

        double B1 = B*cos(phiB);
        double B2 = B*sin(phiB);
        double rb1 = rb*cos(phirb);
        double rb2 = rb*sin(phirb);

        double b1 = B1+db1/2.0;
        double b2 = B2+db2/2.0;
        double bb1 = B1-db1/2.0;
        double bb2 = B2-db2/2.0;

        double jacobian = r*db*B*rb*(dbmax-dbmin)*R_MAX*B_MAX*R_MAX*16.0*PI*PI*PI*PI;

        ff[0] = jacobian * Incoherent::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, A_integrand_params->Q, A_integrand_params->z, A_integrand_params->Delta);

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


    double dsigma_dt (CubaConfig* c_config, IntegrationConfig* integration_config) // TODO include coherent here, this only calculates incoherent right now
    {
        c_config->num_dims = 8;

        integration_config->min = (double*)alloca(sizeof(double));
        integration_config->max = (double*)alloca(sizeof(double));

        double ret;

        ret = IntegrationRoutines::cuba_integrate_one_bessel(Incoherent::integrand, c_config, integration_config);
        ret *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }


    int integrand_cubature (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);

        double db1 = xx[0]*cos(xx[1]);
        double db2 = xx[0]*sin(xx[1]);
        double r1 = xx[2]*cos(xx[3]);
        double r2 = xx[2]*sin(xx[3]);

        double B1 = xx[4]*cos(xx[5]);
        double B2 = xx[4]*sin(xx[5]);
        double rb1 = xx[6]*cos(xx[7]);
        double rb2 = xx[6]*sin(xx[7]);

        double b1 = B1+db1/2.0;
        double b2 = B2+db2/2.0;
        double bb1 = B1-db1/2.0;
        double bb2 = B2-db2/2.0;

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Incoherent::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, A_integrand_params->Q, A_integrand_params->z, A_integrand_params->Delta);

        return 0;
    }


    double dsigma_dt_cubature (double Q, double Delta)
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


    double dsigma_dt_cubature (CubatureConfig* c_config, IntegrationConfig* i_config) // TODO include coherent here, this only calculates incoherent right now
    {
        c_config->num_dims = 8;

        i_config->min = (double*)alloca(8*sizeof(double));
        i_config->max = (double*)alloca(8*sizeof(double));

        i_config->min[1] = 0.0;
        i_config->max[1] = 2.0*M_PI;

        i_config->min[2] = 0.0;
        i_config->max[2] = R_MAX;

        i_config->min[3] = 0.0;
        i_config->max[3] = 2.0*M_PI;

        i_config->min[4] = 0.0;
        i_config->max[4] = B_MAX;

        i_config->min[5] = 0.0;
        i_config->max[5] = 2.0*M_PI;

        i_config->min[6] = 0.0;
        i_config->max[6] = R_MAX;

        i_config->min[7] = 0.0;
        i_config->max[7] = 2.0*M_PI;

        double ret;

        ret = IntegrationRoutines::cubature_integrate_one_bessel(Incoherent::integrand_cubature, c_config, i_config);
        ret *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}