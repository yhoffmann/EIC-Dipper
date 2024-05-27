#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "../include/Incoherent.hpp"
#include "../include/Coherent.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/constants.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/utilities.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace Incoherent
{
    double A_integrand_function_factor (double Q)
    {
        return NRPhoton::wave_function_factor(Q) / (4.0*PI);
    }


    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta1, double Delta2)
    {
        return gsl_sf_bessel_J0((b1-bb1)*Delta1 + (b2-bb2)*Delta2) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q) * SaturationModel::dsigma_d2b_sqr_reduced(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, bb1+rb1/2.0, bb2+rb2/2.0, bb1-rb1/2.0, bb2-rb2/2.0);
    }


    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
    {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* p = (AIntegrandParams*)(integration_config->integrand_params);

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

        ff[0] = jacobian * Incoherent::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, p->Q, p->Delta1, p->Delta2);

        return 0;
    }


    double dsigmadt (double Q, double Delta, double phi)
    {
        CubaConfig c_config;
        c_config.progress_monitor = true;

        IntegrationConfig i_config;

        AIntegrandParams params;
        params.Q = Q;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);

        i_config.integrand_params = &params;

        return dsigmadt(&c_config, &i_config);
    }


    double dsigmadt (CubaConfig* c_config, IntegrationConfig* integration_config) // TODO include coherent here, this only calculates incoherent right now
    {
        c_config->num_dims = 8;

        integration_config->min = (double*)alloca(sizeof(double));
        integration_config->max = (double*)alloca(sizeof(double));

        double ret;

        ret = sqr(A_integrand_function_factor(((AIntegrandParams*)(integration_config->integrand_params))->Q)) * IntegrationRoutines::cuba_integrate_one_bessel(Incoherent::integrand, c_config, integration_config);
        ret *= (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);

        return ret;
    }


    int integrand_cubature (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

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

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Incoherent::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, p->Q, p->Delta1, p->Delta2);

        return 0;
    }


    double dsigmadt_cubature (double Q, double Delta, double phi)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = progress_monitor_global;
        c_config.max_eval = 1e6;

        IntegrationConfig i_config;

        AIntegrandParams params;
        params.Q = Q;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);
        params.phi = phi;
        params.is_incoherent = true;

        i_config.integrand_params = &params;

        return dsigmadt_cubature(&c_config, &i_config);
    }


    double dsigmadt_cubature (CubatureConfig* c_config, IntegrationConfig* i_config)
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

        double ret = sqr(A_integrand_function_factor(((AIntegrandParams*)(i_config->integrand_params))->Q)) * IntegrationRoutines::cubature_integrate_zeros(Incoherent::integrand_cubature, c_config, i_config, &gsl_sf_bessel_zero_J0);
        ret *= (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);

        return ret;
    }
}


namespace Incoherent { namespace Demirci
{
    double one_connected_factor (double Q)
    {
        return sqr( NRPhoton::wave_function_factor(Q)/(4.0*PI) ) / (16.0*PI) * sqr(g2mu02) * (sqr(Nc)-1.0) / (2.0*sqr(PI*Nc)) * 3.0;
    }


    double two_connected_factor (double Q)
    {
        return one_connected_factor(Q) * 2.0;
    }


    double K_integrand_function (double k, double phik, double kb, double phikb, double A, double Delta, double m2, double epsilon2)
    {
        return exp( -A*(k*k+kb*kb+2.0*k*kb*cos(phik-phikb)) ) / ( (sqr(k*k+Delta*Delta/4.0+m2)-sqr(k*Delta*cos(phik))) * (sqr(kb*kb+Delta*Delta/4.0+m2)-sqr(kb*Delta*cos(phikb))) ) * ( 1.0/(epsilon2+Delta*Delta/4.0)-1.0/(epsilon2+k*k) ) * ( 1.0/(epsilon2+Delta*Delta/4.0)-1.0/(epsilon2+kb*kb) );
    }


    int one_connected_integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);
        ff[0] = xx[0]*xx[2]*K_integrand_function(xx[0], xx[1], xx[2], xx[3], rH_sqr, p->Delta1, m*m, sqr(NRPhoton::epsilon(p->Q)) );

        return 0;
    }


    int two_connected_integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);
        ff[0] = xx[0]*xx[2]*K_integrand_function(xx[0], xx[1], xx[2], xx[3], R_sqr+rH_sqr, p->Delta1, m*m, sqr(NRPhoton::epsilon(p->Q)) );

        return 0;
    }


    double one_disconnected (double Q, double Delta)
    {
        return Coherent::Demirci::dsigmadt(Q, Delta) / exp(-sqr(Delta)*(rH_sqr+(3.0-1.0)/3.0*R_sqr)) / 3.0 * ( exp(-sqr(Delta)*rH_sqr)-exp(-sqr(Delta)*(rH_sqr+(3.0-1.0)/3.0*R_sqr)) );
    }


    double two_disconnected (double Q, double Delta)
    {
        return Coherent::Demirci::dsigmadt(Q, Delta) / exp(-sqr(Delta)*(rH_sqr+(3.0-1.0)/3.0*R_sqr)) / 3.0 * 2.0 * ( exp(-sqr(Delta)*(rH_sqr+R_sqr))-exp(-sqr(Delta)*(rH_sqr+(3.0-1.0)/3.0*R_sqr)) );
    }


    double color_fluctuations (double Q, double Delta)
    {
        AIntegrandParams p;
        p.Q = Q;
        p.Delta1 = Delta;

        CubatureConfig c_config;
        c_config.progress_monitor = progress_monitor_global;

        IntegrationConfig i_config;

        i_config.integrand_params = &p;

        c_config.max_eval = 1e8;

        return color_fluctuations(&c_config, &i_config);
    }


    double color_fluctuations (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 4;

        double* min = (double*)alloca(4*sizeof(double));
        double* max = (double*)alloca(4*sizeof(double));

        min[0] = 0.0;
        max[0] = 40.0;

        min[1] = 0.0;
        max[1] = 2.0*PI;

        min[2] = 0.0;
        max[2] = 40.0;

        min[3] = 0.0;
        max[3] = 2.0*PI;

        i_config->min = min;
        i_config->max = max;

        double one_connected_result = one_connected_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * IntegrationRoutines::cubature_integrate(one_connected_integrand, c_config, i_config);
        double two_connected_result = two_connected_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * IntegrationRoutines::cubature_integrate(two_connected_integrand, c_config, i_config);

        return (one_connected_result+two_connected_result) * sqr(GeVm1_to_fm) * fm2_to_nb;
    }

    
    double hotspot_fluctuations (double Q, double Delta)
    {
        return one_disconnected(Q, Delta) + two_disconnected(Q, Delta);
    }


    double dsigmadt (double Q, double Delta)
    {
        return color_fluctuations(Q, Delta) + hotspot_fluctuations(Q, Delta);
    }
} }


namespace Incoherent { namespace Sampled
{
    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta1, double Delta2, const HotspotNucleus* h_nucleus)
    {
        return gsl_sf_bessel_J0((b1-bb1)*Delta1 + (b2-bb2)*Delta2) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q) * SaturationModel::Sampled::dsigma_d2b_sqr_reduced(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, bb1+rb1/2.0, bb2+rb2/2.0, bb1-rb1/2.0, bb2-rb2/2.0, h_nucleus);
    }

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

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

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Incoherent::Sampled::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, p->Q, p->Delta1, p->Delta2, p->h_nucleus);

        return 0;
    }


    double dsigmadt_single_event (double Q, double Delta, double phi, const HotspotNucleus& h_nucleus)
    {
        CubatureConfig c_config;
        c_config.max_eval = 1e7;
        c_config.progress_monitor = progress_monitor_global;

        IntegrationConfig i_config;
        AIntegrandParams params;
        
        i_config.integrand_params = &params;
        
        params.Q = Q;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);
        params.phi = phi;
        params.is_incoherent = true;
        params.h_nucleus = &h_nucleus;

        return dsigmadt_single_event(&c_config, &i_config);
    }


    double dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config)
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

        double factor = sqr(Incoherent::A_integrand_function_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q )) * sqr(GeVm1_to_fm) * fm2_to_nb / (16.0*PI);

        return factor*IntegrationRoutines::cubature_integrate_zeros(Incoherent::Sampled::integrand, c_config, i_config, &gsl_sf_bessel_zero_J0);
    }
} }