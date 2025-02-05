#include "../include/Incoherent.hpp"
#include "../include/Coherent.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/constants.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/utilities.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>


namespace Incoherent
{
    double A_integrand_function_factor (double Q)
    {
        return sqr(Coherent::A_integrand_function_factor(Q));
    }


    double A_real (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta)
    {
        return gsl_sf_bessel_J0(std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q) * SaturationModel::dsigma_d2b_sqr_reduced(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5, bb1+rb1*0.5, bb2+rb2*0.5, bb1-rb1*0.5, bb2-rb2*0.5);
    }


    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        double db1 = xx[0]*cos(xx[1]);
        double db2 = xx[0]*sin(xx[1]);

        double B1 = xx[4]*cos(xx[5]);
        double B2 = xx[4]*sin(xx[5]);

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Incoherent::A_real(B1+db1*0.5, B2+db2*0.5, xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), B1-db1*0.5, B2-db2*0.5, xx[6]*cos(xx[7]), xx[6]*sin(xx[7]), p->Q, p->Delta);

        return 0;
    }


    double dsigmadt (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = g_monitor_progress;
        c_config.max_eval = 1e5;

        IntegrationConfig i_config;

        AIntegrandParams params;
        params.Q = Q;
        params.Delta = Delta;
        params.is_incoherent = true;

        i_config.integrand_params = &params;

        return dsigmadt(&c_config, &i_config);
    }


    void reset_integration_range (double min[], double max[])
    {
        min[1] = 0.0;
        max[1] = 2.0*M_PI;

        min[2] = 0.0;
        max[2] = R_MAX;

        min[3] = 0.0;
        max[3] = 2.0*M_PI;

        min[4] = 0.0;
        max[4] = B_MAX;

        min[5] = 0.0;
        max[5] = 2.0*M_PI;

        min[6] = 0.0;
        max[6] = R_MAX;

        min[7] = 0.0;
        max[7] = 2.0*M_PI;
    }

    double dsigmadt (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 8;

        i_config->min = (double*)alloca(8*sizeof(double));
        i_config->max = (double*)alloca(8*sizeof(double));

        reset_integration_range(i_config->min, i_config->max);
        double ret = Incoherent::A_integrand_function_factor(((AIntegrandParams*)(i_config->integrand_params))->Q) * IntegrationRoutines::cubature_integrate_zeros(Incoherent::integrand, c_config, i_config, &gsl_sf_bessel_zero_J0);

        ret *= (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);
        
        return ret;
    }
}


namespace Incoherent { namespace Demirci
{
    double one_connected_factor (double Q)
    {
        return Incoherent::A_integrand_function_factor(Q) / (16.0*PI) * sqr(g_g2mu02) * (sqr(Nc)-1.0) / (2.0*sqr(PI*Nc)) * 3.0;
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
        AIntegrandParams* p = (AIntegrandParams*)userdata;
        ff[0] = xx[0]*xx[2]*K_integrand_function(xx[0], xx[1], xx[2], xx[3], rH_sqr, p->Delta, m*m, sqr(NRPhoton::epsilon(p->Q)) );

        return 0;
    }


    int two_connected_integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;
        ff[0] = xx[0]*xx[2]*K_integrand_function(xx[0], xx[1], xx[2], xx[3], R_sqr+rH_sqr, p->Delta, m*m, sqr(NRPhoton::epsilon(p->Q)) );

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
        p.Delta = Delta;

        CubatureConfig c_config;
        c_config.progress_monitor = g_monitor_progress;

        IntegrationConfig i_config;

        i_config.integrand_params = &p;

        c_config.max_eval = 1e7;

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
    double A (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta, const HotspotNucleus* h_nucleus)
    {
        return gsl_sf_bessel_J0(std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q) * SaturationModel::Sampled::dsigma_d2b_sqr_reduced(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5, bb1+rb1*0.5, bb2+rb2*0.5, bb1-rb1*0.5, bb2-rb2*0.5, h_nucleus);
    }

    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        double db1 = xx[0]*cos(xx[1]);
        double db2 = xx[0]*sin(xx[1]);

        double B1 = xx[4]*cos(xx[5]);
        double B2 = xx[4]*sin(xx[5]);

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Incoherent::Sampled::A(B1+db1*0.5, B2+db2*0.5, xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), B1-db1*0.5, B2-db2*0.5, xx[6]*cos(xx[7]), xx[6]*sin(xx[7]), p->Q, p->Delta, p->h_nucleus);

        return 0;
    }


    double dsigmadt_single_event (double Q, double Delta, const HotspotNucleus& h_nucleus)
    {
        CubatureConfig c_config;
        c_config.max_eval = 5e7;
#ifdef _TEST
        c_config.max_eval = 1e3;
#endif
        c_config.progress_monitor = g_monitor_progress;

        IntegrationConfig i_config;
        AIntegrandParams params;

        i_config.integrand_params = &params;
        
        params.Q = Q;
        params.Delta = Delta;
        params.is_incoherent = true;
        params.h_nucleus = &h_nucleus;

        return dsigmadt_single_event(&c_config, &i_config);
    }


    double dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 8;

        i_config->min = (double*)alloca(8*sizeof(double));
        i_config->max = (double*)alloca(8*sizeof(double));

        double ret = Incoherent::A_integrand_function_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * sqr(GeVm1_to_fm) * fm2_to_nb / (16.0*PI);
        
        reset_integration_range(i_config->min, i_config->max);
        ret *= IntegrationRoutines::cubature_integrate_zeros(Incoherent::Sampled::integrand, c_config, i_config, &gsl_sf_bessel_zero_J0);

        return ret;
    }
} }


namespace Incoherent { namespace InternalHotspotAvg
{
    // A, including the new hotspot averaged dipole cross sections directly
    double A_tilde (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta)
    {
        return gsl_sf_bessel_J0(std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q);
    }

    // < <sigma sigmabar>c >h
    double A_tilde_ssbch (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta)
    {
        return A_tilde(b1, b2, r1, r2, bb1, bb2, rb1, rb2, Q, Delta) * SaturationModel::InternalHotspotAvg::ssbch(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5, bb1+rb1*0.5, bb2+rb2*0.5, bb1-rb1*0.5, bb2-rb2*0.5);
    }

    // < <sigma sigmabar>c >h
    int integrand_ssbch (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        double db1 = xx[0]*cos(xx[1]);
        double db2 = xx[0]*sin(xx[1]);

        double B1 = xx[4]*cos(xx[5]);
        double B2 = xx[4]*sin(xx[5]);

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * A_tilde_ssbch(B1+db1*0.5, B2+db2*0.5, xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), B1-db1*0.5, B2-db2*0.5, xx[6]*cos(xx[7]), xx[6]*sin(xx[7]), p->Q, p->Delta);

        return 0;
    }

    // < <sigma sigmabar>c >h
    double dsdt_ssbch (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.max_eval = 5e7;
#ifdef _TEST
        c_config.max_eval = 1e3;
#endif
        c_config.progress_monitor = g_monitor_progress;

        IntegrationConfig i_config;
        AIntegrandParams params;

        i_config.integrand_params = &params;
        
        params.Q = Q;
        params.Delta = Delta;
        params.is_incoherent = true;
        params.h_nucleus = nullptr;

        return dsdt_ssbch(&c_config, &i_config);
    }

    // < <sigma sigmabar>c >h
    double dsdt_ssbch (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 8;

        i_config->min = (double*)alloca(8*sizeof(double));
        i_config->max = (double*)alloca(8*sizeof(double));

        double ret = A_integrand_function_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * sqr(GeVm1_to_fm) * fm2_to_nb / (16.0*PI);
        
        reset_integration_range(i_config->min, i_config->max);
        ret *= IntegrationRoutines::cubature_integrate_zeros(integrand_ssbch, c_config, i_config, &gsl_sf_bessel_zero_J0);

        return ret;
    }

    // < <sigma>c <sigmabar>c >h
    double A_tilde_scsbch (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta)
    {
        return A_tilde(b1, b2, r1, r2, bb1, bb2, rb1, rb2, Q, Delta) * SaturationModel::InternalHotspotAvg::scsbch(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5, bb1+rb1*0.5, bb2+rb2*0.5, bb1-rb1*0.5, bb2-rb2*0.5);
    }

    // < <sigma>c <sigmabar>c >h
    int integrand_scsbch (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        double db1 = xx[0]*cos(xx[1]);
        double db2 = xx[0]*sin(xx[1]);

        double B1 = xx[4]*cos(xx[5]);
        double B2 = xx[4]*sin(xx[5]);

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * A_tilde_scsbch(B1+db1*0.5, B2+db2*0.5, xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), B1-db1*0.5, B2-db2*0.5, xx[6]*cos(xx[7]), xx[6]*sin(xx[7]), p->Q, p->Delta);

        return 0;
    }

    // < <sigma>c <sigmabar>c >h
    double dsdt_scsbch (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.max_eval = 5e7;
#ifdef _TEST
        c_config.max_eval = 1e3;
#endif
        c_config.progress_monitor = g_monitor_progress;

        IntegrationConfig i_config;
        AIntegrandParams params;

        i_config.integrand_params = &params;

        params.Q = Q;
        params.Delta = Delta;
        params.is_incoherent = true;
        params.h_nucleus = nullptr;

        return dsdt_scsbch(&c_config, &i_config);
    }

    // < <sigma>c <sigmabar>c >h
    double dsdt_scsbch (CubatureConfig* c_config, IntegrationConfig* i_config)
    {
        c_config->num_dims = 8;

        i_config->min = (double*)alloca(8*sizeof(double));
        i_config->max = (double*)alloca(8*sizeof(double));

        double ret = A_integrand_function_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * sqr(GeVm1_to_fm) * fm2_to_nb / (16.0*PI);
        
        reset_integration_range(i_config->min, i_config->max);
        ret *= IntegrationRoutines::cubature_integrate_zeros(integrand_scsbch, c_config, i_config, &gsl_sf_bessel_zero_J0);

        return ret;
    }

    // these functions are slow, use Output::dsdt_nucleus_internal_avg() or call the needed contributing functions directly to combine manually later
    double color_fluc (double Q, double Delta)
    {
        return dsdt_ssbch(Q, Delta) - dsdt_scsbch(Q, Delta);
    }

    double hotspot_fluc (double Q, double Delta)
    {
        return dsdt_scsbch(Q, Delta) - Coherent::InternalHotspotAvg::dsdt_sch_sbch(Q, Delta);
    }

    double dsdt (double Q, double Delta)
    {
        return InternalHotspotAvg::dsdt_ssbch(Q, Delta) - Coherent::InternalHotspotAvg::dsdt(Q, Delta);
    }
} }