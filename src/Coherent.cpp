#include "../include/Coherent.hpp"
#include "../include/constants.hpp"
#include "../include/IntegrationRoutines.hpp"
#include "../include/GBWModel.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/SaturationModel.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include <math.h>
#include <gsl/gsl_math.h>


namespace Coherent
{
    double A_integrand_function_factor (double Q)
    {
        return  NRPhoton::wave_function_factor(Q) / (4.0*PI);
    }


    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta1, double Delta2)
    {
        return gsl_sf_bessel_J0(b1*Delta1 + b2*Delta2) * NRPhoton::wave_function(r1, r2, Q) * NRPhoton::wave_function(rb1, rb2, Q) * SaturationModel::dsigma_d2b(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5) * SaturationModel::dsigma_d2b(bb1+rb1*0.5, bb2+rb2*0.5, bb1-rb1*0.5, bb2-rb2*0.5);
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

        // double b1 = xx[0]*cos(xx[1]);
        // double b2 = xx[0]*sin(xx[1]);
        // double r1 = xx[2]*cos(xx[3]);
        // double r2 = xx[2]*sin(xx[3]);

        // double bb1 = xx[4]*cos(xx[5]);
        // double bb2 = xx[4]*sin(xx[5]);
        // double rb1 = xx[6]*cos(xx[7]);
        // double rb2 = xx[6]*sin(xx[7]);

        ff[0] = xx[0]*xx[2]*xx[4]*xx[6] * Coherent::A_integrand_function(b1, b2, r1, r2, bb1, bb2, rb1, rb2, p->Q, p->Delta1, p->Delta2);

        return 0;
    }


    double dsigmadt (double Q, double Delta, double phi)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = g_monitor_progress;
        c_config.abs_err = 1.0e-7;
        c_config.rel_err = 1.0e-12;
        c_config.max_eval = 1e5;

        IntegrationConfig i_config;
        AIntegrandParams params;
        params.Q = Q;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);
        params.phi = phi;
        params.is_incoherent = true;

        i_config.integrand_params = &params;

        return dsigmadt(&c_config, &i_config);
    }


    double dsigmadt (CubatureConfig* c_config, IntegrationConfig* i_config)
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

        double ret = sqr(A_integrand_function_factor(((AIntegrandParams*)(i_config->integrand_params))->Q)) * IntegrationRoutines::cubature_integrate_zeros(Coherent::integrand, c_config, i_config, &gsl_sf_bessel_zero_J0);
        ret *= (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);

        return ret;
    }


    double A_real (double b1, double b2, double r1, double r2, double Q, double Delta1, double Delta2)
    {
        return NRPhoton::wave_function(r1, r2, Q) * sin(b1*Delta1 + b2*Delta2) * SaturationModel/*::HotspotAverage*/::dsigma_d2b(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5);
    }

    double A_imag (double b1, double b2, double r1, double r2, double Q, double Delta1, double Delta2)
    {
        return NRPhoton::wave_function(r1, r2, Q) * cos(b1*Delta1 + b2*Delta2) * SaturationModel/*::HotspotAverage*/::dsigma_d2b(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5);
    }

    int integrand_real (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        ff[0] = xx[0] * xx[2] * Coherent::A_real(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta1, p->Delta2);

        return 0;
    }

    int integrand_imag (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        ff[0] = xx[0] * xx[2] * Coherent::A_imag(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta1, p->Delta2);

        return 0;
    }

    double dsigmadt_test (double Q, double Delta, double phi)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = g_monitor_progress;
        c_config.abs_err = 1.0e-7;
        c_config.rel_err = 1.0e-12;
        c_config.max_eval = 1e6;

        IntegrationConfig i_config;
        AIntegrandParams params;
        params.Q = Q;
        params.Delta = Delta;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);
        params.phi = phi;
        params.is_incoherent = false;

        i_config.integrand_params = &params;

        c_config.num_dims = 4;

        i_config.min = (double*)alloca(4*sizeof(double));
        i_config.max = (double*)alloca(4*sizeof(double));

        i_config.min[1] = 0.0;
        i_config.max[1] = 2.0*PI;

        i_config.min[2] = 0.0;
        i_config.max[2] = R_MAX;

        i_config.min[3] = 0.0;
        i_config.max[3] = 2.0*PI;

        double ret = sqr(IntegrationRoutines::cubature_integrate_zeros(Coherent::integrand_real, &c_config, &i_config, &sin_zeros));
        i_config.max[2] = R_MAX;
        ret += sqr(IntegrationRoutines::cubature_integrate_zeros(Coherent::integrand_imag, &c_config, &i_config, &cos_zeros));

        ret *= sqr(A_integrand_function_factor(Q)) * (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);

        return ret;
    }

    /*namespace Test
    {
        double A_real (double bII, double bT, double r1, double r2, double Q, double Delta, double dII, double dT)
        {
            double b1 = bII*dII - bT*dT;
            double b2 = bII*dT + bT*dII;
            return NRPhoton::wave_function(r1, r2, Q) * sin(bII*Delta) * SaturationModel::dsigma_d2b(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5);
        }

        double A_imag (double bII, double bT, double r1, double r2, double Q, double Delta, double dII, double dT)
        {
            double b1 = bII*dII - bT*dT;
            double b2 = bII*dT + bT*dII;
            return NRPhoton::wave_function(r1, r2, Q) * cos(bII*Delta) * SaturationModel::dsigma_d2b(b1+r1*0.5, b2+r2*0.5, b1-r1*0.5, b2-r2*0.5);
        }

        int integrand_real (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
        {
            AIntegrandParams* p = (AIntegrandParams*)userdata;

            ff[0] = xx[2] * Coherent::Test::A_real(xx[0], xx[1], xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta, p->Delta1, p->Delta2);

            return 0;
        }

        int integrand_imag (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
        {
            AIntegrandParams* p = (AIntegrandParams*)userdata;

            ff[0] = xx[2] * Coherent::Test::A_imag(xx[0], xx[1], xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta, p->Delta1, p->Delta2);

            return 0;
        }

        double dsigmadt_test (double Q, double Delta, double phi)
        {
            CubatureConfig c_config;
            c_config.progress_monitor = g_monitor_progress;
            c_config.abs_err = 1.0e-1;
            c_config.rel_err = 1.0e-1;
            c_config.max_eval = 1e5;

            IntegrationConfig i_config;
            AIntegrandParams params;
            params.Q = Q;
            params.Delta = Delta;
            params.Delta1 = cos(phi);
            params.Delta2 = sin(phi);
            params.phi = phi;
            params.is_incoherent = false;
            params.is_cartesian = true;

            i_config.integrand_params = &params;

            c_config.num_dims = 4;

            i_config.min = (double*)alloca(4*sizeof(double));
            i_config.max = (double*)alloca(4*sizeof(double));

            i_config.min[1] = -B_MAX;
            i_config.max[1] = B_MAX;

            i_config.min[2] = 0.0;
            i_config.max[2] = R_MAX;

            i_config.min[3] = 0.0;
            i_config.max[3] = 2.0*PI;

            double ret = sqr(IntegrationRoutines::cubature_integrate_zeros(Coherent::Test::integrand_real, &c_config, &i_config, &sin_zeros));

            i_config.min[1] = -B_MAX;
            i_config.max[1] = B_MAX;

            i_config.min[2] = 0.0;
            i_config.max[2] = R_MAX;

            i_config.min[3] = 0.0;
            i_config.max[3] = 2.0*PI;
            // i_config.min[2] = -R_MAX;
            // i_config.max[2] = R_MAX;

            // i_config.min[3] = -R_MAX;
            // i_config.max[3] = R_MAX;

            ret += sqr(IntegrationRoutines::cubature_integrate_zeros(Coherent::Test::integrand_imag, &c_config, &i_config, &cos_zeros));

            ret *= sqr(A_integrand_function_factor(Q)) * (GeVm1_to_fm*GeVm1_to_fm*fm2_to_nb)/(16.0*PI);

            return ret;
        }
    }*/
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
        return r * bessel_K_safe(0, epsilon*r) * ( 2.0*cos(0.5*Delta*r*cos(phi))*Psi0(Delta, m_sqr) - Psi(r, phi, lambda, Delta, m_sqr) );
    }


    double Z_integrand_function_factor (double Q, double Delta)
    {
        return NRPhoton::wave_function_factor(Q) / (4.0*PI) * g_g2mu02 * NH * CF / PI * exp(-RC_sqr*sqr(Delta)/2.0);
    }


    double Z_integrand_function(double r, double phi, double lambda, double Q, double Delta)
    {
        return  Z(r, phi, lambda, Delta, sqr(m), NRPhoton::epsilon(Q));
    }


    int integrand (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* params = (AIntegrandParams*)userdata;

        ff[0] = Z_integrand_function(xx[0], xx[1], xx[2], params->Q, params->Delta1);

        return 0;
    }


    double dsigmadt (double Q, double Delta)
    {
        CubatureConfig c_config;
        c_config.progress_monitor = false;
        c_config.rel_err = 1.0e-14;
        c_config.max_eval = 2e6;

        IntegrationConfig i_config;

        AIntegrandParams params;
        params.Q = Q;
        params.Delta1 = Delta;

        i_config.integrand_params = &params;

        return dsigmadt(&c_config, &i_config);
    }


    double dsigmadt (CubatureConfig* c_config, IntegrationConfig* i_config)
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

        AIntegrandParams* params = (AIntegrandParams*)(i_config->integrand_params);

        double ret = Z_integrand_function_factor(params->Q, params->Delta1) * IntegrationRoutines::cubature_integrate(Demirci::integrand, c_config, i_config);
        ret = sqr(ret*GeVm1_to_fm)*fm2_to_nb/(16.0*PI);

        return ret;
    }
} }


namespace Coherent { namespace Sampled
{
    double A_real (double b1, double b2, double r1, double r2, double Q, double Delta1, double Delta2, const HotspotNucleus* nucleus)
    {
        return NRPhoton::wave_function(r1, r2, Q) * sin(b1*Delta1 + b2*Delta2) * SaturationModel::Sampled::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, nucleus);
    }


    double A_imag (double b1, double b2, double r1, double r2, double Q, double Delta1, double Delta2, const HotspotNucleus* nucleus)
    {
        return NRPhoton::wave_function(r1, r2, Q) * cos(b1*Delta1 + b2*Delta2) * SaturationModel::Sampled::dsigma_d2b(b1+r1/2.0, b2+r2/2.0, b1-r1/2.0, b2-r2/2.0, nucleus);
    }


    int integrand_real (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        ff[0] = xx[0]*xx[2]*Coherent::Sampled::A_real(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta1, p->Delta2, p->h_nucleus);

        return 0;
    }

    int integrand_imag (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        AIntegrandParams* p = (AIntegrandParams*)userdata;

        ff[0] = xx[0]*xx[2]*Coherent::Sampled::A_imag(xx[0]*cos(xx[1]), xx[0]*sin(xx[1]), xx[2]*cos(xx[3]), xx[2]*sin(xx[3]), p->Q, p->Delta1, p->Delta2, p->h_nucleus);

        return 0;
    }


    std::tuple<double,double> sqrt_dsigmadt_single_event (double Q, double Delta, double phi, const HotspotNucleus& h_nucleus)
    {
        CubatureConfig c_config;
        c_config.max_eval = 2e6;
#ifdef _TEST
        c_config.max_eval = 1e3;
#endif
        c_config.progress_monitor = g_monitor_progress;

        IntegrationConfig i_config;
        AIntegrandParams params;

        i_config.integrand_params = &params;

        params.Q = Q;
        params.Delta = Delta;
        params.Delta1 = Delta*cos(phi);
        params.Delta2 = Delta*sin(phi);
        params.phi = phi;
        params.h_nucleus = &h_nucleus;

        return sqrt_dsigmadt_single_event(&c_config, &i_config);
    }


    std::tuple<double,double> sqrt_dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config)
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

        double factor = Coherent::A_integrand_function_factor( ((AIntegrandParams*)(i_config->integrand_params))->Q ) * GeVm1_to_fm * std::sqrt(fm2_to_nb) / (4.0*M_SQRTPI);

        double real = IntegrationRoutines::cubature_integrate_zeros(Coherent::Sampled::integrand_real, c_config, i_config, &sin_zeros);
        i_config->max[2] = R_MAX;
        double imag = IntegrationRoutines::cubature_integrate_zeros(Coherent::Sampled::integrand_imag, c_config, i_config, &cos_zeros);

        return {factor*real, factor*imag};
    }
} }
