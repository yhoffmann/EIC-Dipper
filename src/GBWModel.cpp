#include "../include/GBWModel.hpp"
#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include <math.h>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "../include/IntegrationRoutines.hpp"


namespace GBWModel {
    double T_times_sigma0 (double b1, double b2)
    {
        return exp( -( sqr(b1)+sqr(b2) )/(2.0*BG) );
    }


    double Q_s_sqr (double b1, double b2)
    {
        return sqr(Qs0) * T_times_sigma0(b1,b2);
    }


    double G_old (double x1, double x2, double y1, double y2) // not needed
    {
        return -0.25 * Q_s_sqr( (x1+y1)/2.0, (x2+y2)/2.0 ) * ( sqr(x1-y1) + sqr(x2-y2) );
    }


    double Kmod (double x1, double x2, double y1, double y2) // not needed
    {
        return 4.0* ( 1.0-m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))*bessel_K_safe(1,m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))) ) / sqr(m) / (sqr(x1-y1)+sqr(x2-y2));
    }


    double G_mod (double x1, double x2, double y1, double y2)
    {
        return -0.25 * exp( -(sqr(x1)+sqr(x2))/(2*BG) ) * exp( -(sqr(y1)+sqr(y2))/(2*BG) ) * ( sqr(x1-y1) + sqr(x2-y2) ) * Kmod(x1,x2,y1,y2);
    }


    double G_integrand_function (double u, double v, double x1, double x2, double y1, double y2)
    {
        if (u==0 && v==0) return 0.0;
        
        volatile double inverse_divisor = 1.0/(u*v+BG/2.0*(u+v));

        return Nq * CF * g2mu02_demirci * exp(-sqr(m)*(u+v)) * inverse_divisor / (16.0*PI*PI) * ( exp( -0.25 * ( u*(sqr(y1)+sqr(y2)) + v*(sqr(x1)+sqr(x2)) + BG/2.0*(sqr(x1-y1)+sqr(x2-y2)) ) * inverse_divisor ) - 0.5 * exp( -0.25 * (sqr(x1)+sqr(x2)) * (u+v) * inverse_divisor ) - 0.5 * exp( -0.25 * (sqr(y1)+sqr(y2)) * (u+v) * inverse_divisor ) );
    }


    int G_integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
    {
        GIntegrandParams* params = (GIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);

        // integration from 0, so no umin needed
        double umax = 20.0/m;

        double u = umax*xx[0];
        double v = umax*xx[1];

        double jacobian = umax*umax;

        ff[0] = jacobian*G_integrand_function(u,v,params->x1,params->x2,params->y1,params->y2);

        return 0;
    }


    int G_integrand_cubature (unsigned ndim, const double* xx, void* userdata, unsigned fdim, double* ff)
    {
        GIntegrandParams* params = (GIntegrandParams*)(((IntegrationConfig*)userdata)->integrand_params);

        ff[0] = G_integrand_function(xx[0],xx[1],params->x1,params->x2,params->y1,params->y2);

        return 0;
    }


    double G_by_integration (double x1, double x2, double y1, double y2)
    {
        CubatureConfig cubature_config;
        cubature_config.num_dims = 2;
        cubature_config.max_eval = 1e6;

        IntegrationConfig integration_config;
        GIntegrandParams G_integrand_params
        {
            x1,
            x2,
            y1,
            y2
        };
        integration_config.integrand_params = &G_integrand_params;

        //double ret = IntegrationRoutines::cuba_integrate(G_integrand,&cuba_config,&integration_config); // TODO check unit and factors

        integration_config.min = (double*)alloca(2*sizeof(double));
        integration_config.max = (double*)alloca(2*sizeof(double));
        
        integration_config.min[0] = 0.0;
        integration_config.min[1] = 0.0;
        
        integration_config.max[0] = 20.0/sqr(m);
        integration_config.max[1] = integration_config.max[0];

        double ret_cubature = IntegrationRoutines::cubature_integrate(G_integrand_cubature,&cubature_config,&integration_config);

        //std::cout << ret << " " << ret_cubature << "\n";

        return ret_cubature;
    }


    double G_wrapper (double r, double rb, double theta)
    {
        double x1 = r;
        double x2 = 0.0;

        double y1 = rb * cos(theta);
        double y2 = rb * sin(theta);

        double ret = G_by_integration(x1,x2,y1,y2);

        return ret;
    }


    Interpolator3D G_ip;


    double G (double x1, double x2, double y1, double y2)
    {
        double r = std::sqrt(sqr(x1)+sqr(x2));
        double rb = std::sqrt(sqr(y1)+sqr(y2));
        double arg = (r!=0.0 && rb!=0.0) ? (x1*y1+x2*y2)/(r*rb) : 0.0;
        if (arg<-1.0)
            arg = -1.0;
        if (arg>1.0)
            arg = 1.0;
        double theta = acos(arg);
        
        return G_ip.get_interp_value_tricubic(r,rb,theta);
    }
}