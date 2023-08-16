#include "../include/GBWModel.h"
#include "../include/utilities.h"
#include "../include/constants.h"
#include <math.h>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "../include/IntegrationRoutines.h"


namespace GBWModel {
    double T_times_sigma0 (double b1, double b2) {
        return exp( -( sqr(b1)+sqr(b2) )/(2.0*BG) );
    }

    double Q_s_sqr (double b1, double b2) {
        return sqr(Qs0) * T_times_sigma0(b1,b2);
    }

    double G_old (double x1, double x2, double y1, double y2) { // not needed
        return -0.25 * Q_s_sqr( (x1+y1)/2.0, (x2+y2)/2.0 ) * ( sqr(x1-y1) + sqr(x2-y2) );
    }

    double Kmod (double x1, double x2, double y1, double y2) { // not needed
        return 4.0* ( 1.0-m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))*bessel_K_safe(1,m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))) ) / sqr(m) / (sqr(x1-y1)+sqr(x2-y2));
    }

    double G_mod (double x1, double x2, double y1, double y2) {
        return -0.25 * exp( -(sqr(x1)+sqr(x2))/(2*BG) ) * exp( -(sqr(y1)+sqr(y2))/(2*BG) ) * ( sqr(x1-y1) + sqr(x2-y2) ) * Kmod(x1,x2,y1,y2);
    }

    double G_integrand_function (double u, double v, double x1, double x2, double y1, double y2) {
        if (u==0 && v==0) { u = 1.0e-20; v = u; } // TODO check integration if this returns 0
        
        return 0.25 * (u*v+BG/2.0*(u+v)) * exp(-sqr(m)*(u+v)) * ( exp( -0.25 * ( u*(sqr(y1)+sqr(y2)) + v*(sqr(x1)+sqr(x2)) + BG/2.0*(sqr(x1-y1)+sqr(x2-y2)) ) / ( u*v+BG/2.0*(u+v) ) ) - 0.5 * exp( -0.25 * (sqr(x1)+sqr(x2)) / (u*v/(u+v)+BG/2.0) ) - 0.5 * exp( -0.25 * (sqr(y1)+sqr(y2)) / (u*v/(u+v)+BG/2.0) ) );
    }

    int G_integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata) {
        IntegrationConfig* integration_config = (IntegrationConfig*) userdata;
        
        GIntegrandParams* params = (GIntegrandParams*) integration_config->integrand_params;

        // integration from 0, so no umin needed
        double umax = 40.0/sqr(m);

        // integration from 0, so no vmin needed
        double vmax = umax;

        double u = umax*xx[0];
        double v = vmax*xx[1];

        double jacobian = umax*vmax;

        ff[0] = jacobian*G_integrand_function(u,v,params->x1,params->x2,params->y1,params->y2);

        return 0;
    }

    double G_by_integration (double x1, double x2, double y1, double y2) {
        CubaConfig cuba_config;
        cuba_config.num_of_dims = 2;
        cuba_config.maxeval = 100000;

        IntegrationConfig integration_config;
        
        GIntegrandParams G_integrand_params;
        G_integrand_params.x1 = x1;
        G_integrand_params.x2 = x2;
        G_integrand_params.y1 = y1;
        G_integrand_params.y2 = y2;

        integration_config.integrand_params = &G_integrand_params;

        double ret;

        ret = IntegrationRoutines::cuba_integrate(G_integrand,&cuba_config,&integration_config); // TODO check unit and factors

        return ret;
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
        if (arg<-1.0) arg = -1.0;
        if (arg>1.0) arg = 1.0;
        double theta = acos(arg);
        
        double ret = G_ip.get_interp_value_tricubic(r,rb,theta);

        return ret;
    }
}