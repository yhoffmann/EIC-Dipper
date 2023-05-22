#pragma once


#include<cuba.h>
#include"integration_routines.cpp"
#include"constants.cpp"
#include"models_and_functions.cpp"
#include"utilities.cpp"

/*
static int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
    IntegrandConfig* c = (IntegrandConfig*) userdata;

    double r_range_factor = get_r_range_factor();
    double r_significant_range = 1.0/Photon::epsilon(c->params.Q,c->params.z);

    double bmin = c->lower_bound;
    double bmax = c->upper_bound;

    double phibmin = 0.0;
    double phibmax = 2*PI;

    double rmin = 0.0;
    double rmax = r_range_factor*r_significant_range;

    double phirmin = 0.0;
    double phirmax = 2*PI;

    double b = bmin+(bmax-bmin)*xx[0];
    double phib = phibmin+(phibmax-phibmin)*xx[1];

    double r = rmin+(rmax-rmin)*xx[1];
    double phir = phirmin+(phirmax-phirmin)*xx[3];

    double b1 = b*cos(phib);
    double b2 = b*sin(phib);

    double r1 = r*cos(phir);
    double r2 = r*sin(phir);

    double Jacobian = r*b*(bmax-bmin)*(phibmax-phibmin)*(rmax-rmin)*(phirmax-phirmin);

    if (c->func_conf.co_or_inco == CoOrInco::Co) {
        if (c->func_conf.t_or_l == TOrL::T) {
            ff[0] = Jacobian * Coherent::Trans::AIntegrandFunction(b1,b2,r1,r2,c->params.Q,c->params.z);
        }
        if (c->func_conf.t_or_l == TOrL::L) {
            ff[0] = Jacobian * Coherent::Longi::AIntegrandFunction(b1,b2,r1,r2,c->params.Q,c->params.z);
        }
    }

    if (c->func_conf.co_or_inco == CoOrInco::Inco) {
        double Bmin = 0.0;
        double Bmax = BMAX;

        double phiBmin = 0.0;
        double phiBmax = 2*PI;

        double Rmin = 0.0;
        double Rmax = r_range_factor*r_significant_range;

        double phiRmin = 0.0;
        double phiRmax = 2*PI;

        double B = Bmin+(Bmax-Bmin)*xx[4];
        double phiB = phiBmin+(phiBmax-phiBmin)*xx[5];

        double R = Rmin+(Rmax-Rmin)*xx[6];
        double phiR = phiRmin+(phiRmax-phiRmin)*xx[7];

        double B1 = B*cos(phiB);
        double B2 = B*sin(phiB);

        double R1 = R*cos(phiR);
        double R2 = R*sin(phiR);

        Jacobian *= R*B*(Bmax-Bmin)*(phiBmax-phiBmin)*(Rmax-Rmin)*(phiRmax-phiRmin);

        b1 = B1+b1/2.0;
        b2 = B2+b2/2.0;

        B1 = B1-b1/2.0;
        B2 = B2-b2/2.0;

        if (c->func_conf.t_or_l == TOrL::T) {
            ff[0] = Jacobian * Incoherent::Trans::AIntegrandFunction(b1,b2,B1,B2,r1,r2,R1,R2,c->params.Q,c->params.z);
        }
        if (c->func_conf.t_or_l == TOrL::L) {
            ff[0] = Jacobian * Incoherent::Longi::AIntegrandFunction(b1,b2,B1,B2,r1,r2,R1,R2,c->params.Q,c->params.z);
        }
    }


    
    return 0;
}
*/


int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
    IntegrandParams* i_params = (IntegrandParams*) userdata;
    
    double args[8]{0.0};

    double r_range_factor = get_r_range_factor();
    double r_significant_range = 1.0/Photon::epsilon(i_params->Q,i_params->z);

    double lower_bound[8];
    lower_bound[0] = i_params->min; //bmin
    lower_bound[1] = 0.0; //phibmin
    lower_bound[2] = 0.0; //rmin
    lower_bound[3] = 0.0; //phirmin
    lower_bound[4] = 0.0; //Bmin
    lower_bound[5] = 0.0; //phiBmin
    lower_bound[6] = 0.0; //Rmin
    lower_bound[7] = 0.0; //phiRmin

    double upper_bound[8];
    upper_bound[0] = i_params->max; //bmax
    upper_bound[1] = 2.0*PI; //phibmax
    upper_bound[2] = r_range_factor*r_significant_range; //rmax
    upper_bound[3] = 2.0*PI; //phirmax
    upper_bound[4] = BMAX; //Bmax
    upper_bound[5] = 2.0*PI; //phiBmax
    upper_bound[6] = r_range_factor*r_significant_range; //Rmax
    upper_bound[7] = 2.0*PI; //phiRmax

    double jacobian = 1.0;

    for (int n=0; n<*ndim; n++) {
        args[n] = lower_bound[n]+(upper_bound[n]-lower_bound[n])*xx[n];
        
        jacobian *= (upper_bound[n]-lower_bound[n]);
        
        if (n%2==0) {
            jacobian *= args[n];
        }
    }

    ff[0] = jacobian * i_params->func(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],i_params->Q,i_params->z);

    return 0;
}