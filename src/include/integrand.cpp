#ifndef INTEGRAND_CPP
#define INTEGRAND_CPP


#include<cuba.h>
#include"integration_routines.cpp"
#include"constants.cpp"
#include"models_and_functions.cpp"


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


#endif