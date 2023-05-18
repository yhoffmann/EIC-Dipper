#ifndef INTEGRATION_CPP
#define INTEGRATION_CPP


#define SIX_D double,double,double,double,double,double
#define TEN_D double,double,double,double,double,double,double,double,double,double

#include<vector>
#include<cuba.h>
#include<iostream>
#include<gsl/gsl_math.h>
#include"utilities.cpp"
#include"coherent.cpp"
#include"incoherent.cpp"


enum Model {
    GBW
};


enum CoOrInco {
    Co,
    Inco
};


enum TOrL {
    T,
    L
};


enum Integrator {
    C,
    S
};


struct CubaConfig {
    Integrator integrator;
    double epsrel = 1.0e-6;
    double epsabs = 1.0e-8;
    luint mineval = 2e5;
    luint maxeval = 1e6;
    int flags1 = 0;
    int flasg2 = 4;
    int seed = 0;
    double bessel_tolerance = 1e-5;
};


struct FunctionParameters {
    double Q;
    double z;
    double Dx;
    double Dy;
};


struct FunctionConfig {
    Model model = Model::GBW;
    CoOrInco co_or_inco;
    TOrL t_or_l;
};


struct IntegrandConfig {
    FunctionParameters params;
    
    double lower_bound;
    double upper_bound;
};


struct IntegrationConfig {
    IntegrandConfig int_conf;
    FunctionConfig func_conf;
    CubaConfig cuba_config;
};


int get_num_of_dim (FunctionConfig& func_conf) {
    switch (func_conf.co_or_inco) {
        case CoOrInco::Co:
            return 4;
        case CoOrInco::Inco:
            return 8;
    }
}


struct IntegrandWrapper {
    double (*func)(TEN_D);

    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
        IntegrandConfig* c = (IntegrandConfig*) userdata;
        
        double args[8]{0.0};

        double r_range_factor = get_r_range_factor();
        double r_significant_range = 1.0/Photon::epsilon(c->params.Q,c->params.z);

        double lower_bound[8];
        lower_bound[0] = c->lower_bound; //bmin
        lower_bound[1] = 0.0; //phibmin
        lower_bound[2] = 0.0; //rmin
        lower_bound[3] = 0.0; //phirmin
        lower_bound[4] = 0.0; //Bmin
        lower_bound[5] = 0.0; //phiBmin
        lower_bound[6] = 0.0; //Rmin
        lower_bound[7] = 0.0; //phiRmin

        double upper_bound[8];
        upper_bound[0] = c->upper_bound; //bmax
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

        ff[0] = jacobian * func(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],c->params.Q,c->params.z);

        return 0;
    }
};


void get_integrand_function (double (*func)(TEN_D), FunctionConfig& func_conf) {
    if (func_conf.co_or_inco == Co) {
        if (func_conf.t_or_l == T) {
            func = &Coherent::Trans::AIntegrandWrapper;
        } else {
            func = &Coherent::Longi::AIntegrandWrapper;
        }
    } else {
        if (func_conf.t_or_l == T) {
            func = &Incoherent::Trans::AIntegrandWrapper;
        } else {
            func = &Incoherent::Longi::AIntegrandWrapper;
        }
    }
}


std::vector<double> cuba_integrate (IntegrationConfig conf) {
    std::vector<double> ret;

    int num_of_regions(0), num_of_evals(0), error_status(0);

    int num_of_integrals = 1;



    return ret;
}


std::vector<double> cuba_bessel_integrate (IntegrationConfig conf) {
    std::vector<double> ret;

    IntegrandWrapper integrand_wrapper;
    get_integrand_function(integrand_wrapper.func, conf.func_conf);

    

    return ret;
}


#endif