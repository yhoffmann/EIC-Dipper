#include "../include/Coherent.h"
#include <math.h>
#include "../include/constants.h"
#include "../include/IntegrationRoutines.h"
#include "../include/GBWModel.h"
#include "../include/NRPhoton.h"
#include "../include/SaturationModel.h"
#include <gsl/gsl_math.h>
#include <iostream>


namespace Coherent {
    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal) {
        double x1 = b1+r1/2.0;
        double x2 = b2+r2/2.0;
        double y1 = b1+r1/2.0;
        double y2 = b2+r2/2.0;

        return 1.0/(8.0*PI*PI) * NRPhoton::wave_function(r1,r2,Q,z,transverse_or_longitudinal) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta ) * SaturationModel::dsigma_d2b(x1,x2,y1,y2);
    }


    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata) {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(integration_config->integrand_params);

        double bmin = integration_config->min;
        double bmax = integration_config->max;
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

        ff[0] = jacobian * Coherent::A_integrand_function(b1,b2,r1,r2,A_integrand_params->Q,A_integrand_params->z,A_integrand_params->Delta,A_integrand_params->transverse_or_longitudinal);

        return 0;
    }


    std::vector<double> dsigma_dt (CubaConfig* c_config, IntegrationConfig* integration_config) {
        std::vector<double> ret(2);
        c_config->num_of_dims = 4;

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = T;

        ret[0] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[0] = sqr(ret[0]*GeVm1Tofm)*fm2TonB/(16.0*PI);

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = L;

        ret[1] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[1] = sqr(ret[1]*GeVm1Tofm)*fm2TonB/(16.0*PI);

        return ret;
    }
}