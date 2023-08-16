#include "../include/Coherent.h"
#include <vector>
#include <math.h>
#include <cuba.h>
#include "../include/constants.h"
#include "../include/IntegrationRoutines.h"
#include "../include/GBWModel.h"
#include "../include/NRPhoton.h"
#include "../include/SaturationModel.h"
#include <gsl/gsl_math.h>
#include <iostream>


namespace Coherent {
    double A_integrand_function (double x1, double x2, double y1, double y2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal) {        
        return 1.0/(4.0*PI) * NRPhoton::wave_function(x1-y1,x2-y2,Q,z,transverse_or_longitudinal) * gsl_sf_bessel_J0( 0.25*std::sqrt(sqr(x1+y1)+sqr(x2+y2))*Delta )/(2.0*PI) * SaturationModel::dsigma_d2b(x1,x2,y1,y2);
    }


    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
        IntegrationConfig* integration_config = (IntegrationConfig*) userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*) integration_config->integrand_params;

        double r_range_factor = R_RANGE_FACTOR;
        double r_significant_range = 1.0/NRPhoton::epsilon(A_integrand_params->Q,A_integrand_params->z);

        double xmin = integration_config->min;
        double xmax = integration_config->max;
        //double phixmin = 0.0;
        double phixmax = 2*PI;

        //double ymin = 0.0;
        double ymax = r_range_factor*r_significant_range;
        //double phiymin = 0.0;
        double phiymax = 2.0*PI;

        double x = xmax*xx[0];
        double phix = phixmax*xx[1];
        double y = ymax*xx[2];
        double phiy = phiymax*xx[3];

        double x1 = x*cos(phix);
        double x2 = x*sin(phix);
        double y1 = y*cos(phiy);
        double y2 = y*sin(phiy);

        double jacobian = x*y*xmax*phixmax*ymax*phiymax;

        ff[0] = jacobian * Coherent::A_integrand_function(x1,y2,y1,y2,A_integrand_params->Q,A_integrand_params->z,A_integrand_params->Delta,A_integrand_params->transverse_or_longitudinal);

        return 0;
    }


    std::vector<double> dsigma_dt (CubaConfig* c_config, IntegrationConfig* integration_config) {
        std::vector<double> ret(2);
        c_config->num_of_dims = 4;

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = T;

        ret[0] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[0] = sqr(ret[0])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = L;

        ret[1] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[1] = sqr(ret[1])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}