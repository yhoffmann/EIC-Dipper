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
    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double z, double Delta, bool t_not_l) {        
        return 1.0/(4.0*PI) * NRPhoton::wave_function(r1,r2,Q,z,t_not_l) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta )/(2.0*PI) * SaturationModel::dsigma_d2b(b1+r1/2.0,b2+r2/2.0,b1-r1/2.0,b2-r2/2.0);
    }


    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
        IntegrationConfig* integration_config = (IntegrationConfig*) userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*) integration_config->integrand_params;

        double r_range_factor = R_RANGE_FACTOR;
        double r_significant_range = 1.0/NRPhoton::epsilon(A_integrand_params->Q,A_integrand_params->z);

        double bmin = integration_config->min;
        double bmax = integration_config->max;
        double phibmin = 0.0;
        double phibmax = 2*PI;

        double rmin = 0.0;
        double rmax = r_range_factor*r_significant_range;
        double phirmin = 0.0;
        double phirmax = 2.0*PI;

        double b = bmin + (bmax-bmin)*xx[0];
        double phib = phibmin + (phibmax-phibmin)*xx[1];
        double r = rmin + (rmax-rmin)*xx[2];
        double phir = phirmin + (phirmax-phirmin)*xx[3];

        double b1 = b*cos(phib);
        double b2 = b*sin(phib);
        double r1 = r*cos(phir);
        double r2 = r*sin(phir);

        double jacobian = r*b*(bmax-bmin)*(phibmax-phibmin)*(rmax-rmin)*(phirmax-phirmin);

        ff[0] = jacobian * Coherent::A_integrand_function(b1,b2,r1,r2,A_integrand_params->Q,A_integrand_params->z,A_integrand_params->Delta,A_integrand_params->t_not_l);

        return 0;
    }


    std::vector<double> dsigma_dt (CubaConfig c_config, IntegrationConfig integration_config) {
        std::vector<double> ret(2);
        c_config.num_of_dims = 4;

        AIntegrandParams A_integrand_params = *(AIntegrandParams*) integration_config.integrand_params;
        integration_config.integrand_params = &A_integrand_params;

        A_integrand_params.t_not_l = true;

        ret[0] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[0] = sqr(ret[0])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        A_integrand_params.t_not_l = false;

        ret[1] = IntegrationRoutines::cuba_integrate_one_bessel(Coherent::integrand,c_config,integration_config);
        ret[1] = sqr(ret[1])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}