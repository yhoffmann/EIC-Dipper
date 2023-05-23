#include "../include/coherent.h"
#include <vector>
#include <math.h>
#include <cuba.h>
#include "../include/constants.h"
#include "../include/integration_routines.h"
#include "../include/models_and_functions.h"


namespace Coherent {
    double dsigma_d2b (double b1, double b2, double r1, double r2) {
        return 2.0*(1.0 - exp(GBWModel::G( x(b1,r1),x(b2,r2),y(b1,r1),y(b2,r2) )));
    }

    double A_integrand_function (double b1, double b2, double r1, double r2, double Q, double z, double Delta, bool t_not_l) {
        return 1.0/(4.0*PI) * NRPhoton::wave_function(r1,r2,Q,z,t_not_l) * gsl_sf_bessel_J0( std::sqrt(sqr(b1)+sqr(b2))*Delta )/(2.0*PI) * dsigma_d2b(b1,b2,r1,r2);
    }


    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
        IntegrandParams* i_params = (IntegrandParams*) userdata;

        double r_range_factor = get_r_range_factor();
        double r_significant_range = 1.0/NRPhoton::epsilon(i_params->Q,i_params->z);

        double bmin = i_params->min;
        double bmax = i_params->max;
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

        ff[0] = jacobian * Coherent::A_integrand_function(b1,b2,r1,r2,i_params->Q,i_params->z,i_params->Delta,i_params->t_not_l);

        return 0;
    }


    std::vector<double> calculate_dsigma_dt (CubaConfig c_config, IntegrandParams i_params) {
        std::vector<double> ret(2);
        c_config.num_of_dims = 4;

        i_params.t_not_l = true;

        ret[0] = Routines::cuba_integrate_one_bessel(Coherent::integrand,c_config,i_params);
        ret[0] = sqr(ret[0])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        i_params.t_not_l = false;

        ret[1] = Routines::cuba_integrate_one_bessel(Coherent::integrand,c_config,i_params);
        ret[1] = sqr(ret[1])*(GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}