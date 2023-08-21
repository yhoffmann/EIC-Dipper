#include "../include/Incoherent.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "../include/NRPhoton.h"
#include "../include/constants.h"
#include "../include/SaturationModel.h"
#include "../include/utilities.h"


namespace Incoherent {
    double A_integrand_function_simple (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal) {
        return 1.0/(16.0*PI*PI) * gsl_sf_bessel_J0( std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta )/(4.0*PI*PI) * NRPhoton::wave_function(r1,r2,Q,z,transverse_or_longitudinal) * NRPhoton::wave_function(rb1,rb2,Q,z,transverse_or_longitudinal) * SaturationModel::dsigma_d2b(b1,b2,r1,r2) * SaturationModel::dsigma_d2b(bb1,bb2,rb1,rb2);
    }

    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, TransverseOrLongitudinal transverse_or_longitudinal) {
        double x1 = b1+r1/2.0;
        double x2 = b2+r2/2.0;
        double y1 = b1-r1/2.0;
        double y2 = b2-r2/2.0;
        double xb1 = bb1+rb1/2.0;
        double xb2 = bb2+rb2/2.0;
        double yb1 = bb1-rb1/2.0;
        double yb2 = bb2-rb2/2.0;

        return 1.0/(64.0*PI*PI*PI*PI) * gsl_sf_bessel_J0( std::sqrt(sqr(b1-bb1)+sqr(b2-bb2))*Delta ) * NRPhoton::wave_function(r1,r2,Q,z,transverse_or_longitudinal) * NRPhoton::wave_function(rb1,rb2,Q,z,transverse_or_longitudinal) * SaturationModel::dsigma_d2b_sqr(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
    }

    int integrand (const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata) {
        IntegrationConfig* integration_config = (IntegrationConfig*)userdata;

        AIntegrandParams* A_integrand_params = (AIntegrandParams*)(integration_config->integrand_params);

        double dbmin = integration_config->min;
        double dbmax = integration_config->max;
        //double phidbmin = 0.0;
        //double phidbmax = 2*PI;

        //double rmin = 0.0;
        //double rmax = R_MAX;
        //double phirmin = 0.0;
        //double phirmax = 2.0*PI;

        //double Bmin = 0.0;
        //double Bmax = B_MAX;
        //double phiBmin = 0.0;
        //double phiBmax = 2*PI;

        //double rbmin = 0.0;
        //double rbmax = rmax;
        //double phirbmin = 0.0;
        //double phirbmax = 2.0*PI;

        double db = dbmin + (dbmax-dbmin)*xx[0];
        double phib = 2.0*PI*xx[1];
        double r = R_MAX*xx[2];
        double phir = 2.0*PI*xx[3];

        double B = B_MAX*xx[4];
        double phiB = 2.0*PI*xx[5];
        double rb = R_MAX*xx[6];
        double phirb = 2.0*PI*xx[7];

        double db1 = db*cos(phib);
        double db2 = db*sin(phib);
        double r1 = r*cos(phir);
        double r2 = r*sin(phir);

        double B1 = B*cos(phiB);
        double B2 = B*sin(phiB);
        double rb1 = rb*cos(phirb);
        double rb2 = rb*sin(phirb);

        double b1 = B1+db1/2.0;
        double b2 = B2+db2/2.0;
        double bb1 = B1-db1/2.0;
        double bb2 = B2-db2/2.0;

        double jacobian = r*db*B*rb*(dbmax-dbmin)*R_MAX*B_MAX*R_MAX*16.0*PI*PI*PI*PI;

        ff[0] = jacobian * Incoherent::A_integrand_function(b1,b2,r1,r2,bb1,bb2,rb1,rb2,A_integrand_params->Q,A_integrand_params->z,A_integrand_params->Delta,A_integrand_params->transverse_or_longitudinal);

        return 0;
    }

    std::vector<double> dsigma_dt (CubaConfig* c_config, IntegrationConfig* integration_config) {
        std::vector<double> ret(2);
        c_config->num_of_dims = 8;

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = T;

        ret[0] = IntegrationRoutines::cuba_integrate_one_bessel(Incoherent::integrand,c_config,integration_config);
        ret[0] *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        ((AIntegrandParams*)(integration_config->integrand_params))->transverse_or_longitudinal = L;

        ret[1] = IntegrationRoutines::cuba_integrate_one_bessel(Incoherent::integrand,c_config,integration_config);
        ret[1] *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}