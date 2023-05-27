#include "../include/incoherent.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "../include/constants.h"
#include "../include/models_and_functions.h"
#include "../include/coherent.h"
#include "../include/utilities.h"


namespace Incoherent {
    double DD (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        double a = CorrelationMatrixElements::a(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
        double b = CorrelationMatrixElements::b(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
        double c = CorrelationMatrixElements::c(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
        double d = CorrelationMatrixElements::d(b1,b2,r1,r2,bb1,bb2,rb1,rb2);

        double Sqrt = std::sqrt(4.0*b*c+std::pow(a-d,2.0));
        double aplusd = a+d;
        double factor = (a-d+2.0/Nc*b)/Sqrt;

        //double ret = 0.5* ( (1.0+factor)*exp((aplusd+Sqrt)/2.0) + (1.0-factor)*exp((aplusd-Sqrt)/2.0) );

        double ret = exp(aplusd/2.0) * ( cosh(0.5*Sqrt) + factor*sinh(0.5*Sqrt) );

        if (gsl_isnan(ret)) {
            std::cerr << ret << " " << b1 << " " << b2 << " " << bb1 << " " << bb2 << " " << r1 << " " << r2 << " " << rb1 << " " << rb2 << std::endl;
            std::cout << 4.0*b*c+std::pow(a-d,2.0) << " " << aplusd << " " << factor;
            std::cerr << "DipoleDipole is nan. Aborting" << std::endl;
            exit(0);
        }

        return ret;
    }


    double dsigma_d2b_sqr (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        using namespace CorrelationMatrixElements;
        
        return 4.0 * (1.0 -exp(G_xy(b1,b2,r1,r2,bb1,bb2,rb1,rb2)) -exp(G_xbyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2)) + DD(b1,b2,r1,r2,bb1,bb2,rb1,rb2));
    }

    double A_integrand_function_simple (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, bool t_not_l) {
        return 1.0/(16.0*PI*PI) * gsl_sf_bessel_J0(std::sqrt( sqr(b1-bb1)+sqr(b2-bb2) )*Delta)/(4.0*PI*PI) * NRPhoton::wave_function(r1,r2,Q,z,t_not_l) * NRPhoton::wave_function(rb1,rb2,Q,z,t_not_l) * Coherent::dsigma_d2b(b1,b2,r1,r2) * Coherent::dsigma_d2b(bb1,bb2,rb1,rb2);
    }

    double A_integrand_function (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z, double Delta, bool t_not_l) {
        return 1.0/(16.0*PI*PI) * gsl_sf_bessel_J0(std::sqrt( sqr(b1-bb1)+sqr(b2-bb2) )*Delta)/(4.0*PI*PI) * NRPhoton::wave_function(r1,r2,Q,z,t_not_l) * NRPhoton::wave_function(rb1,rb2,Q,z,t_not_l) * dsigma_d2b_sqr(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
    }

    int integrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata) {
        IntegrandParams* i_params = (IntegrandParams*) userdata;

        double r_range_factor = get_r_range_factor();
        double r_significant_range = 1.0/NRPhoton::epsilon(i_params->Q,i_params->z);

        double dbmin = i_params->min;
        double dbmax = i_params->max;
        double phidbmin = 0.0;
        double phidbmax = 2*PI;

        double rmin = 0.0;
        double rmax = r_range_factor*r_significant_range;
        double phirmin = 0.0;
        double phirmax = 2.0*PI;

        double Bmin = 0.0;
        double Bmax = BMAX;
        double phiBmin = 0.0;
        double phiBmax = 2*PI;

        double rbmin = 0.0;
        double rbmax = rmax;
        double phirbmin = 0.0;
        double phirbmax = 2.0*PI;

        double db = dbmin + (dbmax-dbmin)*xx[0];
        double phib = phidbmin + (phidbmax-phidbmin)*xx[1];
        double r = rmin + (rmax-rmin)*xx[2];
        double phir = phirmin + (phirmax-phirmin)*xx[3];

        double B = Bmin + (Bmax-Bmin)*xx[4];
        double phiB = phiBmin + (phiBmax-phiBmin)*xx[5];
        double rb = rbmin + (rbmax-rbmin)*xx[6];
        double phirb = phirbmin + (phirbmax-phirbmin)*xx[7];

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

        double jacobian = r*db*B*rb*(dbmax-dbmin)*(phidbmax-phidbmin)*(rmax-rmin)*(phirmax-phirmin)*(Bmax-Bmin)*(phiBmax-phiBmin)*(rbmax-rbmin)*(phirbmax-phirbmin);

        ff[0] = jacobian * Incoherent::A_integrand_function(b1,b2,r1,r2,bb1,bb2,rb1,rb2,i_params->Q,i_params->z,i_params->Delta,i_params->t_not_l);

        return 0;
    }

    std::vector<double> calculate_dsigma_dt (CubaConfig c_config, IntegrandParams i_params) {
        std::vector<double> ret(2);
        c_config.num_of_dims = 8;

        i_params.t_not_l = true;

        ret[0] = Routines::cuba_integrate_one_bessel(Incoherent::integrand,c_config,i_params);
        ret[0] *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        i_params.t_not_l = false;

        ret[1] = Routines::cuba_integrate_one_bessel(Incoherent::integrand,c_config,i_params);
        ret[1] *= (GeVm1Tofm*GeVm1Tofm*fm2TonB)/(16.0*PI);

        return ret;
    }
}