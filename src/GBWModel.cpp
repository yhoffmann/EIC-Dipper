#include "../include/utilities.h"
#include "../include/constants.h"
#include <math.h>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>


namespace GBWModel {
    double T_times_sigma0 (double b1, double b2) {
        return exp( -( sqr(b1)+sqr(b2) )/(2.0*BG) );
    }

    double Q_s_sqr (double b1, double b2) {
        return sqr(Qs0) * T_times_sigma0(b1,b2);
    }

    double G_old (double x1, double x2, double y1, double y2) {
        return -0.25 * Q_s_sqr( (x1+y1)/2.0, (x2+y2)/2.0 ) * ( sqr(x1-y1) + sqr(x2-y2) );
    }

    double Kmod (double x1, double x2, double y1, double y2) {
        return 4.0* ( 1.0-m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))*bessel_K_safe(1,m*std::sqrt(sqr(x1-y1)+sqr(x2-y2))) ) / sqr(m) / (sqr(x1-y1)+sqr(x2-y2));
    }

    double G (double x1, double x2, double y1, double y2) {
        return -0.25 * exp( -(sqr(x1)+sqr(x2))/(2*BG) ) * exp( -(sqr(y1)+sqr(y2))/(2*BG) ) * ( sqr(x1-y1) + sqr(x2-y2) ) * Kmod(x1,x2,y1,y2);
    }
}