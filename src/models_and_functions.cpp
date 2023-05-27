#include "../include/models_and_functions.h"
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

    double G (double x1, double x2, double y1, double y2) {
        return -0.25 * Q_s_sqr( (x1+y1)/2.0, (x2+y2)/2.0 ) * ( sqr(x1-y1) + sqr(x2-y2) );
    }
}


namespace NRPhoton {
    double epsilon (double Q, double z) {
        return std::sqrt( sqr(Q)*z*(1.0-z) + sqr(m_Q_c) );
    }

    double wave_function (double r1, double r2, double Q, double z, bool t_not_l) {
        if (sqr(r1)+sqr(r2) < 2.0e-40) { r1 = 1.0e-20; r2 = 1.0e-20;}
        
        if (t_not_l) {
            return -A_Q * sqrt_2m_c_Nc * e * e_Q * gsl_sf_bessel_K0( epsilon(Q,z) * std::sqrt( sqr(r1)+sqr(r2) ) );
        } else {
            return -2.0 / m_Q_c * Q * z * (1.0-z) * A_Q * sqrt_2m_c_Nc * e * e_Q * gsl_sf_bessel_K0( epsilon(Q,z) * std::sqrt( sqr(r1)+sqr(r2) ) );
        }
    }
}


namespace CorrelationMatrixElements
{
    using namespace GBWModel;

    double G_xy (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(x(b1,r1),x(b2,r2),y(b1,r1),y(b2,r2));
    }

    double G_xyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(x(b1,r1),x(b2,r2),y(bb1,rb1),y(bb2,rb2));
    }

    double G_xby (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(x(bb1,rb1),x(bb2,rb2),y(b1,r1),y(b2,r2));
    }

    double G_xbyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(x(bb1,rb1),x(bb2,rb2),y(bb1,rb1),y(bb2,rb2));
    }

    double G_xxb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(x(b1,r1),x(b2,r2),x(bb1,rb1),x(bb2,rb2));
    }

    double G_yyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G(y(b1,r1),y(b2,r2),y(bb1,rb1),y(bb2,rb2));
    }

    double T1 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G_xyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xby(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xxb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_yyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
    }

    double T2 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G_xy(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xbyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xxb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_yyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
    }


    double a (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G_xy(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xbyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - (T1(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(Nc*Nc-1.0);
    }

    double b (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return (T1(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(2.0*CF);
    }

    double c (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {       
        return (T2(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(2.0*CF);
    }

    double d (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
        return G_xyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xby(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - (T2(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(Nc*Nc-1.0);
    }
}