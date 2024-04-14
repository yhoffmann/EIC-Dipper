#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include "../include/SaturationModel.hpp"
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"


namespace SaturationModel
{
    inline double D(double G)
    {
#ifndef _DILUTE
        return exp(G);
#else
        return 1.0+G;
#endif
    }


    double dsigma_d2b(double x1, double x2, double y1, double y2)
    {
        return 2.0*( 1.0-D(NH*GBWModel::G(x1, x2, y1, y2)) );
    }


    double dsigma_d2b_sqr (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xy = NH*GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = NH*GBWModel::G(xb1, xb2, yb1, yb2);
        double G_xxb = NH*GBWModel::G(x1, x2, xb1, xb2);
        double G_xyb = NH*GBWModel::G(x1, x2, yb1, yb2);
        double G_xby = NH*GBWModel::G(xb1, xb2, y1, y2);
        double G_yyb = NH*GBWModel::G(y1, y2, yb1, yb2);

#ifndef _DILUTE
        double T_xy_ybxb = G_xxb + G_yyb - G_xyb - G_xby;
        double T_xxb_yby = G_xy + G_xbyb - G_xyb - G_xby;
        
        double a = G_xy + G_xbyb - T_xy_ybxb*Ncsqrm1_inverse;
        double b = T_xy_ybxb*twoCF_inverse;
        double c = T_xxb_yby*twoCF_inverse;
        double d = G_xxb + G_yyb - T_xxb_yby*Ncsqrm1_inverse;

        double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
        double factor = ( a-d+2.0*b ) / Sqrt;

        double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );
        
        return 4.0 * (1.0-D(G_xy)-D(G_xbyb)+DD);
#else
        return 4.0*(G_xy*G_xbyb+sqr(G_xxb+G_yyb-G_xyb-G_xby)/2.0);
#endif
    }


    double dsigma_d2b_sqr_reduced (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
#ifndef _DILUTE
        double G_xy = NH*GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = NH*GBWModel::G(xb1, xb2, yb1, yb2);
#endif
        double G_xxb = NH*GBWModel::G(x1, x2, xb1, xb2);
        double G_xyb = NH*GBWModel::G(x1, x2, yb1, yb2);
        double G_xby = NH*GBWModel::G(xb1, xb2, y1, y2);
        double G_yyb = NH*GBWModel::G(y1, y2, yb1, yb2);

#ifndef _DILUTE
        double T_xy_ybxb = G_xxb + G_yyb - G_xyb - G_xby;
        double T_xxb_yby = G_xy + G_xbyb - G_xyb - G_xby;
        
        double a = G_xy + G_xbyb - T_xy_ybxb*Ncsqrm1_inverse;
        double b = T_xy_ybxb*twoCF_inverse;
        double c = T_xxb_yby*twoCF_inverse;
        double d = G_xxb + G_yyb - T_xxb_yby*Ncsqrm1_inverse;

        double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
if (Sqrt == 0.0) std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << " " << xb1 << " " << xb2 << " " << yb1 << " " << yb2 << "\n" << G_xy << " " << G_xbyb << " " << G_xxb << " " << G_xyb << " " << G_xby << " " << G_yyb << "\n" << a << " " << b << " " << c << " " << d << std::endl;
        double factor = ( a-d+2.0*b ) / Sqrt;

        double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );
        
        return 4.0 * ( DD-D(G_xy)*D(G_xbyb) );
#else
        return 2.0*sqr(G_xxb+G_yyb-G_xyb-G_xby);
#endif
    }


    namespace Sampled
    {
        double dsigma_d2b (double x1, double x2, double y1, double y2, const HotspotNucleus* nucleus)
        {
            double G_sum = 0.0;
            
            for (uint n=0, a=nucleus->get_atomic_num(); n<a; ++n)
                for (uint i=0, nh=nucleus->get_num_hotspots_per_nucleon(); i<nh; ++i)
                {
                    HotspotPos b0 = *nucleus->get_hotspot_pos(n, i);

                    G_sum += GBWModel::G(x1-b0.x, x2-b0.y, y1-b0.x, y2-b0.y);
                }

            return 2.0*( 1.0-D(G_sum) );
        }


        double dsigma_d2b_sqr (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus)
        {
            double G_xy_sum = 0.0;
            double G_xbyb_sum = 0.0;
            double G_xxb_sum = 0.0;
            double G_xyb_sum = 0.0;
            double G_xby_sum = 0.0;
            double G_yyb_sum = 0.0;

            for (uint n=0, a=nucleus->get_atomic_num(); n<a; ++n)
                for (uint i=0, nh=nucleus->get_num_hotspots_per_nucleon(); i<nh; ++i)
                {
                    HotspotPos b0 = *nucleus->get_hotspot_pos(n, i);

                    double x1_mod = x1-b0.x;
                    double x2_mod = x2-b0.y;
                    double y1_mod = y1-b0.x;
                    double y2_mod = y2-b0.y;
                    double xb1_mod = xb1-b0.x;
                    double xb2_mod = xb2-b0.y;
                    double yb1_mod = yb1-b0.x;
                    double yb2_mod = yb2-b0.y;

                    G_xy_sum += GBWModel::G(x1_mod, x2_mod, y1_mod, y2_mod);
                    G_xbyb_sum += GBWModel::G(xb1_mod, xb2_mod, yb1_mod, yb2_mod);
                    G_xxb_sum += GBWModel::G(x1_mod, x2_mod, xb1_mod, xb2_mod);
                    G_xyb_sum += GBWModel::G(x1_mod, x2_mod, yb1_mod, yb2_mod);
                    G_xby_sum += GBWModel::G(xb1_mod, xb2_mod, y1_mod, y2_mod);
                    G_yyb_sum += GBWModel::G(y1_mod, y2_mod, yb1_mod, yb2_mod);
                }
    #ifndef _DILUTE
            double T_xy_ybxb = G_xxb_sum + G_yyb_sum - G_xyb_sum - G_xby_sum;
            double T_xxb_yby = G_xy_sum + G_xbyb_sum - G_xyb_sum - G_xby_sum;
            
            double a = G_xy_sum + G_xbyb_sum - T_xy_ybxb*Ncsqrm1_inverse;
            double b = T_xy_ybxb*twoCF_inverse;
            double c = T_xxb_yby*twoCF_inverse;
            double d = G_xxb_sum + G_yyb_sum - T_xxb_yby*Ncsqrm1_inverse;

            double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
            double factor = ( a-d+2.0*b ) / Sqrt;

            double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );
            
            return 4.0 * ( 1.0-D(G_xy_sum)-D(G_xbyb_sum)+DD );
    #else
            return 4.0*(G_xy_sum*G_xbyb_sum+sqr(G_xxb_sum+G_yyb_sum-G_xyb_sum-G_xby_sum)/2.0);
    #endif
        }


        double dsigma_d2b_sqr_reduced (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus)
        {
    #ifndef _DILUTE
            double G_xy_sum = 0.0;
            double G_xbyb_sum = 0.0;
    #endif
            double G_xxb_sum = 0.0;
            double G_xyb_sum = 0.0;
            double G_xby_sum = 0.0;
            double G_yyb_sum = 0.0;

            for (uint n=0, a=nucleus->get_atomic_num(); n<a; ++n)
                for (uint i=0, nh=nucleus->get_num_hotspots_per_nucleon(); i<nh; ++i)
                {
                    HotspotPos b0 = *nucleus->get_hotspot_pos(n, i);

                    double x1_mod = x1-b0.x;
                    double x2_mod = x2-b0.y;
                    double y1_mod = y1-b0.x;
                    double y2_mod = y2-b0.y;
                    double xb1_mod = xb1-b0.x;
                    double xb2_mod = xb2-b0.y;
                    double yb1_mod = yb1-b0.x;
                    double yb2_mod = yb2-b0.y;
            #ifndef _DILUTE
                    G_xy_sum += GBWModel::G(x1_mod, x2_mod, y1_mod, y2_mod);
                    G_xbyb_sum += GBWModel::G(xb1_mod, xb2_mod, yb1_mod, yb2_mod);
            #endif
                    G_xxb_sum += GBWModel::G(x1_mod, x2_mod, xb1_mod, xb2_mod);
                    G_xyb_sum += GBWModel::G(x1_mod, x2_mod, yb1_mod, yb2_mod);
                    G_xby_sum += GBWModel::G(xb1_mod, xb2_mod, y1_mod, y2_mod);
                    G_yyb_sum += GBWModel::G(y1_mod, y2_mod, yb1_mod, yb2_mod);
                }
    #ifndef _DILUTE
            double T_xy_ybxb = G_xxb_sum + G_yyb_sum - G_xyb_sum - G_xby_sum;
            double T_xxb_yby = G_xy_sum + G_xbyb_sum - G_xyb_sum - G_xby_sum;

            double a = G_xy_sum + G_xbyb_sum - T_xy_ybxb*Ncsqrm1_inverse;
            double b = T_xy_ybxb*twoCF_inverse;
            double c = T_xxb_yby*twoCF_inverse;
            double d = G_xxb_sum + G_yyb_sum - T_xxb_yby*Ncsqrm1_inverse;

            double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
            double factor = ( a-d+2.0*b ) / Sqrt;

            double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );
            
            return 4.0 * ( DD-D(G_xy_sum)*D(G_xbyb_sum) );
    #else
            return 2.0*sqr(G_xxb_sum+G_yyb_sum-G_xyb_sum-G_xby_sum);
    #endif
        }
    }
}