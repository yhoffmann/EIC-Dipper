#include "../include/SaturationModel.hpp"
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>


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
        double G_yxb = NH*GBWModel::G(xb1, xb2, y1, y2);
        double G_yyb = NH*GBWModel::G(y1, y2, yb1, yb2);
#ifndef _DILUTE
        double T_xy_xbyb = G_xyb + G_yxb - G_xxb - G_yyb;
        double T_xyb_xby = G_xy + G_xbyb - G_xxb - G_yyb;
        
        double a = G_xy + G_xbyb - T_xy_xbyb*Ncsqrm1_inverse;
        double b = T_xy_xbyb*twoCF_inverse;
        double c = T_xyb_xby*twoCF_inverse;
        double d = G_xyb + G_yxb - T_xyb_xby*Ncsqrm1_inverse;

        double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
        double DD;
        if (Sqrt == 0.0)
            DD = exp(0.5*(a+d));
        else
        {
            double factor = ( a-d+2.0/Nc*b ) / Sqrt;
            DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)*0.5 ) + (1.0-factor) * exp( (a+d-Sqrt)*0.5 ) );
        }

        return 4.0 * (1.0-D(G_xy)-D(G_xbyb)+DD);
#else
        return 4.0*(G_xy*G_xbyb+sqr(G_xxb+G_yyb-G_xyb-G_yxb)/2.0); // probably wrong (but also not in use but) needs to be checked (factor of 1/8 probably missing)
#endif
    }

    double dsigma_d2b_sqr_reduced (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xxb = NH*GBWModel::G(x1, x2, xb1, xb2);
        double G_xyb = NH*GBWModel::G(x1, x2, yb1, yb2);
        double G_yxb = NH*GBWModel::G(xb1, xb2, y1, y2);
        double G_yyb = NH*GBWModel::G(y1, y2, yb1, yb2);
#ifndef _DILUTE
        double G_xy = NH*GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = NH*GBWModel::G(xb1, xb2, yb1, yb2);

        double T_xy_xbyb = G_xyb + G_yxb - G_xxb - G_yyb;
        double T_xyb_xby = G_xy + G_xbyb - G_xxb - G_yyb;
        
        double a = G_xy + G_xbyb - T_xy_xbyb*Ncsqrm1_inverse;
        double b = T_xy_xbyb*twoCF_inverse;
        double c = T_xyb_xby*twoCF_inverse;
        double d = G_xyb + G_yxb - T_xyb_xby*Ncsqrm1_inverse;

        double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
        double DD;
        if (Sqrt == 0.0)
            DD = exp(0.5*(a+d));
        else
        {
            double factor = ( a-d+2.0/Nc*b ) / Sqrt;
            DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)*0.5 ) + (1.0-factor) * exp( (a+d-Sqrt)*0.5 ) );
        }

        return 4.0 * (DD - D(G_xy)*D(G_xbyb));
#else
        return 4.0*0.5*Ncsqrm1_inverse*sqr(G_xxb+G_yyb-G_xyb-G_yxb);
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
            double G_xy = 0.0;
            double G_xbyb = 0.0;
            double G_xxb = 0.0;
            double G_xyb = 0.0;
            double G_yxb = 0.0;
            double G_yyb = 0.0;

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

                    G_xy += GBWModel::G(x1_mod, x2_mod, y1_mod, y2_mod);
                    G_xbyb += GBWModel::G(xb1_mod, xb2_mod, yb1_mod, yb2_mod);
                    G_xxb += GBWModel::G(x1_mod, x2_mod, xb1_mod, xb2_mod);
                    G_xyb += GBWModel::G(x1_mod, x2_mod, yb1_mod, yb2_mod);
                    G_yxb += GBWModel::G(xb1_mod, xb2_mod, y1_mod, y2_mod);
                    G_yyb += GBWModel::G(y1_mod, y2_mod, yb1_mod, yb2_mod);
                }
    #ifndef _DILUTE
            double T_xy_xbyb = G_xyb + G_yxb - G_xxb - G_yyb;
            double T_xyb_xby = G_xy + G_xbyb - G_xxb - G_yyb;
            
            double a = G_xy + G_xbyb - T_xy_xbyb*Ncsqrm1_inverse;
            double b = T_xy_xbyb*twoCF_inverse;
            double c = T_xyb_xby*twoCF_inverse;
            double d = G_xyb + G_yxb - T_xyb_xby*Ncsqrm1_inverse;

            double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
            double DD;
            if (Sqrt == 0.0)
                DD = exp(0.5*(a+d));
            else
            {
                double factor = ( a-d+2.0*b ) / Sqrt;
                DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)*0.5 ) + (1.0-factor) * exp( (a+d-Sqrt)*0.5 ) );
            }
            
            return 4.0 * (1.0 - D(G_xy) - D(G_xbyb) + DD);
    #else
            return 0.5*(G_xy*G_xbyb + sqr(G_xxb+G_yyb-G_xyb-G_yxb)*0.5); // probably wrong (but also not in use but) needs to be checked (factor of 1/8 probably missing)
    #endif
        }


        double dsigma_d2b_sqr_reduced (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus)
        {
    #ifndef _DILUTE
            double G_xy = 0.0;
            double G_xbyb = 0.0;
    #endif
            double G_xxb = 0.0;
            double G_xyb = 0.0;
            double G_yxb = 0.0;
            double G_yyb = 0.0;

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
                    G_xy += GBWModel::G(x1_mod, x2_mod, y1_mod, y2_mod);
                    G_xbyb += GBWModel::G(xb1_mod, xb2_mod, yb1_mod, yb2_mod);
            #endif
                    G_xxb += GBWModel::G(x1_mod, x2_mod, xb1_mod, xb2_mod);
                    G_xyb += GBWModel::G(x1_mod, x2_mod, yb1_mod, yb2_mod);
                    G_yxb += GBWModel::G(y1_mod, y2_mod, xb1_mod, xb2_mod);
                    G_yyb += GBWModel::G(y1_mod, y2_mod, yb1_mod, yb2_mod);
                }
    #ifndef _DILUTE
            double T_xy_xbyb = G_xyb + G_yxb - G_xxb - G_yyb;
            double T_xyb_xby = G_xy + G_xbyb - G_xxb - G_yyb;
            
            double a = G_xy + G_xbyb - T_xy_xbyb*Ncsqrm1_inverse;
            double b = T_xy_xbyb*twoCF_inverse;
            double c = T_xyb_xby*twoCF_inverse;
            double d = G_xyb + G_yxb - T_xyb_xby*Ncsqrm1_inverse;

            double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
            double DD;
            if (Sqrt == 0.0)
                DD = exp(0.5*(a+d));
            else
            {
                double factor = ( a-d+2.0/Nc*b ) / Sqrt;
                DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)*0.5 ) + (1.0-factor) * exp( (a+d-Sqrt)*0.5 ) );
            }
            
            return 4.0 * (DD - D(G_xy)*D(G_xbyb));
    #else
            return 4.0*0.5*Ncsqrm1_inverse*sqr(G_xxb+G_yyb-G_xyb-G_yxb);
    #endif
        }
    }


    namespace HotspotAverage
    {
        void sample (uint A, uint H, uint num, uint start_seed)
        {
            HotspotAverage::num = num;
            inv_num = 1.0/double(num);

            hn = (HotspotNucleus*)malloc(sizeof(HotspotNucleus)*num);

            for (uint i=0; i<num; ++i)
                new (&hn[i]) HotspotNucleus(A, H, start_seed + i);
        }

        void clear()
        {
            for (uint i=0; i<HotspotAverage::num; ++i)
                hn[i].~HotspotNucleus();

            free(hn);
        }

        double dsigma_d2b (double x1, double x2, double y1, double y2)
        {
            double avg = 0.0;
            for (uint i=0; i<HotspotAverage::num; ++i)
                avg += Sampled::dsigma_d2b(x1, x2, y1, y2, &hn[i]);
            
            return avg * inv_num;
        }

        double dsigma_d2b_sqr_reduced (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double avg = 0.0;
            for (uint i=0; i<HotspotAverage::num; ++i)
                avg += Sampled::dsigma_d2b_sqr_reduced(x1, x2, y1, y2, xb1, xb2, yb1, yb2, &hn[i]);
            
            return avg * inv_num;
        }
    }
}