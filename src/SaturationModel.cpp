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
    double D(double G)
    {
        return exp(G);
    }


    double D_dilute (double G)
    {
        return 1.0+G;
    }


    double DD(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double a = DDCorrelationMatrixElements::a(x1, x2, y1, y2, xb1, xb2, yb1, yb2);
        double b = DDCorrelationMatrixElements::b(x1, x2, y1, y2, xb1, xb2, yb1, yb2);
        double c = DDCorrelationMatrixElements::c(x1, x2, y1, y2, xb1, xb2, yb1, yb2);
        double d = DDCorrelationMatrixElements::d(x1, x2, y1, y2, xb1, xb2, yb1, yb2);

        double Sqrt = std::sqrt(4.0 * b * c + sqr(a - d));
        double factor = (a - d + 2.0 / Nc * b) / Sqrt;

        double ret = 0.5 * ((1.0 + factor) * exp((a + d + Sqrt) / 2.0) + (1.0 - factor) * exp((a + d - Sqrt) / 2.0));

        return ret;
    }


    double DD_dilute (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xy = GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = GBWModel::G(xb1, xb2, yb1, yb2);

        return 1.0+G_xy+G_xbyb;
    }


    double dsigma_d2b(double x1, double x2, double y1, double y2)
    {
        return -2.0*GBWModel::G(x1, x2, y1, y2);//2.0 * (1.0 - exp(GBWModel::G(x1, x2, y1, y2)));//
    }


    double dsigma_d2b_sqr_old (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xy = GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = GBWModel::G(xb1, xb2, yb1, yb2);

        return 4.0 * (1.0 - D(G_xy) - D(G_xbyb) + DD(x1, x2, y1, y2, xb1, xb2, yb1, yb2));
    }

    double dsigma_d2b_sqr (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xxb = GBWModel::G(x1, x2, xb1, xb2);
        double G_xy = GBWModel::G(x1, x2, y1, y2);
        double G_xyb = GBWModel::G(x1, x2, yb1, yb2);
        double G_xby = GBWModel::G(xb1, xb2, y1, y2);
        double G_xbyb = GBWModel::G(xb1, xb2, yb1, yb2);
        double G_yyb = GBWModel::G(y1, y2, yb1, yb2);

        double T_xy_ybxb = G_xxb + G_yyb - G_xyb - G_xby;
        double T_xxb_yby = G_xy + G_xbyb - G_xyb - G_xby;
        
        double a = G_xy + G_xbyb - T_xy_ybxb / (sqr(Nc)-1.0);
        double b = T_xy_ybxb / (2.0*CF);
        double c = T_xxb_yby / (2.0*CF);
        double d = G_xxb + G_yyb - T_xxb_yby / (sqr(Nc)-1.0);

        double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
        double factor = ( a-d+2.0*b ) / Sqrt;

        double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );

        return 4.0 * ( 1.0-exp(G_xy)-exp(G_xbyb)+DD );
    }


    namespace GeometryAverage
    {
        double dsigma_d2b (double x1, double x2, double y1, double y2, const HotspotNucleus* nucleus)
        {
            double G_sum = 0.0;
            for (uint n=0, a=nucleus->get_atomic_num(); n<a; ++n)
                for (uint i=0, nh=nucleus->get_num_hotspots_per_nucleon(); i<nh; ++i)
                {
                    const double* B0_ptr = nucleus->get_hotspot_pos(n, i);

                    double B0[2] = {B0_ptr[0], B0_ptr[1]};

                    G_sum += GBWModel::G(x1-B0[0], x2-B0[1], y1-B0[0], y2-B0[1]);
                }

            return 2.0 * ( 1.0-exp(G_sum) );//-2.0*G_sum;//
        }


        double dsigma_d2b_sqr (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus)
        {
            double G_xxb_sum = 0.0;
            double G_xy_sum = 0.0;
            double G_xyb_sum = 0.0;
            double G_xby_sum = 0.0;
            double G_xbyb_sum = 0.0;
            double G_yyb_sum = 0.0;

            for (uint n=0, a=nucleus->get_atomic_num(); n<a; ++n)
                for (uint i=0, nh=nucleus->get_num_hotspots_per_nucleon(); i<nh; ++i)
                {
                    const double* B0_ptr = nucleus->get_hotspot_pos(n, i);

                    double B0[2] = {B0_ptr[0], B0_ptr[1]};

                    double x1_mod = x1-B0[0];
                    double x2_mod = x2-B0[1];
                    double y1_mod = y1-B0[0];
                    double y2_mod = y2-B0[1];
                    double xb1_mod = xb1-B0[0];
                    double xb2_mod = xb2-B0[1];
                    double yb1_mod = yb1-B0[0];
                    double yb2_mod = yb2-B0[1];

                    G_xxb_sum += GBWModel::G(x1_mod, x2_mod, xb1_mod, xb2_mod);
                    G_xy_sum += GBWModel::G(x1_mod, x2_mod, y1_mod, y2_mod);
                    G_xyb_sum += GBWModel::G(x1_mod, x2_mod, yb1_mod, yb2_mod);
                    G_xby_sum += GBWModel::G(xb1_mod, xb2_mod, y1_mod, y2_mod);
                    G_xbyb_sum += GBWModel::G(xb1_mod, xb2_mod, yb1_mod, yb2_mod);
                    G_yyb_sum += GBWModel::G(y1_mod, y2_mod, yb1_mod, yb2_mod);
                }
            
            double T_xy_ybxb = G_xxb_sum + G_yyb_sum - G_xyb_sum - G_xby_sum;
            double T_xxb_yby = G_xy_sum + G_xbyb_sum - G_xyb_sum - G_xby_sum;
            
            double a = G_xy_sum + G_xbyb_sum - T_xy_ybxb / (sqr(Nc)-1.0);
            double b = T_xy_ybxb / (2.0*CF);
            double c = T_xxb_yby / (2.0*CF);
            double d = G_xxb_sum + G_yyb_sum - T_xxb_yby / (sqr(Nc)-1.0);

            double Sqrt = std::sqrt( 4.0*b*c + sqr(a-d) );
            double factor = ( a-d+2.0*b ) / Sqrt;

            double DD = 0.5 * ( (1.0+factor) * exp( (a+d+Sqrt)/2.0 ) + (1.0-factor) * exp( (a+d-Sqrt)/2.0 ) );

            return 4.0 * ( 1.0-exp(G_xy_sum)-exp(G_xbyb_sum)+DD );
        }
    }


    namespace DDCorrelationMatrixElements
    {
        using namespace GBWModel;

        double T_xy_ybxb(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double G_xxb = G(x1, x2, xb1, xb2);
            double G_yyb = G(y1, y2, yb1, yb2);
            double G_xyb = G(x1, x2, yb1, yb2);
            double G_yxb = G(y1, y2, xb1, xb2);

            return G_xxb + G_yyb - G_xyb - G_yxb;
        }

        double T_xxb_yby(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double G_xy = G(x1, x2, y1, y2);
            double G_xbyb = G(xb1, xb2, yb1, yb2);
            double G_xyb = G(x1, x2, yb1, yb2);
            double G_ybx = G(y1, y2, xb1, xb2);

            return G_xy + G_xbyb - G_xyb - G_ybx;
        }

        double a(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double G_xy = G(x1, x2, y1, y2);
            double G_ybxb = G(xb1, xb2, yb1, yb2);
            double T = T_xy_ybxb(x1, x2, y1, y2, xb1, xb2, yb1, yb2);

            return G_xy + G_ybxb - T / (sqr(Nc) - 1.0);
        }

        double b(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double T = T_xy_ybxb(x1, x2, y1, y2, xb1, xb2, yb1, yb2);

            return T / (2.0 * CF);
        }

        double c(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double T = T_xxb_yby(x1, x2, y1, y2, xb1, xb2, yb1, yb2);

            return T / (2.0 * CF);
        }

        double d(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double G_xxb = G(x1, x2, xb1, xb2);
            double G_yyb = G(y1, y2, yb1, yb2);
            double T = T_xxb_yby(x1, x2, y1, y2, xb1, xb2, yb1, yb2);

            return G_xxb + G_yyb - T / (sqr(Nc) - 1.0);
        }
    }
}