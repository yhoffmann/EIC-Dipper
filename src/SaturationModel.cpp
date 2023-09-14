#include "../include/SaturationModel.hpp"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"

// #include "../../../libs/EigenLib/Eigen/Eigen"
// #include "../../../libs/EigenLib/unsupported/Eigen/MatrixFunctions"

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
        double G_xy = GBWModel::G(x1,x2,y1,y2);
        double G_xbyb = GBWModel::G(xb1,xb2,yb1,yb2);

        return 1.0+G_xy+G_xbyb;
    }


    double dsigma_d2b(double x1, double x2, double y1, double y2)
    {
        return -2.0*GBWModel::G(x1,x2,y1,y2);//2.0 * (1.0 - exp(GBWModel::G(x1, x2, y1, y2)));
    }


    double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
    {
        double G_xy = GBWModel::G(x1, x2, y1, y2);
        double G_xbyb = GBWModel::G(xb1, xb2, yb1, yb2);

        return 4.0 * (1.0 - D(G_xy) - D(G_xbyb) + DD(x1, x2, y1, y2, xb1, xb2, yb1, yb2));
    }

    namespace DDCorrelationMatrixElements
    {
        using namespace GBWModel;

        double T_xy_ybxb(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2)
        {
            double G_xxb = G(x1, x2, xb1, xb2);
            double G_yyb = G(y1, y2, yb1, yb2);
            double G_xyb = G(x1, x2, yb1, yb2);
            double G_ybx = G(y1, y2, xb1, xb2);

            return G_xxb + G_yyb - G_xyb - G_ybx;
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