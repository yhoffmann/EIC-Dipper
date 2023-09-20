#pragma once

namespace SaturationModel
{
    double D(double G);

    double D_dilute(double G);

    double DD(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    double DD_dilute(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    double dsigma_d2b(double x1, double x2, double y1, double y2);

    namespace DDCorrelationMatrixElements
    {

        double T_xy_ybxb(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double T_xxb_yby(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double a(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double b(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double c(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double d(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);
    }
}