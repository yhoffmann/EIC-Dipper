#pragma once

namespace SaturationModel {
    double D (double x);

    double DD (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);


    double dsigma_d2b_sqr (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

    double dsigma_d2b (double b1, double b2, double r1, double r2);

    /*namespace CorrelationMatrixElements {

        double G_xy (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double G_xyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double G_xby (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double G_xbyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double G_xxb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double G_yyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double T1 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double T2 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double a (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double b (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double c (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);

        double d (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2);
    }*/

    namespace DDCorrelationMatrixElements {
        double T_xy_ybxb (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double T_xxb_yby (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double a (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double b (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double c (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double d (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);
    }
}