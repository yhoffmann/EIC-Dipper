#pragma once


inline double x (double b, double r) {
    return b+r/2.0;
}

inline double y (double b, double r) {
    return b-r/2.0;
}

namespace GBWModel {
    double T_times_sigma0 (double b1, double b2);

    double Q_s_sqr (double b1, double b2);

    double G (double x1, double x2, double y1, double y2);
}

namespace NRPhoton {
    double epsilon (double Q, double z);

    double wave_function (double r1, double r2, double Q, double z, bool t_not_l);
}

namespace CorrelationMatrixElements{
    using namespace GBWModel;

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
}