#pragma once


#include <cuba.h>
#include "../external/Interpolation3D/include/Interpolator3D.hpp"


namespace GBWModel 
{
    extern Interpolator3D G_ip;
    
    double T_times_sigma0(double b1, double b2);

    double Q_s_sqr(double b1, double b2);

    double G_old(double x1, double x2, double y1, double y2);

    double Kmod(double x1, double x2, double y1, double y2);

    double G_mod(double x1, double x2, double y1, double y2);

    double G_integrand_function(double u, double v, double x1, double x2, double y1, double y2);

    int G_integrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);

    double G_by_integration(double x1, double x2, double y1, double y2);

    double G_wrapper(double r, double rb, double theta);

    double G(double x1, double x2, double y1, double y2);
}