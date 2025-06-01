#pragma once
#ifndef EIC_GBWMODEL_HPP_
#define EIC_GBWMODEL_HPP_

#include "../external/Interpolation3D/include/Interpolator3D.hpp"

inline std::string interpolator_filepath = "";

namespace GBWModel {
extern Interpolator3D G_ip;

double T_times_sigma0(double b1, double b2);

double Q_s_sqr(double b1, double b2);

double G_old(double x1, double x2, double y1, double y2);
double Kmod(double x1, double x2, double y1, double y2);
double G_mod(double x1, double x2, double y1, double y2);

double G_integrand_function(double u, double v, double x1, double x2, double y1,
                            double y2);
double G_by_integration(double x1, double x2, double y1, double y2);
double G_wrapper(double r, double rb, double phi);

double G(double x1, double x2, double y1, double y2);
}  // namespace GBWModel

#endif  // EIC_GBWMODEL_HPP_
