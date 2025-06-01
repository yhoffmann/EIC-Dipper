#include "../include/GBWModel.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "../include/IntegrationRoutines.hpp"
#include "../include/constants.hpp"
#include "../include/utilities.hpp"

namespace GBWModel {
double T_times_sigma0(double b1, double b2) {
  return exp(-(sqr(b1) + sqr(b2)) / (2.0 * rH_sqr));
}

double G_integrand_function(double u, double v, double x1, double x2, double y1,
                            double y2) {
  if (u == 0 && v == 0) return 0.0;

  double inverse_divisor = 1.0 / (u * v + rH_sqr / 2.0 * (u + v));

  return /*missing factor is in G()*/ exp(-sqr(m) * (u + v)) * inverse_divisor *
         (exp(-0.25 *
              (u * (sqr(y1) + sqr(y2)) + v * (sqr(x1) + sqr(x2)) +
               rH_sqr / 2.0 * (sqr(x1 - y1) + sqr(x2 - y2))) *
              inverse_divisor) -
          0.5 * exp(-0.25 * (sqr(x1) + sqr(x2)) * (u + v) * inverse_divisor) -
          0.5 * exp(-0.25 * (sqr(y1) + sqr(y2)) * (u + v) * inverse_divisor));
}

int G_integrand_cubature([[maybe_unused]] unsigned ndim, const double* xx,
                         void* userdata, [[maybe_unused]] unsigned fdim,
                         double* ff) {
  GIntegrandParams* params = (GIntegrandParams*)userdata;

  ff[0] = G_integrand_function(xx[0], xx[1], params->x1, params->x2, params->y1,
                               params->y2);

  return 0;
}

double G_by_integration(double x1, double x2, double y1, double y2) {
  CubatureConfig cubature_config;
  cubature_config.num_dims = 2;
  cubature_config.max_eval = 5e6;
  cubature_config.abs_err = 0.0;

  IntegrationConfig integration_config;
  GIntegrandParams G_integrand_params{x1, x2, y1, y2};
  integration_config.integrand_params = &G_integrand_params;

  integration_config.min = (double*)alloca(2 * sizeof(double));
  integration_config.max = (double*)alloca(2 * sizeof(double));

  integration_config.min[0] = 0.0;
  integration_config.max[0] = 20.0 / sqr(m);

  integration_config.min[1] = 0.0;
  integration_config.max[1] = integration_config.max[0];

  return /*CF * t_g2mu02 / (16.0*PI*PI) **/ IntegrationRoutines::
      cubature_integrate(G_integrand_cubature, &cubature_config,
                         &integration_config);
}

double G_wrapper(double r, double rb, double phi) {
  return G_by_integration(r, 0.0, rb * cos(phi), rb * sin(phi));
}

Interpolator3D G_ip;
constexpr const double CF_div_16pipi = CF / (16.0 * PI * PI);

double G(double x1, double x2, double y1, double y2) {
  double r = std::sqrt(sqr(x1) + sqr(x2));
  double rb = std::sqrt(sqr(y1) + sqr(y2));
  float arg = (r == 0.0 || rb == 0.0) ? 0.0 : (x1 * y1 + x2 * y2) / (r * rb);

  return CF_div_16pipi * t_g2mu02 * G_ip(r, rb, acos(arg));
}
}  // namespace GBWModel
