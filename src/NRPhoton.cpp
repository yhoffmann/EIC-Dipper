#include "../include/NRPhoton.hpp"

#include "../include/constants.hpp"
#include "../include/utilities.hpp"

#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdlib.h>

namespace NRPhoton {
double epsilon() { return std::sqrt(sqr(Q) * 0.25 + sqr(m_Q)); }

double wave_function_factor_T() {
  return -A_Q * std::sqrt(2.0 * m_Q * Nc) * e * e_Q;
}

double wave_function_factor_L() {
  return 0.5 * wave_function_factor_T() * Q / m_Q;
}

double wave_function_factor_combined() {
  return std::sqrt(sqr(wave_function_factor_T()) +
                   sqr(wave_function_factor_L()));
}

double wave_function(double r1, double r2) {
  double rsqr = sqr(r1) + sqr(r2);

  if (rsqr == 0.0) rsqr = 2.0e-15;

  return wave_function_factor_combined() *
         gsl_sf_bessel_K0(epsilon() * std::sqrt(rsqr));
}
}  // namespace NRPhoton
