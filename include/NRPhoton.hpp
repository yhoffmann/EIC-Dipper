#pragma once
#ifndef EIC_NRPHOTON_CPP_
#define EIC_NRPHOTON_CPP_

namespace NRPhoton {
double epsilon();

double wave_function_factor_T();

double wave_function_factor_L();

double wave_function_factor_combined();

double wave_function(double r1, double r2);
}  // namespace NRPhoton

#endif  // EIC_NRPHOTON_CPP_
