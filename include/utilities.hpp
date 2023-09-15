#pragma once


#define LOG(x) std::cout << x << std::endl


#include <gsl/gsl_sf.h>
#include <string>
#include "../external/Interpolation3D/include/Interpolator3D.hpp"


typedef long unsigned int luint;


inline double sqr (double x)
{
    return x*x;
}


double bessel_K_safe (int n, double x);


void set_import_filepath_by_m (std::string& filepath, DataGenerationConfig* config);