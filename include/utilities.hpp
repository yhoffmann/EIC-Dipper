#pragma once


#define _LOG(x) std::cout << x << std::endl


#include <gsl/gsl_sf.h>
#include <string>
#include <random>
#include "../external/Interpolation3D/include/Interpolator3D.hpp"


typedef long unsigned int luint;


inline double sqr (double x)
{
    return x*x;
}


double bessel_K_safe (int n, double x);


void set_import_filepath_by_m (std::string& filepath, DataGenerationConfig* config);


inline std::string filepath_global;


void set_parameters(int argc, char** argv);


inline double rng01()
{
    return drand48();
}


inline double rng (double min, double max)
{
    return min+(max-min)*2.0*(drand48()-0.5);
}