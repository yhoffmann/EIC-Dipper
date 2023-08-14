#pragma once


#include <gsl/gsl_sf.h>
#include <string>


typedef long unsigned int luint;


inline double sqr (double x)
{
    return x*x;
}


double bessel_K_safe (int n, double x);


void set_import_filepath_by_m (std::string& filepath);