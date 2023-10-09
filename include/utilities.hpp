#pragma once


#define _COUT(x) std::cout << x << std::endl


#include <gsl/gsl_sf.h>
#include <string>
#include <random>
#include <thread>
#include "../external/Interpolation3D/include/Interpolator3D.hpp"


typedef long unsigned int luint;


double sqr (double x);


double bessel_K_safe (int n, double x);


void set_import_filepath_by_m (std::string& filepath, DataGenerationConfig* config);


inline std::string filepath_global;


void set_parameters(int argc, char** argv);


inline std::mt19937 rng;
inline std::uniform_real_distribution<double> dist_01 = std::uniform_real_distribution<double>();


void create_info_file (const std::string& filepath);