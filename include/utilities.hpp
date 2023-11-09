#pragma once


#define _COUT(x) std::cout << x << std::endl


#include <gsl/gsl_sf.h>
#include <string>
#include <random>
#include <thread>
#include "constants.hpp"
#include "../external/Interpolation3D/include/Interpolator3D.hpp"


const DataGenerationConfig default_data_generation_config { .x_max = 8.0/m, .y_max = 8.0/m, .n_z = 40 };


inline bool progress_monitor_global = false;


typedef long unsigned int luint;


double sqr (double x);


double bessel_K_safe (int n, double x);


void set_import_filepath_by_parameters (std::string& filepath, const DataGenerationConfig* config = &default_data_generation_config);


inline std::string filepath_global = "";


void set_parameters (int argc, char** argv);


inline std::mt19937 rng;
inline std::uniform_real_distribution<double> dist_01 = std::uniform_real_distribution<double>();


void print_infos (std::ofstream& out);