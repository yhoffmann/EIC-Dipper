#pragma once


#define _COUT(x) std::cout << x << std::endl


#include <gsl/gsl_sf.h>
#include <string>
#include <random>
#include <thread>
#include "constants.hpp"
#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


inline DataGenerationConfig default_data_generation_config
{
    .nx = 400,
    .x_max = 8.0/m,
    .ny = 400,
    .y_max = 8.0/m,
    .nz = 60,
    .z_exp_grid_spacing_parameter = 0.0
};

inline bool progress_monitor_global = false;

inline uint seed = 147541768;

typedef long unsigned int luint;

double sqr(double x);

double bessel_K_safe(int n, double x);

void import_interp_data_by_params(std::string& filepath, const DataGenerationConfig* config = &default_data_generation_config);

inline std::string filepath_global = "";

void set_parameters(int argc, char** argv);

// inline std::mt19937 rng;
// inline std::uniform_real_distribution<double> dist_01 = std::uniform_real_distribution<double>();

void print_infos(std::ofstream& out, uint seed);

void print_infos(std::ofstream& out, uint seed, const HotspotNucleus& nucleus);

uint get_unique_seed();

enum class Quark : char
{
    Charm = 'c', Bottom = 'b'
};
inline Quark quark_config = Quark::Charm;

std::string get_default_filepath_from_parameters();

double sin_zeros(uint n);
double cos_zeros(uint n);