#pragma once
#ifndef EIC_UTILITIES_HPP_
#define EIC_UTILITIES_HPP_

#include <gsl/gsl_sf.h>

#include <chrono>
#include <string>

#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "constants.hpp"

#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#define GET_PROCESS_ID() getpid()

#elif (defined(_WIN32) || defined(_WIN64))
#include <windows.h>
#define GET_PROCESS_ID() GetCurrentProcessId()

#else
#error Unkown operating system: Unable to define GET_PROCESS_ID() function

#endif

inline std::chrono::high_resolution_clock::time_point g_time_program_start;

#ifdef _TEST
#define TEST_LOG(x)                                                         \
  {                                                                         \
    auto now = std::chrono::high_resolution_clock::now();                   \
    double time_ms = std::chrono::duration_cast<std::chrono::milliseconds>( \
                         now - g_time_program_start)                        \
                         .count();                                          \
    std::cerr << "(TEST LOG " << time_ms << "ms) " << x << std::endl;       \
  }
#else
#define TEST_LOG(x) ;
#endif

void init(int argc, char* argv[]);

inline Interpolator3D::DataGenerationConfig default_data_generation_config{
    .nx = 400,
    .x_max = 8.0 / m,
    .ny = 400,
    .y_max = 8.0 / m,
    .nz = 60,
    .z_exp_grid_spacing_parameter = 0.0};

inline bool g_monitor_progress = false;
inline double g_Delta_single = 1.0;
inline bool g_Delta_single_set = false;

inline uint g_seed = 147541768;

inline size_t g_num_threads = 10;

typedef long unsigned int luint;

inline double sqr(double x) { return x * x; }

inline double bessel_K_safe(int n, double x) {
  if (x == 0) x = 1.0e-20;

  return gsl_sf_bessel_Kn(n, x);
}

void import_interp_data_by_params(const Interpolator3D::DataGenerationConfig*
                                      config = &default_data_generation_config);

inline std::string g_filepath = "";

void set_parameters(int argc, char** argv);

void print_infos(std::ofstream& out, uint seed);

void print_infos(std::ofstream& out, uint seed, const HotspotNucleus& nucleus);

uint get_unique_seed();

enum class Quark : char { Charm = 'c', Bottom = 'b' };
inline Quark quark_output_config = Quark::Charm;
std::string get_default_filepath_from_parameters();

double sin_zeros(uint n);
double cos_zeros(uint n);

#endif  // EIC_UTILITIES_HPP_
