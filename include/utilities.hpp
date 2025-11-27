#pragma once
#ifndef EIC_UTILITIES_HPP_
#define EIC_UTILITIES_HPP_

#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "constants.hpp"

#include <chrono>
#include <gsl/gsl_sf.h>
#include <mutex>
#include <string>

#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#define GET_PROCESS_ID() getpid()

#elif (defined(_WIN32) || defined(_WIN64))
#include <windows.h>
#define GET_PROCESS_ID() GetCurrentProcessId()

#else
#error Unkown operating system: Unable to define GET_PROCESS_ID() function

#endif

// error definitions
#define EIC_ERROR_COULDNT_OPEN 20
#define EIC_ERROR_NO_INTERP_DATA 21
#define EIC_ERROR_INTERP_USER_NO 22
#define EIC_ERROR_UNKNOWN_OPT 23
#define EIC_ERROR_BAD_ALLOC 24

inline std::chrono::high_resolution_clock::time_point g_time_program_start;
double get_time_ms_since_program_start();
#ifdef _TEST
#define TEST_LOG(x)                                                          \
  {                                                                          \
    std::cerr << "(TEST LOG " << get_time_ms_since_program_start() << "ms) " \
              << x << std::endl;                                             \
  }
#else
#define TEST_LOG(x) ;
#endif

inline std::mutex cout_mutex;

void init(int argc, char* argv[]);

inline Interpolator3D::DataGenerationConfig default_data_generation_config{
    .nx = 400,
    .x_max = 8.0 / m,
    .ny = 400,
    .y_max = 8.0 / m,
    .nz = 60,
    .z_exp_grid_spacing_parameter = 0.0};

enum ProgressLogLevel : size_t { None = 0, Converged, All };

inline ProgressLogLevel g_progress_log_level = Converged;
inline double g_Delta_single = 1.0;
inline bool g_Delta_single_set = false;

inline uint g_seed = 0;

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
