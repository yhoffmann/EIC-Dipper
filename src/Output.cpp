#include "../include/Output.hpp"

#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "../external/thread-pool/include/ThreadPool.hpp"
#include "../include/Coherent.hpp"
#include "../include/DipoleModel.hpp"
#include "../include/Incoherent.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/constants.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

/*#define _STD_ARRAY(type, name, ...) \
    constexpr size_t name##_size = sizeof( (type[]){ __VA_ARGS__
})/sizeof(type); \ std::array<type, name##_size> name( { __VA_ARGS__ } );

#define _GET_STD_ARRAY(type, name) \
    std::array<type, name##_size> get_##name()

_STD_ARRAY(double, default_Q_vec, 0.001, 0.01, 0.04, 0.08, 0.12, 0.2, 0.3, 0.4,
0.5, 0.6, 0.7, 0.8,
0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.3, 3.0);
_GET_STD_ARRAY(double, default_Q_vec)
{
    return default_Q_vec;
}*/

namespace Output {
inline double get_default_Q() { return std::sqrt(0.1); }

std::vector<double> get_default_Delta_vec() {
  return std::vector<
      double>{0.001, 0.01, 0.04, 0.08, 0.12, 0.2, 0.3, 0.4, 0.5,
              0.6,   0.7,  0.8,  0.9,  1.0,  1.1, 1.2, 1.3, 1.4,
              1.5,   1.6,  1.7,  1.8,  2.0,  2.3, 3.0};  //{0.001, 0.002,
                                                         // 0.005, 0.007,
                                                         // 0.01, 0.03,
                                                         // 0.05, 0.07,
                                                         // 0.08, 0.09,
                                                         // 0.12, 0.16, 0.2,
                                                         // 0.3, 0.4, 0.5,
                                                         // 0.6, 0.7, 0.8,
                                                         // 0.9, 1.0, 1.2, 1.4,
                                                         // 1.6,
                                                         // 1.8, 2.0, 2.3, 2.7,
                                                         // 3.0,
                                                         // 3.3, 3.7, 4.0};
                                                         // //
}

std::vector<double> get_default_g2mu02_factor_vec() {
  return std::vector<double>{0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                             0.9,  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                             1.8,  1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
}

std::vector<double> get_default_phi_vec() {
  constexpr uint num_angles = 48;
  std::vector<double> default_phi_vec;
  default_phi_vec.reserve(num_angles);
  for (uint i = 0; i < num_angles; i++)
    default_phi_vec.push_back(
        PI * double(i) /
        double(num_angles));  // no num_angles-1 because we do not want to reach
                              // the end // NOTE if the calculation of the
                              // coherent cross section ever changes to not use
                              // sin and cos anymore but something more
                              // complicated, this needs to be adjusted to 2pi
                              // again and also the writing of data to file
                              // needs to be changed (currently I am using that
                              // cos and sin are symmetric/antisymmetric to
                              // print 2 results for one calculation/angle)

  return default_phi_vec;
}

void dsigmadt(bool do_coherent, bool do_incoherent) {
  std::vector<double> default_Delta_vec = get_default_Delta_vec();
  std::vector<double> default_phi_vec =
      std::vector<double>{0.0};  // get_default_phi_vec();//

  dsigmadt(do_coherent, do_incoherent, default_Delta_vec, default_phi_vec);
}

void dsigmadt([[maybe_unused]] bool do_coherent,
              [[maybe_unused]] bool do_incoherent,
              std::vector<double> value_vec, std::vector<double> phi_vec) {
  TEST_LOG(
      "In function Output::dsigmadt(bool, bool, double, std::vector<double>, "
      "std::vector<double>)")
  if (phi_vec.size() == 0) phi_vec = std::vector<double>{0.0};

  std::vector<std::vector<double>> coherent_results(
      value_vec.size(), std::vector<double>(phi_vec.size()));
  std::vector<double> coherent_avg_results(value_vec.size());

  std::vector<double> incoherent_results(value_vec.size());
  TEST_LOG("Starting ThreadPool")
  ThreadPool pool(g_num_threads);
  TEST_LOG("Queueing incoherent jobs")
  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
#ifndef _G2MU02
    double Delta = value_vec[value_index];
#else
    double Delta = g_Delta_single;
#endif
    pool.enq_job([value_index, &value_vec, Delta, &incoherent_results] {
#ifdef _G2MU02
      t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
      incoherent_results[value_index] = Incoherent::dsigmadt(Delta);
    });
  }
  TEST_LOG("Queueing coherent jobs")
  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
#ifndef _G2MU02
    double Delta = value_vec[value_index];
#else
    double Delta = g_Delta_single;
#endif
    for (uint phi_index = 0, phi_size = phi_vec.size(); phi_index < phi_size;
         ++phi_index) {
      double phi = phi_vec[phi_index];
      pool.enq_job(
          [value_index, &value_vec, phi_index, Delta, phi, &coherent_results] {
#ifdef _G2MU02
            t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
            coherent_results[value_index][phi_index] =
                Coherent::dsigmadt_test(Delta, phi);
          });
    }
  }
  TEST_LOG("Finished queueing, waiting for jobs to finish")
  pool.await();
  pool.stop();

  for (uint value_index = 0; value_index < value_vec.size(); ++value_index) {
    double coherent_avg = 0.0;

    for (uint phi_index = 0; phi_index < phi_vec.size(); ++phi_index) {
      coherent_avg += coherent_results[value_index][phi_index];
    }
    coherent_avg_results[value_index] = coherent_avg / double(phi_vec.size());
  }

  std::ofstream out;
  std::string filepath = g_filepath;

  if (filepath == "") {
    filepath = "data/dsigmadt_avg.dat";
  }
  out.open(filepath);

  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);

  out << std::setprecision(10);
  out << "#Q, Delta, Coher, Incoher; " << std::endl;
  out << "#1, 2,           3,4" << std::endl;

  for (uint i = 0; i < value_vec.size(); i++)
    out << value_vec[i] << " " << Q << " " << coherent_avg_results[i] << " "
        << incoherent_results[i] << std::endl;
  out << std::endl;

  out.flush();
  out.close();

  out.open(filepath + ".all");
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);

  for (uint i = 0; i < value_vec.size(); i++) {
    out << value_vec[i] << " " << Q << "   ";
    for (uint j = 0, jmax = phi_vec.size(); j < jmax; j++)
      out << phi_vec[j] << " " << coherent_results[i][j] << " "
          << incoherent_results[i] << "   ";
    out << std::endl;
  }
  out << std::endl;
  TEST_LOG(
      "Returning from function Output::dsigmadt(bool, bool, double, "
      "std::vector<double>, std::vector<double>)")
}

void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed) {
  std::vector<double> value_vec;
  std::vector<double> phi_vec = get_default_phi_vec();

#ifndef _G2MU02
  if (!g_Delta_single_set)
    value_vec = get_default_Delta_vec();

  else
    value_vec = {g_Delta_single};
#else
  value_vec = get_default_g2mu02_factor_vec();
#endif

  dsigmadt_nucleus(atomic_num, num_hotspots, seed, value_vec, phi_vec);
}

void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed,
                      std::vector<double> value_vec,
                      std::vector<double> phi_vec) {
  TEST_LOG(
      "In function Output::dsigmadt_nucleus(uint, uint, uint, double, "
      "std::vector<double>, std::vector<double>)")
  std::vector<std::vector<double>> coherent_results_real(
      value_vec.size(), std::vector<double>(phi_vec.size()));
  std::vector<std::vector<double>> coherent_results_imag(
      value_vec.size(), std::vector<double>(phi_vec.size()));
  std::vector<double> incoherent_results_real(value_vec.size());
  std::vector<double> incoherent_results_imag(value_vec.size());

  if (seed == 0) seed = get_unique_seed();

  HotspotNucleus nucleus(seed, atomic_num, num_hotspots, std::sqrt(R_sqr),
                         std::sqrt(rH_sqr));
  /* uncomment this line for legacy mode (a change to the
   * Nucleus class removed the need for constructing, then setting the size,
   * then sampling again to have the changed size take effect; only one
   * construction with correct parameters possible now; this means that all old
   * data will contain twice-sampled Nucleus objects, and event IDs are not
   * compatible anymore)
   */
  // nucleus.sample();
  // nucleus.sample_hotspot_weights();

  TEST_LOG("Starting ThreadPool")
  ThreadPool pool(g_num_threads);
  TEST_LOG("Queueing incoherent jobs")
  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
#ifndef _G2MU02
    double Delta = value_vec[value_index];
#else
    double Delta = g_Delta_single;
#endif
    pool.enq_job([value_index, &value_vec, Delta, &nucleus,
                  &incoherent_results_real, &incoherent_results_imag] {
#ifdef _G2MU02
      t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
      incoherent_results_real[value_index] =
          Incoherent::Sampled::dsigmadt_single_event(Delta, nucleus);
      incoherent_results_imag[value_index] = 0.0;
    });
  }
  TEST_LOG("Queueing coherent jobs")
  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
#ifndef _G2MU02
    double Delta = value_vec[value_index];
#else
    double Delta = g_Delta_single;
#endif
    for (uint phi_index = 0, phi_size = phi_vec.size(); phi_index < phi_size;
         ++phi_index) {
      double phi = phi_vec[phi_index];
      pool.enq_job([value_index, &value_vec, phi_index, Delta, phi, &nucleus,
                    &coherent_results_real, &coherent_results_imag] {
#ifdef _G2MU02
        t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
        auto [coh_real, coh_imag] =
            Coherent::Sampled::sqrt_dsigmadt_single_event(Delta, phi, nucleus);

        coherent_results_real[value_index][phi_index] = coh_real;
        coherent_results_imag[value_index][phi_index] = coh_imag;
      });
    }
  }
  TEST_LOG("Finished queueing, waiting for jobs to finish")
  pool.await();
  pool.stop();

  std::string filepath =
      (g_filepath != "")
          ? g_filepath
          : get_default_filepath_from_parameters() + "_Amplitude.dat";
  TEST_LOG("Opening file " << filepath)
  std::ofstream out(filepath);
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);
  TEST_LOG("File opened")
  TEST_LOG("Printing to file")
  out << std::setprecision(10);
  print_infos(out, seed, nucleus);

  out << "##Delta,   Q,        A Co real,Co imag,  A2 Inco\n";

  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
    out << value_vec[value_index] << " " << Q << "   ";
    for (uint phi_index = 0, phi_size = phi_vec.size(); phi_index < phi_size;
         ++phi_index)
      out << phi_vec[phi_index] << " "
          << coherent_results_real[value_index][phi_index] << " "
          << coherent_results_imag[value_index][phi_index] << " "
          << incoherent_results_real[value_index] << " "
          << incoherent_results_imag[value_index] << "   ";
    for (uint phi_index = 0, phi_size = phi_vec.size(); phi_index < phi_size;
         ++phi_index)
      out << phi_vec[phi_index] + PI << " "
          << -coherent_results_real[value_index][phi_index] << " "
          << coherent_results_imag[value_index][phi_index] << " "
          << incoherent_results_real[value_index] << " "
          << incoherent_results_imag[value_index] << "   ";
    out << std::endl;
  }
  out << std::endl;
  TEST_LOG("Finished priting to file")
  out.close();
  TEST_LOG(
      "Returning from function Output::dsigmadt_nucleus(uint, uint, uint, "
      "double, std::vector<double>, std::vector<double>)")
}

void dsdt_nucleus_internal_avg(uint atomic_num, uint num_hotspots, uint seed,
                               std::vector<double> value_vec) {
  std::vector<double> dsdt_sch_sbch(value_vec.size());
  std::vector<double> dsdt_ssbch(value_vec.size());
  std::vector<double> dsdt_scsbch(value_vec.size());

  if (0 == seed) seed = get_unique_seed();

  SaturationModel::InternalHotspotAvg::init(atomic_num, num_hotspots, 256,
                                            seed);

  ThreadPool pool(g_num_threads);

  // enqueueing <<sigma>c>h <<sigmabar>c>h (coherent) jobs
  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
#ifndef _G2MU02
    double Delta = value_vec[value_index];
#else
    double Delta = g_Delta_single;
#endif
    // enqueueing  <<sigma>c>h <<sigmabar>c>h  (coherent) jobs
    pool.enq_job([value_index, &value_vec, Delta, &dsdt_sch_sbch] {
#ifdef _G2MU02
      t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
      dsdt_sch_sbch[value_index] =
          Coherent::InternalHotspotAvg::dsdt_sch_sbch(Delta);
    });

    // enqueueing  < <sigma sigmabar>c >h  (color fluc first term) jobs
    pool.enq_job([value_index, &value_vec, Delta, &dsdt_ssbch] {
#ifdef _G2MU02
      t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
      dsdt_ssbch[value_index] =
          Incoherent::InternalHotspotAvg::dsdt_ssbch(Delta);
    });

    // enqueueing  < <sigma>c <sigmabar>c >h  (hotspot fluc first term) jobs
    pool.enq_job([value_index, &value_vec, Delta, &dsdt_scsbch] {
#ifdef _G2MU02
      t_g2mu02 = G2MU02_DEMIRCI * value_vec[value_index];
#endif
      dsdt_scsbch[value_index] =
          Incoherent::InternalHotspotAvg::dsdt_scsbch(Delta);
    });
  }
  pool.await();
  pool.stop();

  for (uint value_index = 0, value_size = value_vec.size();
       value_index < value_size; ++value_index) {
    std::cout << dsdt_sch_sbch[value_index] << " " << dsdt_ssbch[value_index]
              << " " << dsdt_scsbch[value_index] << "\n";
  }
}

void dsdt_nucleus_internal_avg(uint atomic_num, uint num_hotspots, uint seed) {
  std::vector<double> value_vec;
  std::vector<double> phi_vec = get_default_phi_vec();

#ifndef _G2MU02
  if (!g_Delta_single_set)
    value_vec = get_default_Delta_vec();

  else
    value_vec = {g_Delta_single};
#else
  value_vec = get_default_g2mu02_factor_vec();
#endif

  dsdt_nucleus_internal_avg(atomic_num, num_hotspots, seed, value_vec);
}

void dsigmadt_demirci(std::string filepath) {
  constexpr uint Delta_vec_size = 50;
  std::vector<double> Delta_vec;
  Delta_vec.reserve(Delta_vec_size);
  for (uint i = 0; i < Delta_vec_size; ++i)
    Delta_vec.push_back(
        std::sqrt(25.0) * double(i) / double(Delta_vec_size - 1) + 0.001);

  double coherent[Delta_vec_size];
  double color_fluc[Delta_vec_size];
  double hotspot_fluc[Delta_vec_size];

#pragma omp parallel for
  for (uint i = 0; i < Delta_vec_size; ++i) {
    coherent[i] = Coherent::Demirci::dsigmadt(Delta_vec[i]);
    color_fluc[i] = Incoherent::Demirci::color_fluctuations(Delta_vec[i]);
    hotspot_fluc[i] = Incoherent::Demirci::hotspot_fluctuations(Delta_vec[i]);

    if (g_progress_log_level) std::cout << i << std::endl;
  }

  std::ofstream out(filepath);
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);

  out << "t Co Inco Color Hotspot" << std::endl;
  for (uint i = 0; i < Delta_vec_size; ++i)
    out << sqr(Delta_vec[i]) << " " << coherent[i] << " "
        << color_fluc[i] + hotspot_fluc[i] << " " << color_fluc[i] << " "
        << hotspot_fluc[i] << std::endl;

  out.close();
}

// void dsdt_nucleus_avg_test (uint atomic_num, uint num_hotspots, uint seed,
// std::vector<double> Delta_vec, std::vector<double> phi_vec)
// {
//     DipoleModel::G_ip.import_data(interpolator_filepath);

//     SaturationModel::HotspotAverage::sample(A, H, 32, time(0));

//     ThreadPool pool(num_threads);

//     SaturationModel::HotspotAverage::clear();
// }

void G(uint num_points, std::string filepath) {
  const uint num_angles = 8;
  std::vector<double> angles(num_angles);
  for (uint i = 0; i < num_angles; ++i)
    angles[i] = double(i) / double(num_angles) * PI;

  double r = 1.0;

  std::vector<double> b_vec(num_points);

  std::vector<std::vector<std::array<double, 2>>> results(
      num_points, std::vector<std::array<double, 2>>(num_angles));

  // #pragma omp parallel for schedule(dynamic, 1)
  for (uint i = 0; i < num_points; ++i) {
    double x1, x2, y1, y2;
    double b = 100.0 * double(i) / double(num_points - 1);

    b_vec[i] = b;
    for (uint j = 0; j < num_angles; ++j) {
      double b1 = 0.0;
      double b2 = b;

      double r1 = r * cos(angles[j]);
      double r2 = r * sin(angles[j]);

      x1 = b1 + 0.5 * r1;
      x2 = b2 + 0.5 * r2;
      y1 = b1 - 0.5 * r1;
      y2 = b2 - 0.5 * r2;

      results[i][j][0] = DipoleModel::G(x1, x2, y1, y2);
      results[i][j][1] = DipoleModel::G_by_integration(x1, x2, y1, y2);
    }
    std::cout << i << std::endl;
  }

  std::ofstream out;
  out.open(filepath);
  if (!out.is_open()) {
#ifndef _QUIET
    std::cout << "couldnt open " << filepath << std::endl;
#endif
    exit(EIC_ERROR_COULDNT_OPEN);
  }

  out << std::setprecision(10);

  out << "x ";
  for (uint i = 0; i < num_angles; ++i) {
    out << "angle" << i << "interp angle" << i << "integration ";
  }
  out << std::endl;

  for (unsigned int i = 0; i < num_points; i++) {
    out << b_vec[i] << " ";
    for (uint j = 0; j < num_angles; ++j) {
      out << results[i][j][0] << " " << results[i][j][1] << " ";
    }
    out << std::endl;
  }

  out.close();
}

void hotspot_nucleus_thickness_1d(uint atomic_num,
                                  uint num_hotspots_per_nucleon,
                                  uint num_samples, uint num_points, uint seed,
                                  std::string filepath) {
  if (seed == 0) seed = get_unique_seed();

  HotspotNucleus hn(seed, atomic_num, num_hotspots_per_nucleon);

  std::vector<double> thickness(num_points, 0.0);
  std::vector<double> x(num_points);

  double x_max =
      hn.get_mean_bulk_radius() + 10.0 * hn.get_mean_surface_diffusiveness();
  double x_min = -x_max;
  double inverse_x_divisor = 1.0 / double(num_points - 1);
  for (uint i = 0; i < num_points; i++)
    x[i] = x_min + (x_max - x_min) * double(i) * inverse_x_divisor;

  for (uint i = 0; i < num_samples; i++) {
    for (uint j = 0; j < num_points; j++)
      thickness[j] += hn.get_hotspot_thickness(x[j], 0.0);
    hn.sample();
  }

  double inverse_thickness_divisor = 1.0 / double(num_samples);
  for (uint i = 0; i < num_points; i++)
    thickness[i] *= inverse_thickness_divisor;

  if (filepath == "") filepath = "Data/hotspot_nucleus_thickness_1d.dat";

  std::ofstream out(filepath);
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);

  for (uint i = 0; i < num_points; i++)
    out << x[i] << " " << thickness[i] << std::endl;

  out.close();
}

// Avg over many different configurations of hotspots/nucleons inside nuclei and
// write avg and stddev to file
void hotspot_nucleus_thickness_avg(uint atomic_num,
                                   uint num_hotspots_per_nucleon,
                                   uint start_seed, uint num_events,
                                   std::string filepath) {
  TEST_LOG("In function: hotspot_nucleus_thickness_avg")
  const uint size_x = 1e2;
  const uint size_y = 1e2;

  const double xmin = -40.0;
  const double xmax = 40.0;
  const double ymin = -40.0;
  const double ymax = 40.0;

  const double inverse_num = 1.0 / double(num_events);

  std::vector<std::vector<double>> thickness(size_y,
                                             std::vector<double>(size_x, 0.0));
  std::vector<std::vector<double>> thickness_stddev(
      size_y, std::vector<double>(size_x, 0.0));

  std::vector<double> x(size_x);
  std::vector<double> y(size_y);
  TEST_LOG("Allocating HotspotPos array")
  HotspotPos* pos =
      new HotspotPos[atomic_num * num_hotspots_per_nucleon * num_events];

  for (uint j = 0; j < size_y; ++j)
    y[j] = ymin + (ymax - ymin) * double(j) / double(size_y - 1);

  for (uint k = 0; k < size_x; ++k)
    x[k] = xmin + (xmax - xmin) * double(k) / double(size_x - 1);

  HotspotNucleus hn(start_seed, atomic_num, num_hotspots_per_nucleon,
                    std::sqrt(R_sqr), std::sqrt(rH_sqr));
  for (uint i = 0; i < num_events; ++i) {
    hn.seed(start_seed + i);
    TEST_LOG("Sampling HotspotNucleus with seed " << start_seed + i)
    hn.sample();
    hn.sample_hotspot_weights();
    hn.reset_hotspot_weights();
    TEST_LOG("Finished sampling that HotspotNucleus")

    for (uint n = 0; n < atomic_num; ++n)
      for (uint h = 0; h < num_hotspots_per_nucleon; ++h)
        pos[i * atomic_num * num_hotspots_per_nucleon +
            n * num_hotspots_per_nucleon + h] = *hn.get_hotspot_pos(n, h);

    for (uint j = 0; j < size_y; ++j)
      for (uint k = 0; k < size_x; ++k)
        thickness[j][k] += hn.get_hotspot_thickness(x[k], y[j]) * inverse_num;
  }

  for (uint i = 0; i < num_events; ++i) {
    hn.seed(start_seed + i);
    hn.sample();

    for (uint j = 0; j < size_y; ++j)
      for (uint k = 0; k < size_x; ++k)
        thickness_stddev[j][k] +=
            sqr(hn.get_hotspot_thickness(x[k], y[j]) - thickness[j][k]) *
            inverse_num;
  }

  for (uint j = 0; j < size_y; ++j)
    for (uint k = 0; k < size_x; ++k)
      thickness_stddev[j][k] = std::sqrt(thickness_stddev[j][k]);
  TEST_LOG("Writing to file")
  std::ofstream out(filepath);
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);
  // out << "x y thickness" << std::endl;
  for (uint j = 0; j < size_y; ++j)
    for (uint k = 0; k < size_x; ++k)
      out << x[k] << " " << y[j] << " " << thickness[j][k] << std::endl;
  out.close();

  out.open(filepath + ".stddev");
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);
  for (uint j = 0; j < size_y; ++j)
    for (uint k = 0; k < size_x; ++k)
      out << x[k] << " " << y[j] << " " << thickness_stddev[j][k] << std::endl;
  out.close();

  out.open(filepath + ".pos");
  if (!out.is_open()) exit(EIC_ERROR_COULDNT_OPEN);
  out << "x y" << std::endl;
  for (uint i = 0; i < num_events; ++i)
    for (uint n = 0; n < atomic_num; ++n)
      for (uint h = 0; h < num_hotspots_per_nucleon; ++h)
        out << pos[i * atomic_num * num_hotspots_per_nucleon +
                   n * num_hotspots_per_nucleon + h]
                   .x
            << " "
            << pos[i * atomic_num * num_hotspots_per_nucleon +
                   n * num_hotspots_per_nucleon + h]
                   .y
            << std::endl;

  delete[] pos;
}
}  // namespace Output
