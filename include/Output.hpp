#pragma once
#ifndef EIC_OUTPUT_HPP_
#define EIC_OUTPUT_HPP_

#include <string>
#include <vector>

namespace Output {
void dsigmadt(bool do_coherent, bool do_incoherent);
void dsigmadt(bool do_coherent, bool do_incoherent,
              std::vector<double> Delta_vec, std::vector<double> phi_vec);

void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed);
void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed,
                      std::vector<double> value_vec,
                      std::vector<double> phi_vec);

void dsdt_nucleus_internal_avg(uint atomic_num, uint num_hotspots, uint seed);
void dsdt_nucleus_internal_avg(uint atomic_num, uint num_hotspots, uint seed,
                               std::vector<double> value_vec);

// void dsdt_nucleus_avg_test(uint atomic_num, uint num_hotspots, uint seed,
// uint num_events); void dsdt_nucleus_avg_test(uint atomic_num, uint
// num_hotspots, uint seed, std::vector<double> Delta_vec);

void dsigmadt_demirci(std::string filepath);

void G(uint num_points, std::string filepath);

void hotspot_nucleus_thickness_1d(uint atomic_num,
                                  uint num_hotspots_per_nucleon,
                                  uint num_samples, uint num_points,
                                  uint seed = 0,
                                  std::string filepath = std::string(""));

void hotspot_nucleus_thickness_avg(uint atomic_num,
                                   uint num_hotspots_per_nucleon,
                                   uint start_seed, uint num_events,
                                   std::string filepath);
}  // namespace Output

#endif  // EIC_OUTPUT_HPP_
