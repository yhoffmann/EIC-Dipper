#pragma once


#include <string>
#include <vector>


namespace Output
{
    void dsigmadt(bool do_coherent, bool do_incoherent, std::string output_file);

    void dsigmadt(bool do_coherent, bool do_incoherent, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec, std::string output_file = "");

    void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed, std::string filepath = std::string(""));

    void dsigmadt_nucleus(uint atomic_num, uint num_hotspots, uint seed, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec, std::string filepath = std::string(""));

    void G(unsigned int num_points, std::string filepath);

    void hotspot_nucleus_thickness_1d(uint atomic_num, uint num_hotspots_per_nucleon, uint num_samples, uint num_points, uint seed = 0, std::string filepath = std::string(""));
}