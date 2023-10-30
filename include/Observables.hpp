#pragma once


#include <string>
#include <vector>


namespace Observables
{
    void calculate_dsigma_dt(bool do_coherent, bool do_incoherent, std::string output_file);

    void calculate_dsigma_dt(bool do_coherent, bool do_incoherent, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string output_file = "");

    void calculate_dsigma_dt_nucleus (uint atomic_num, uint num_hotspots, uint seed, std::string filepath = std::string(""));

    void calculate_dsigma_dt_nucleus (uint atomic_num, uint num_hotspots, uint seed, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath = std::string(""));
}