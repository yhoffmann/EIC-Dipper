#pragma once


#include <string>
#include <vector>


namespace Observables
{
    void calculate_dsigma_dt(bool do_coherent, bool do_incoherent, std::string output_file);

    void calculate_dsigma_dt(bool do_coherent, bool do_incoherent, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string output_file = "");

    void calculate_dsigma_dt_nucleus(uint event_id, uint atomic_num, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath = "");
}