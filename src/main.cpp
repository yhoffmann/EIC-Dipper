#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <tuple>
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Interpolation3D/external/easy-progress-monitor/include/ProgressMonitor.hpp"
#include "../include/Output.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/NRPhoton.hpp"


int main (int argc, char** argv)
{
#ifndef _QUIET
#ifndef _DILUTE
    std::cout << "Dense" << std::endl;
#else
    std::cout << "Dilute" << std::endl;
#endif
#endif
    set_parameters(argc, argv);
    
    std::string interpolator_filepath = "";

    set_import_filepath_by_parameters(interpolator_filepath);

    GBWModel::G_ip.import_data(interpolator_filepath);

    //Output::dsigmadt_nucleus(A, H, 0);

    std::mt19937 rng(getpid());
    HotspotNucleus hn(A, H, rng);
    double avg = 0.0;
    for (uint i=0, imax=1e6; i<imax; i++)
    {
        avg += SaturationModel::Sampled::dsigma_d2b(1.0, 2.0, 3.0, -0.5, &hn)/double(imax);
        hn.sample_hotspot_pos();
    }
    std::cout << "hotspot averaged " << avg << std::endl;
    std::cout << SaturationModel::dsigma_d2b(0.5, 1.0, 2.0, -0.5) << std::endl;
    double Q = std::sqrt(0.1);
    double Delta = 0.01;
    //std::cout << Q << " " << Delta << "\n" << std::get<1>(Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta, hn)) << std::endl;

    //std::cout << Coherent::Demirci::dsigma_dt(std::sqrt(0.1), 0.001) << std::endl;
    //std::cout << Coherent::dsigma_dt_cubature(std::sqrt(0.1), 0.001) << std::endl;

    //std::vector<double> Q_vec = {std::sqrt(0.1)};
    //std::vector<double> Delta_vec;
    //for (uint i=0, imax=20; i<imax; ++i)
    //{
    //    Delta_vec.push_back( 1.0*double(i)/double(imax-1)+1.0001 );
    //}

    //std::cout << Incoherent::dsigmadt_cubature(std::sqrt(0.1), 0.001) << std::endl;

    //std::cout << SaturationModel::dsigma_d2b_sqr(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) << std::endl;

    //std::ofstream out(filepath_global);
    //for (uint i=0, imax=Delta_vec.size(); i<imax; ++i)
    //    out << Delta_vec[i] << " " << Coherent::Demirci::dsigmadt(Q_vec[0], Delta_vec[i]) << std::endl;
    //out.close();

    //Output::dsigmadt(true, true, filepath_global);

    //Output::hotspot_nucleus_thickness_1d(A, H, 1e5, 1e3, 0, filepath_global);

    return 0;
}