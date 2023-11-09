#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Interpolation3D/external/easy-progress-monitor/include/ProgressMonitor.hpp"
#include "../include/Observables.hpp"
#include "../include/Coherent.hpp"
#include "../include/SaturationModel.hpp"

#include "../include/NRPhoton.hpp"
int main (int argc, char** argv)
{
#ifndef _DILUTE
    std::cout << "Dense" << std::endl;
#else
    std::cout << "Dilute" << std::endl;
#endif
    set_parameters(argc, argv);
    
    std::string interpolator_filepath = "";

    set_import_filepath_by_parameters(interpolator_filepath);

    GBWModel::G_ip.import_data(interpolator_filepath);

    //std::cout << Coherent::Demirci::dsigma_dt(std::sqrt(0.1), 0.001) << std::endl;
    //std::cout << Coherent::dsigma_dt_cubature(std::sqrt(0.1), 0.001) << std::endl;

    std::vector<double> Q_vec = {std::sqrt(0.1)};
    std::vector<double> Delta_vec;
    for (uint i=0, imax=10; i<imax; ++i)
    {
        Delta_vec.push_back( 1.5*double(i)/double(imax-1)+0.01 );
    }

    std::ofstream out(filepath_global);
    for (uint i=0, imax=Delta_vec.size(); i<imax; ++i)
        out << Delta_vec[i] << " " << Coherent::Demirci::dsigmadt(Q_vec[0], Delta_vec[i]) << std::endl;
    out.close();

    //Observables::calculate_dsigma_dt(true, false, Q_vec, Delta_vec, filepath_global);

    return 0;
}