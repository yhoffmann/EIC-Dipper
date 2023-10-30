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


int main (int argc, char** argv)
{
    set_parameters(argc, argv);

    std::string interpolator_filepath = "";

    set_import_filepath_by_m(interpolator_filepath);
    GBWModel::G_ip.import_data(interpolator_filepath);

    std::vector<double> Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0};
    for (uint i=0; i<Delta_vec.size(); i++)
    {
        std::cout << Delta_vec[i] << " " << Coherent::Demirci::dsigma_dt(std::sqrt(0.1), Delta_vec[i]) << std::endl;
    }

    
    
    //std::cout << std::get<0>(Coherent::GeometryAverage::)
    //Observables::calculate_dsigma_dt_nucleus(A, H, 0);
    //GBWModel::G_ip.generate_data(GBWModel::G_wrapper, &config, true);
    //GBWModel::G_ip.export_data(filepath);
/*
    GBWModel::G_ip.import_data(interpolator_filepath);
    GBWModel::G_ip.export_data(interpolator_filepath);
    //std::cout << Coherent::Demirci::dsigma_dt(std::sqrt(0.1), 0.001) << " " << Coherent::dsigma_dt_cubature(std::sqrt(0.1), 0.001) << std::endl;


    std::vector<double> Q_vec = {0.05, std::sqrt(0.1)};
    std::vector<double> Delta_vec;
    
    uint n = 20;
    double min = 0.0;
    double max = 2.5;
    double k = 8.0;
    for (uint i=0; i<n; i++)
    {
        //double Delta = min+(max-min)*(double(i))/(double(n-1));
        double Delta = min+(max-min)*( exp( M_LN2*double(i)/double(n-1)*k )-1.0 )/( std::pow(2.0, k)-1.0 );
        Delta_vec.push_back(Delta);
    }

    std::mt19937 rng;
    HotspotNucleus hn(A, H, rng);
    //std::cout << SaturationModel::GeometryAverage::dsigma_d2b(1.0, 0.1, 0.2, 0.3, &hn) << std::endl;
    //std::cout << std::get<0>(Coherent::GeometryAverage::A(0.31, 0.001, hn)) << std::endl;

*/
    return 0;
}