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

    Output::dsigmadt_nucleus(A, H, 0);

    //for (uint i=0, imax=10; i<imax; i++)
        //std::cout << (SaturationModel::dsigma_d2b_sqr(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0, -double(i)/7.0, 0.0, 0.0, double(i)/4.0) - SaturationModel::dsigma_d2b(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0)*SaturationModel::dsigma_d2b(-double(i)/7.0, 0.0, 0.0, double(i)/4.0))- SaturationModel::dsigma_d2b_sqr_reduced(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0, -double(i)/7.0, 0.0, 0.0, double(i)/4.0) << std::endl;

/*
    uint seed = 580800416;
    std::mt19937 rng(seed);
    HotspotNucleus hn(A, H, rng);

    std::cout << seed << std::endl;
    std::cout << "Hotspot positions" << std::endl;
    for (uint i=0, imax=hn.get_num_hotspots_per_nucleon(); i<imax; i++)
    {
        const double* hp = hn.get_hotspot_pos(0,i);
        // *(double*)(hp) = 0.0;
        // *(double*)(hp+1) = 0.0;
        std::cout << hp[0] << " " << hp[1] << std::endl;
    }
    std::cout << std::endl;
*/

/*
    const uint points = 1e2;
    double xmin = -10.0;
    double xmax = 10.0;
    double dsigmad2b_sampling[points];
    double dsigmad2b[points];

    for (uint n=0; n<points; n++)
    {
        double avg = 0.0;
        double x = xmin+(xmax-xmin)*double(n)/double(points-1);
        for (uint i=0, imax=1e6; i<imax; i++)
        {
            avg += SaturationModel::Sampled::dsigma_d2b(x, 1.0, 0.0, 0.0, &hn)/double(imax);
            hn.sample_hotspot_pos();
        }
        dsigmad2b_sampling[n] = avg;
        dsigmad2b[n] = SaturationModel::dsigma_d2b(x, 1.0, 0.0, 0.0);

        std::cout << x << " " << dsigmad2b[n] << " " << dsigmad2b_sampling[n] << std::endl;
    }

    
    std::cout << Q << " " << Delta << "\ncross section, demirci\n" << Coherent::Demirci::dsigmadt(Q, Delta) << "\nmine\n" << Coherent::dsigmadt_cubature(Q, Delta) << "\nsampled\n" << sqr( std::get<1>(Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta, hn)) ) << std::endl;
*/
/*
    double Q = std::sqrt(0.1);
    double Delta = std::sqrt(1.0);

    uint base_id = 758169;
    std::cout << "id " << base_id << std::endl;
    const uint imax = 300;
    double values[imax];
    double avg = 0.0;
    double avg_of_sqrs = 0.0;

    #pragma omp parallel for
    for (uint i=0; i<imax; i++)
    {
        std::mt19937 rng(base_id+i);
        HotspotNucleus nucleus(1, 3, rng);
        
        values[i] = std::get<1>( Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta, nucleus) );
        avg += values[i]/double(imax);
        avg_of_sqrs += sqr(values[i])/double(imax);

        std::cout << base_id+i << " " << values[i] << std::endl;
    }

    std::cout << "Average of values is " << avg << ". So sqr of avg is " << avg*avg << ".\nAvg of squares is " << avg_of_sqrs << std::endl;
*/
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