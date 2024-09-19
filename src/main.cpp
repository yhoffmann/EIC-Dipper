#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <tuple>
#include <thread>
#include <future>
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
    set_parameters(argc, argv);

    // Output::dsigmadt_demirci("Data/dsdt_demirci.dat");
    // GBWModel::G_ip.import_data("interpolator-data/G_rH2_070_m_022.dat");
    Output::dsigmadt_nucleus(1, 3, seed);

    // Output::hotspot_nucleus_thickness_avg(A, H, seed, 1, "data/thickness-event-avg-1-1.dat");
    // Output::hotspot_nucleus_thickness_avg(A, H, seed+1, 1, "data/thickness-event-avg-1-2.dat");
    // Output::hotspot_nucleus_thickness_avg(A, H, seed+2, 1, "data/thickness-event-avg-1-3.dat");
    // Output::hotspot_nucleus_thickness_avg(A, H, seed, 128, "data/thickness-event-avg-128.dat");
    // Output::hotspot_nucleus_thickness_avg(A, H, seed, 256, "data/thickness-event-avg" + filepath_global);
    // Output::hotspot_nucleus_thickness_avg(A, H, seed, 512, "data/thickness-event-avg-512.dat");

    // std::mt19937 rng(12390);
    // HotspotNucleus hn(1, 3, rng);
    // auto [real, imag] = Coherent::Sampled::sqrt_dsigmadt_single_event(std::sqrt(0.1), 2.3, 1.0, hn);
    // std::cout << real << " " << imag << std::endl;
    // double inco = Incoherent::Sampled::dsigmadt_single_event(std::sqrt(0.1), 2.3, 1.0, hn);
    // std::cout << inco << std::endl;

    // Output::G(1000, "Data/G_m_022.dat");

    // Output::dsigmadt(true, true, filepath_global);

    // std::cout << Coherent::dsigmadt(std::sqrt(0.1), 1.5) << std::endl;
    // std::cout << Coherent::dsigmadt_test(std::sqrt(0.1), 1.5) << std::endl;

    //std::cout << Coherent::dsigmadt(std::sqrt(0.1), 0.0001) << "\n" << Incoherent::dsigmadt_cubature(std::sqrt(0.1), 0.0001) << std::endl;

    // std::cout << SaturationModel::dsigma_d2b_sqr_reduced(0.1, 0.2, 0.3, -0.4, 0.5, 0.6, -0.7, -0.8) << " " << SaturationModel::dsigma_d2b_sqr(0.1, 0.2, 0.3, -0.4, 0.5, 0.6, -0.7, -0.8) - SaturationModel::dsigma_d2b(0.1, 0.2, 0.3, -0.4)*SaturationModel::dsigma_d2b(0.5, 0.6, -0.7, -0.8) << std::endl;

    //Output::hotspot_nucleus_thickness_1d(A, H, 1e5, 1e3, 0, filepath_global);
/*
    std::mt19937 rng(1392);
    HotspotNucleus hn(A, H, rng);
    std::cout << Incoherent::Sampled::dsigmadt_single_event(std::sqrt(0.1), std::sqrt(1.0e-8), hn) << std::endl;

    //for (uint i=0, imax=10; i<imax; i++)
        //std::cout << (SaturationModel::dsigma_d2b_sqr(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0, -double(i)/7.0, 0.0, 0.0, double(i)/4.0) - SaturationModel::dsigma_d2b(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0)*SaturationModel::dsigma_d2b(-double(i)/7.0, 0.0, 0.0, double(i)/4.0))- SaturationModel::dsigma_d2b_sqr_reduced(double(i)/10.0, 0.0, 0.0, -double(i)/10.0-4.0, -double(i)/7.0, 0.0, 0.0, double(i)/4.0) << std::endl;
*/
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

    
    std::cout << Q << " " << Delta << "\ncross section, demirci\n" << Coherent::Demirci::dsigmadt(Q, Delta) << "\nmine\n" << Coherent::dsigmadt(Q, Delta) << "\nsampled\n" << sqr( std::get<1>(Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta, hn)) ) << std::endl;
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


    // double Q = std::sqrt(0.1);
    // std::vector<double> Delta_vec;
    // std::vector<double> phi_vec;

    // uint imax = 64;

    // double inco[imax];
    // double color[imax];
    // double hotspot[imax];

    // for (uint i=0; i<imax; ++i)
    // {
    //     Delta_vec.push_back( std::sqrt(16.0)*double(i)/double(imax-1)+0.0001 );
    // }

    // // Output::dsigmadt(true, false, Q, Delta_vec, phi_vec, filepath_global);

    // // #pragma omp parallel for ordered
    // // for (uint i=0; i<imax; ++i)
    // // {
    // //    results[i] = Coherent::Demirci::dsigmadt(Q_vec[0], Delta_vec[i]);
    // //    results2[i] = Coherent::dsigmadt(Q_vec[0], Delta_vec[i]);
    // //    std::cout << i << "\n";
    // // }

    // #pragma omp parallel for ordered
    // for (uint i=0; i<imax; ++i)
    // {
    //     inco[i] = Incoherent::Demirci::dsigmadt(Q, Delta_vec[i]);
    //     // results2[i] = Incoherent::dsigmadt(Q, Delta_vec[i], {0.0});
        
    //     color[i] = Incoherent::Demirci::color_fluctuations(Q, Delta_vec[i]);
    //     hotspot[i] = Incoherent::Demirci::hotspot_fluctuations(Q, Delta_vec[i]);
        
    //     std::cout << i << "\n";
    // }
    
    // if (filepath_global == "")
    //     filepath_global = "Data/inco-dilute-demirci.dat";
    // std::ofstream out(filepath_global);
    
    // for (uint i=0; i<imax; i++)
    //     out << Q << " " << Delta_vec[i] << " " << inco[i] << " " << color[i] << " " << hotspot[i] << std::endl;
    
    // out.close();

    //std::cout << Incoherent::dsigmadt_cubature(std::sqrt(0.1), std::sqrt(8.0)*double(2)/double(15)+0.0001) << std::endl;



    // double Q = std::sqrt(0.1);
    // double Delta = 1.0;

    // const uint num_runs = 0;
    // double results[num_runs];

    // uint seed = get_unique_seed();

    // #pragma omp parallel for ordered
    // for (uint i=0; i<num_runs; ++i)
    // {
    //     std::mt19937 rng(seed+i);
    //     HotspotNucleus hn(1, 3, rng);

    //     results[i] = Incoherent::Sampled::dsigmadt_single_event(Q, Delta, hn);

    //     std::cout << i << "\n";
    // }

    // double sample_avg = 0.0;
    // for (uint i=0; i<num_runs; ++i)
    //     sample_avg += results[i]/double(num_runs);

    // std::future<double> demirci_fut = std::async
    // (
    //     std::launch::async,
    //     [](double Q, double Delta){ return Incoherent::Demirci::color_fluctuations(Q, Delta); },
    //     Q, Delta
    // );
    // std::future<double> analytical_hotspot_avg_fut = std::async
    // (
    //     std::launch::async,
    //     [](double Q, double Delta) { return Incoherent::dsigmadt_cubature(Q, Delta); },
    //     Q, Delta
    // );

    // double demirci = demirci_fut.get();
    // double analytical_hotspot_avg = analytical_hotspot_avg_fut.get();

    // // std::cout << "seed " << seed << std::endl;
    // std::cout <<
    //     "Color Fluctuations"
    //     /*"\nSampling avg " << sample_avg <<*/
    //     "\nAnalytical hotspot avg " << analytical_hotspot_avg <<
    //     "\nDemirci color fluc " << demirci << std::endl;

    return 0;
}