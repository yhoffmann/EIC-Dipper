#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <unistd.h>
#include "../include/IntegrationRoutines.hpp"
#include "../include/Output.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace Output
{
    void dsigmadt (bool do_coherent, bool do_incoherent, std::string output_file = "")
    {
        std::vector<double> default_Q_vec = {/*0.05, */std::sqrt(0.1)};
        std::vector<double> default_Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0};
        // std::vector<double> DeltaRange = {1.8, 1.9, 2.0, 2.05, 2.1, 2.11, 2.12, 2.13, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.9, 3.0};

        dsigmadt(do_coherent, do_incoherent, default_Q_vec, default_Delta_vec, output_file);
    }
    

    void dsigmadt (bool do_coherent, bool do_incoherent, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath)
    {
        double coherent_results[Q_vec.size()][Delta_vec.size()];
        //double demirci_coherent_results[Q_vec.size()][Delta_vec.size()];
        double incoherent_results[Q_vec.size()][Delta_vec.size()];
#ifndef _QUIET
        std::cout << "Integrating with different parameters..." << std::endl;
#endif
        #pragma omp parallel for schedule(static, 1)
        for (uint j = 0; j < Delta_vec.size(); j++)
        {
            for (uint i = 0; i < Q_vec.size(); i++)
            {
        #ifndef _QUIET
                std::cout << Q_vec[i] << " " << Delta_vec[j] << std::endl;
        #endif
                if (do_coherent)
                    coherent_results[i][j] = Coherent::dsigmadt_test(Q_vec[i], Delta_vec[j]);

                if (do_incoherent)
                    incoherent_results[i][j] = Incoherent::dsigmadt_cubature(Q_vec[i], Delta_vec[j]);
            }
        }

        std::ofstream out_stream;

        if (filepath == "")
            filepath = "Data/dsigma_dt.dat";

        out_stream.open(filepath);
        if (!out_stream.is_open())
            exit(20);

        out_stream << "#Q, Delta, Coher, Incoher; " << Q_vec.size() << " values of Q (for Gnuplot)" << std::endl;
        out_stream << "#1, 2,           3,4" << std::endl;

        for (luint i = 0; i < Q_vec.size(); i++)
        {
            out_stream << "\"Q=" << Q_vec[i] << " Coh.\"" << " " << "\"Q=" << Q_vec[i] << " Incoh.\"" << std::endl;
            for (luint j = 0; j < Delta_vec.size(); j++)
            {
                out_stream << Q_vec[i] << " " << Delta_vec[j] << " " << " " << coherent_results[i][j] << " " << incoherent_results[i][j] << std::endl;
            }
            out_stream << "\n" << std::endl; // Two line breaks for correct gnuplot reading // TODO check if this is necessary
        }

        out_stream.close();
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed, std::string filepath)
    {
        std::vector<double> default_Q_vec = {std::sqrt(0.1)};
        std::vector<double> default_Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1};
        // std::vector<double> DeltaRange = {1.8, 1.9, 2.0, 2.05, 2.1, 2.11, 2.12, 2.13, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.9, 3.0};

        dsigmadt_nucleus(atomic_num, num_hotspots, seed, default_Q_vec, default_Delta_vec, filepath);
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath)
    {
        double coherent_results_real[Q_vec.size()][Delta_vec.size()];
        double coherent_results_imag[Q_vec.size()][Delta_vec.size()];
        double incoherent_results[Q_vec.size()][Delta_vec.size()];

        if (seed==0)
            seed = get_unique_seed();

        std::mt19937 rng(seed);

        HotspotNucleus nucleus(atomic_num, num_hotspots, rng);

        #pragma omp parallel for
        for (uint j = 0; j < Delta_vec.size(); j++)
        {
            for (uint i = 0; i < Q_vec.size(); i++)
            {
                auto [coh_real, coh_imag] = Coherent::Sampled::sqrt_dsigmadt_single_event(Q_vec[i], Delta_vec[j], nucleus);

                coherent_results_real[i][j] = coh_real;
                coherent_results_imag[i][j] = coh_imag;

                incoherent_results[i][j] = Incoherent::Sampled::dsigmadt_single_event(Q_vec[i], Delta_vec[j], nucleus);
            }
        }

        if (filepath == std::string(""))
    #ifndef _DILUTE
            filepath = "Data/raw/dense/" + std::to_string(seed);
    #else
            filepath = "Data/raw/dilute/" + std::to_string(seed);
    #endif
        else filepath += std::to_string(seed);

        std::ofstream out;
        out.open(filepath+"_Amplitude.dat");

        if (!out.is_open())
            exit(20);

        print_infos(out, seed, nucleus);

        out << "##Delta,   Q,        A Co real,Co imag,  A2 Inco\n";
        
        for (uint i = 0; i < Q_vec.size(); i++)
        {
            for (uint j = 0; j < Delta_vec.size(); j++)
            {
                out << std::setprecision(7) << Delta_vec[j] << " " << Q_vec[i] << " " << coherent_results_real[i][j] << " " << coherent_results_imag[i][j] << " " << incoherent_results[i][j] << std::endl;
            }
            out << std::endl;
        }

        out.close();
    }


    void G (unsigned int num_points, std::string filepath)
    {
        double results[num_points][3];
        #pragma omp parallel for schedule(dynamic,1)
        for (unsigned int i=0; i<num_points; i++)
        {
            double x = 100.0*double(i)/double(num_points-1);
            results[i][0] = x;
            results[i][1] = GBWModel::G(x, 0.0, 0.0, 0.0);
            results[i][2] = GBWModel::G_by_integration(x, 0.0, 0.0, 0.0);
            std::cout << i << std::endl;
        }

        std::ofstream out;
        out.open(filepath);
        if (!out.is_open())
        {
#ifndef _QUIET
            std::cout << "couldnt open" << std::endl;
#endif
            exit(20);
        }

        for (unsigned int i=0; i<num_points; i++)
        {
            out << std::setprecision(10) << results[i][0] << " " << results[i][1] << " " << results[i][2] << std::endl;
        }

        out.close();
    }


    void hotspot_nucleus_thickness_1d (uint atomic_num, uint num_hotspots_per_nucleon, uint num_samples, uint num_points, uint seed, std::string filepath)
    {
        if (seed==0)
            seed = get_unique_seed();

        std::mt19937 rng(seed);

        HotspotNucleus hn(atomic_num, num_hotspots_per_nucleon, rng);

        double* thickness = new(std::nothrow) double [num_points];
        double* x = new(std::nothrow) double [num_points];
        if (thickness==nullptr || x==nullptr)
            exit(24);

        double x_max = hn.get_mean_bulk_radius()+10.0*hn.get_mean_surface_diffusiveness(), x_min = -x_max;
        double inverse_x_divisor = 1.0/double(num_points-1);
        for (uint i=0; i<num_points; i++)
        {
            thickness[i] = 0.0;
            x[i] = x_min+(x_max-x_min)*double(i)*inverse_x_divisor;
        }
        
        for (uint i=0; i<num_samples; i++)
        {
            for (uint j=0; j<num_points; j++)
            {
                thickness[j] += hn.get_hotspot_thickness(x[j], 0.0);
            }
            hn.sample_nucleon_pos();
        }

        double inverse_thickness_divisor = 1.0/double(num_samples);
        for (uint i=0; i<num_points; i++)
            thickness[i] *= inverse_thickness_divisor;

        if (filepath=="")
            filepath = "Data/hotspot_nucleus_thickness_1d.dat";

        std::ofstream out(filepath);
        if (!out.is_open())
            exit(20);

        for (uint i=0; i<num_points; i++)
            out << x[i] << " " << thickness[i] << std::endl;

        out.close();

        delete[] thickness;
        delete[] x;
    }
}