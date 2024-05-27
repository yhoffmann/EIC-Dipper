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
    double get_default_Q()
    {
        return std::sqrt(0.1);
    }

    std::vector<double> get_default_Delta_vec()
    {
        return std::vector<double>{0.001, 0.002, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.08, 0.09, 0.12, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0}; // {0.001, 0.01, 0.04, 0.08, 0.12, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.3, 3.0, 4.0}; //
    }

    std::vector<double> get_default_phi_vec()
    {
        std::vector<double> default_phi_vec;
        uint num_angles = 8;
        for (uint i=0; i<num_angles; i++)
            default_phi_vec.push_back(PI*double(i)/double(num_angles)); // no num_angles-1 because we do not want to reach 2pi

        return default_phi_vec;
    }

    void dsigmadt (bool do_coherent, bool do_incoherent, std::string output_file = "")
    {
        double default_Q = get_default_Q();
        std::vector<double> default_Delta_vec = get_default_Delta_vec();
        std::vector<double> default_phi_vec = get_default_phi_vec();

        dsigmadt(do_coherent, do_incoherent, default_Q, default_Delta_vec, default_phi_vec, output_file);
    }
    

    void dsigmadt (bool do_coherent, bool do_incoherent, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec, std::string filepath)
    {
        if (phi_vec.size()==0)
            phi_vec = get_default_phi_vec();

        double coherent_results[Delta_vec.size()][phi_vec.size()];
        double coherent_avg[Delta_vec.size()];

        double incoherent_results[Delta_vec.size()][phi_vec.size()];
        double incoherent_avg[Delta_vec.size()];
#ifndef _QUIET
        std::cout << "Integrating with different parameters..." << std::endl;
#endif
        #pragma omp parallel for schedule(static, 1)
        for (uint i=0; i<Delta_vec.size(); i++)
        {
            coherent_avg[i] = 0.0;
            incoherent_avg[i] = 0.0;

            for (uint j=0, jmax=phi_vec.size(); j<jmax; j++)
            {
                if (progress_monitor_global)
                    std::cout << Q << " " << Delta_vec[i] << " " << phi_vec[j] << std::endl;

                if (do_coherent)
                {
                    coherent_results[i][j] = Coherent::dsigmadt_test(Q, Delta_vec[i], phi_vec[j]);
                    coherent_avg[i] += coherent_results[i][j]/double(jmax);
                }

                if (do_incoherent)
                {
                    incoherent_results[i][j] = Incoherent::dsigmadt_cubature(Q, Delta_vec[i], phi_vec[j]);
                    incoherent_avg[i] += incoherent_results[i][j]/double(jmax);
                }
            }
        }

        std::ofstream out;

        if (filepath == "")
            filepath = "Data/dsigma_dt.dat";

        out.open(filepath);
        if (!out.is_open())
            exit(20);

        out << "#Q, Delta, Coher, Incoher; " << std::endl;
        out << "#1, 2,           3,4" << std::endl;
        out.flush();

        for (uint i=0; i<Delta_vec.size(); i++)
            out << std::setprecision(7) << Delta_vec[i] << " " << Q << " " << coherent_avg[i] << " " << incoherent_avg[i] << std::endl;
        out << std::endl;

        out.close();
        out.open(filepath + ".all");
        if (!out.is_open())
            exit(20);

        for (uint i=0; i<Delta_vec.size(); i++)
        {
            out << std::setprecision(7) << Delta_vec[i] << " " << Q << "   ";
            for (uint j=0; j<phi_vec.size(); j++)
                out << std::setprecision(7) << phi_vec[j] << " " << coherent_results[i][j] << " " << incoherent_results[i][j] << "   ";
            out << std::endl;
        }
        out << std::endl;

        out.close();
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed, std::string filepath)
    {
        double default_Q = get_default_Q();
        std::vector<double> default_Delta_vec = get_default_Delta_vec();
        std::vector<double> default_phi_vec = get_default_phi_vec();

        dsigmadt_nucleus(atomic_num, num_hotspots, seed, default_Q, default_Delta_vec, default_phi_vec, filepath);
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec, std::string filepath)
    {
        double coherent_results_real[Delta_vec.size()][phi_vec.size()];
        double coherent_results_imag[Delta_vec.size()][phi_vec.size()];
        double incoherent_results[Delta_vec.size()][phi_vec.size()];

        if (seed==0)
            seed = get_unique_seed();

        std::mt19937 rng(seed);

        HotspotNucleus nucleus(atomic_num, num_hotspots, rng);

        #pragma omp parallel for
        for (uint i=0; i<Delta_vec.size(); i++)
            for (uint j=0; j<phi_vec.size(); j++)
            {
                auto [coh_real, coh_imag] = Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta_vec[i], phi_vec[j], nucleus);

                coherent_results_real[i][j] = coh_real;
                coherent_results_imag[i][j] = coh_imag;

                incoherent_results[i][j] = Incoherent::Sampled::dsigmadt_single_event(Q, Delta_vec[i], phi_vec[j], nucleus);
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

        for (uint i=0; i<Delta_vec.size(); i++)
        {
            out << std::setprecision(7) << Delta_vec[i] << " " << Q << "   ";
            for (uint j=0; j<phi_vec.size(); j++)
                out << std::setprecision(7) << phi_vec[j] << " " << coherent_results_real[i][j] << " " << coherent_results_imag[i][j] << " " << incoherent_results[i][j] << "   ";
            out << std::endl;
        }
        out << std::endl;

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