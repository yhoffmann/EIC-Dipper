#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <unistd.h>
#include "../include/IntegrationRoutines.hpp"
#include "../include/Observables.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace Observables
{
    void calculate_dsigma_dt (bool do_coherent, bool do_incoherent, std::string output_file = "")
    {
        std::vector<double> default_Q_vec = {0.05, 0.3};
        std::vector<double> default_Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0};
        // std::vector<double> DeltaRange = {1.8, 1.9, 2.0, 2.05, 2.1, 2.11, 2.12, 2.13, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.9, 3.0};

        calculate_dsigma_dt(do_coherent, do_incoherent, default_Q_vec, default_Delta_vec, output_file);
    }
    

    void calculate_dsigma_dt (bool do_coherent, bool do_incoherent, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath)
    {
        double coherent_results[Q_vec.size()][Delta_vec.size()];
        double incoherent_results[Q_vec.size()][Delta_vec.size()];

        std::cout << "Integrating with different parameters..." << std::endl;

        #pragma omp parallel for schedule(static, 1)
        for (uint j = 0; j < Delta_vec.size(); j++)
        {
            for (uint i = 0; i < Q_vec.size(); i++)
            {
                std::cout << Q_vec[i] << " " << Delta_vec[j] << std::endl;

                CubatureConfig c_config;
                IntegrationConfig integration_config;
                AIntegrandParams A_integrand_params;

                integration_config.integrand_params = &A_integrand_params;

                c_config.progress_monitor = true;

                A_integrand_params.Q = Q_vec[i];
                A_integrand_params.Delta = Delta_vec[j];

                if (do_coherent)
                    coherent_results[i][j] = Coherent::dsigma_dt_cubature(&c_config, &integration_config);

                if (do_incoherent)
                    incoherent_results[i][j] = Incoherent::dsigma_dt_cubature(&c_config, &integration_config);
            }
        }

        std::ofstream out_stream;

        out_stream.open(filepath);
        if (!out_stream.is_open())
            exit(0);

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


    void calculate_dsigma_dt_nucleus (uint seed, uint atomic_num, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath)
    {
        double coherent_results_real[Q_vec.size()][Delta_vec.size()];
        double coherent_results_imag[Q_vec.size()][Delta_vec.size()];
        double incoherent_results[Q_vec.size()][Delta_vec.size()];

        if (seed != 0)
        {
            std::thread::id thread_id_temp = std::this_thread::get_id();
            std::hash<std::thread::id> hash;
            uint thread_id = hash(thread_id_temp);

            seed = thread_id+getpid();
        }

        std::mt19937 rng(seed);

        HotspotNucleus nucleus(atomic_num, 1, rng);

        #pragma omp parallel for
        for (uint j = 0; j < Delta_vec.size(); j++)
        {
            for (uint i = 0; i < Q_vec.size(); i++)
            {
                double Q = Q_vec[i];
                double Delta = Delta_vec[j];

                auto [coh_real, coh_imag] = Coherent::GeometryAverage::A(Q, Delta, nucleus);

                coherent_results_real[i][j] = coh_real;
                coherent_results_imag[i][j] = coh_imag;

                incoherent_results[i][j] = Incoherent::GeometryAverage::A(Q, Delta, nucleus);
            }
        }

        if (filepath == std::string(""))
            filepath = "Data/" + std::to_string(seed);

        std::ofstream out;
        out.open(filepath+"_Amplitude");
        if (!out.is_open())
            exit(648002);

        out << "#Delta,   Q,        A Co real,Co imag,  A2 Inco\n";
        
        for (uint j = 0; j < Delta_vec.size(); j++)
        {
            for (uint i = 0; i < Q_vec.size(); i++)
            {
                out << std::setprecision(7) << Delta_vec[j] << " " << Q_vec[i] << " " << coherent_results_real[i][j] << " " << coherent_results_imag[i][j] << " " << incoherent_results[i][j] << std::endl;
            }
            out << std::endl;
        }

        out.close();

        create_info_file(filepath+"_Info");
    }


    void calculate_G (unsigned int num_points, std::string filepath)
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
            std::cout << "couldnt open" << std::endl;
            exit(0);
        }

        for (unsigned int i=0; i<num_points; i++)
        {
            out << std::setprecision(10) << results[i][0] << " " << results[i][1] << " " << results[i][2] << std::endl;
        }

        out.close();
    }
}