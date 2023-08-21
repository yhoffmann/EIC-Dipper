#include "../include/Observables.h"
#include "../include/IntegrationRoutines.h"
#include "../include/Coherent.h"
#include "../include/Incoherent.h"
#include <iostream>
#include <fstream>
#include "../include/constants.h"


namespace Observables
{
    void calculate_dsigma_dt (bool do_coherent, bool do_incoherent, std::string output_file)
    {
        std::vector<double> default_Q_vec = {0.05, 0.3};
        std::vector<double> default_Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0};
        //std::vector<double> DeltaRange = {1.8, 1.9, 2.0, 2.05, 2.1, 2.11, 2.12, 2.13, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.9, 3.0};

        calculate_dsigma_dt(do_coherent,do_incoherent,default_Q_vec,default_Delta_vec,output_file);
    }


    void calculate_dsigma_dt (bool do_coherent, bool do_incoherent, std::vector<double> Q_vec, std::vector<double> Delta_vec, std::string filepath)
    {
        double coherent_results[Q_vec.size()][Delta_vec.size()][2]{0.0};
        double incoherent_results[Q_vec.size()][Delta_vec.size()][2]{0.0};

        std::cout << "Integrating with different parameters..." << std::endl;

        #pragma omp parallel for schedule(static,1)
        for (uint j=0; j<Delta_vec.size(); j++)
        {
            for (uint i=0; i<Q_vec.size(); i++)
            {
                //std::cout << "#" << omp_get_thread_num() << std::endl;
                // Prodedure monitor
                std::cout << Q_vec[i] << " " << Delta_vec[j] << std::endl; 

                CubaConfig cuba_config;
                IntegrationConfig integration_config;
                AIntegrandParams A_integrand_params;

                integration_config.integrand_params = &A_integrand_params;

                cuba_config.progress_monitor = true;
                
                // Saving Q and Delta values in struct
                A_integrand_params.Q = Q_vec[i];
                A_integrand_params.Delta = Delta_vec[j];

                if (do_coherent)
                {
                    std::vector<double> coherent_result = Coherent::dsigma_dt(&cuba_config,&integration_config);
                    coherent_results[i][j][0] = coherent_result[0]; coherent_results[i][j][1] = coherent_result[1];
                }

                if (do_incoherent)
                {
                    std::vector<double> incoherent_result = Incoherent::dsigma_dt(&cuba_config,&integration_config);
                    incoherent_results[i][j][0] = incoherent_result[0]; incoherent_results[i][j][1] = incoherent_result[1];
                }
            }
        }

        std::ofstream out_stream;
        
        out_stream.open(filepath);
        // Closing program if it cannot open the file
        if (!out_stream.is_open()) exit(0);

        out_stream << "#Q, Delta, Coher T,L, Incoher T,L; " << Q_vec.size() << " values of Q (for Gnuplot)" << std::endl;
        out_stream << "#1, 2,           3,4          5,6" << std::endl;

        // Printing results to file
        for (luint i=0; i<Q_vec.size(); i++)
        {
            out_stream << "\"Q=" << Q_vec[i] << " Coherent; "<<B_RANGE_FACTOR<<","<<R_RANGE_FACTOR<<"\" \"Q=" << Q_vec[i] << " Incoherent; "<<B_RANGE_FACTOR<<","<<R_RANGE_FACTOR << std::endl;
            for (luint j=0; j<Delta_vec.size(); j++)
            {
                // Writing values to file
                out_stream << Q_vec[i] << " " << Delta_vec[j] << " " << " " << coherent_results[i][j][0] << " " << coherent_results[i][j][1] << " " << incoherent_results[i][j][0] << " " << incoherent_results[i][j][1] << std::endl;
            }
            out_stream << "\n" << std::endl; // Two line breaks for correct gnuplot reading
        }

        // Cleaning up
        out_stream.close();
    }
}