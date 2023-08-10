#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/include.h"
#include "../Interpolation3D/include/Interpolator3D.h" // TODO change include path to the file in this directory


double rng01() {
    return drand48();
}


int main (int argc, char** argv)
{
    double x1 = std::atof(argv[1]);
    double x2 = std::atof(argv[2]);
    double y1 = std::atof(argv[3]);
    double y2 = std::atof(argv[4]);

    std::string filepath = "InterpolatorData/G.txt";

    CubaConfig cuba_config;
    IntegrationConfig integration_config;
    AIntegrandParams A_integrand_params;

    integration_config.integrand_params = &A_integrand_params;

    cuba_config.progress_monitor = true;
    cuba_config.integrator = 'c';
    A_integrand_params.Delta = 0.001;


    //dc.x_grid_spacing = "log";
    //dc.y_grid_spacing = "log";
    //dc.z_grid_spacing = "linear";

    GBWModel::G_ip.import_data(filepath);
/*
    double counter = 0.0;
    while (true) {
        double bx = (5+counter)*2*(rng01()-0.5);
        double by = (5+counter)*2*(rng01()-0.5);
        double bbx = (5+counter)*2*(rng01()-0.5);
        double bby = (5+counter)*2*(rng01()-0.5);
        double rx = (5+counter)*2*(rng01()-0.5);
        double ry = (5+counter)*2*(rng01()-0.5);
        double rbx = (5+counter)*2*(rng01()-0.5);
        double rby = (5+counter)*2*(rng01()-0.5);

        double printthis = SaturationModel::DD(bx,by,bbx,bby,rx,ry,rbx,rby);

        //std::cout << bvalue << "\t" << bbvalue << "\t" << _a << "\t" << _b << "\t" << _c << "\t" << _d << "\t" << GBWModel::dsigma_dip_d2b(bx,by,rx,ry) << std::endl;
        if (printthis>1.0 || printthis<0.0)
        {
            std::cout << "val " << printthis << "\n" << bx << ", " << by << ", " << bbx << ", " << bby << "," << rx << ", " << ry << ", " << rbx << ", " << rby << std::endl;
        }
        
        //std::this_thread::sleep_for(std::chrono::milliseconds(150));
        //counter += 0.1;
    }
*/
    
/*
    const int imax = 1000;
    double results[imax][3];
    #pragma omp parallel for schedule(dynamic,1)
    for (int i=0; i<imax; i++)
    {
        double x = -10.0+20.0*double(i)/double(imax-1);
        results[i][0] = x;
        results[i][1] = GBWModel::G(x,x2,y1,y2);
        results[i][2] = GBWModel::G_by_integration(x,x2,y1,y2);
        std::cout << i << std::endl;
    }

    std::ofstream out;
    out.open("Data/G_interpolation_vs_integration.dat");
    if (!out.is_open()) { std::cout << "couldnt open" << std::endl; exit(0); }

    for (int i=0; i<imax; i++)
    {
        out << std::setprecision(10) << results[i][0] << " " << results[i][1] << " " << results[i][2] << std::endl;
    }

    out.close();
*/

    std::vector<double> Q_vec = {0.05, 0.3};
    std::vector<double> Delta_vec = {0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.13, 0.17, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.8, 1.9, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0};
    //std::vector<double> DeltaRange = {1.8, 1.9, 2.0, 2.05, 2.1, 2.11, 2.12, 2.13, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.9, 3.0};

    double coherent_results[Q_vec.size()][Delta_vec.size()][2];
    double incoherent_results[Q_vec.size()][Delta_vec.size()][2];

    std::cout << "Integrating with different parameters..." << std::endl;

    #pragma omp parallel for schedule(static,1)
    for (int j=0; j<Delta_vec.size(); j++)
    {
        for (int i=0; i<Q_vec.size(); i++)
        {
            //std::cout << "#" << omp_get_thread_num() << std::endl;
            // Prodedure monitor
            std::cout << Q_vec[i] << " " << Delta_vec[j] << std::endl; 

            CubaConfig cuba_config;
            IntegrationConfig integration_config;
            AIntegrandParams A_integrand_params;

            integration_config.integrand_params = &A_integrand_params;

            cuba_config.progress_monitor = true;
            cuba_config.integrator = 'c';
            
            // Saving Q and Delta values in struct
            A_integrand_params.Q = Q_vec[i];
            A_integrand_params.Delta = Delta_vec[j];

            std::vector<double> coherent_result = Coherent::dsigma_dt(cuba_config,integration_config);

            std::vector<double> incoherent_result = Incoherent::dsigma_dt(cuba_config,integration_config);

            coherent_results[i][j][0] = coherent_result[0]; coherent_results[i][j][1] = coherent_result[1];
            incoherent_results[i][j][0] = incoherent_result[0]; incoherent_results[i][j][1] = incoherent_result[1];
        }
    }

    std::string filename = "Data/dsigma_dt.dat";

    std::ofstream out_stream;
    
    out_stream.open(filename);
    // Closing program if it cannot open the file
    if (!out_stream.is_open()) exit(0);

    out_stream << "#Q, Delta, Coher T,L, Incoher T,L; " << Q_vec.size() << " values of Q (for Gnuplot)" << std::endl;
    out_stream << "#1, 2,           3,4          5,6" << std::endl;

    // Printing results to file
    for (luint i=0; i<Q_vec.size(); i++)
    {
        out_stream << "\"Q=" << Q_vec[i] << " Coherent; "<<get_b_range_factor()<<","<<get_r_range_factor()<<"\" \"Q=" << Q_vec[i] << " Incoherent; "<<get_b_range_factor()<<","<<get_r_range_factor() << std::endl;
        for (luint j=0; j<Delta_vec.size(); j++)
        {
            // Writing values to file
            out_stream << Q_vec[i] << " " << Delta_vec[j] << " " << " " << coherent_results[i][j][0] << " " << coherent_results[i][j][1] << " " << incoherent_results[i][j][0] << " " << incoherent_results[i][j][1] << std::endl;
        }
        out_stream << "\n" << std::endl; // Two line breaks for correct gnuplot reading
    }

    // Cleaning up
    out_stream.close();

    return 0;
}