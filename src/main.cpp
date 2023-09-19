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


int main (int argc, char** argv)
{
    set_parameters(argc, argv);

    std::string filepath = "InterpolatorData/G";
    DataGenerationConfig config;
    config.n_x = 300;
    config.n_y = 300;
    config.n_z = 30;

    config.x_max = 8.0/m;
    config.y_max = 8.0/m;
    config.z_max = M_PI;

    set_import_filepath_by_m(filepath, &config);

    //GBWModel::G_ip.generate_data(GBWModel::G_wrapper, &config, true);
    //GBWModel::G_ip.export_data(filepath);

    GBWModel::G_ip.import_data(filepath);

    //auto [t1,l1] = Coherent::dsigma_dt(0.3,0.4);
    //auto [t2,l2] = Coherent::dsigma_dt_cubature(0.3,0.4);
    //std::cout << t1 << " " << l1 << "\t" << t2 << " " << l2 << "\n";

    std::vector<double> Q_vec = {0.05, 0.3};
    std::vector<double> Delta_vec;
    
    uint imax = 20;
    for (uint i=0; i<imax; i++)
    {
        Delta_vec.push_back(3.0*double(i)/double((imax-1))+0.001);
    }

    Observables::calculate_dsigma_dt(true,false,Q_vec,Delta_vec,"Data/dsigma_dt_m_022_.dat");

/*
    double x1 = std::atof(argv[1]);
    double x2 = std::atof(argv[2]);
    double y1 = std::atof(argv[3]);
    double y2 = std::atof(argv[4]);
    (void)GBWModel::G_by_integration(x1,x2,y1,y2);
*/
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
        double x = 100.0*double(i)/double(imax-1);
        results[i][0] = x;
        results[i][1] = GBWModel::G(x,0.0,0.0,0.0);
        results[i][2] = GBWModel::G_by_integration(x,0.0,0.0,0.0);
        std::cout << i << std::endl;
    }

    std::ofstream out;
    out.open("Data/G_m_022.dat");
    if (!out.is_open()) { std::cout << "couldnt open" << std::endl; exit(0); }

    for (int i=0; i<imax; i++)
    {
        out << std::setprecision(10) << results[i][0] << " " << results[i][1] << " " << results[i][2] << std::endl;
    }

    out.close();
*/
    return 0;
}