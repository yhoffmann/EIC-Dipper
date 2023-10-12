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

    std::string filepath = "InterpolatorData/G_BG_m_022.dat";
    DataGenerationConfig config;
    config.n_x = 300;
    config.n_y = 300;
    config.n_z = 30;

    config.x_max = 8.0/m;
    config.y_max = 8.0/m;
    config.z_max = M_PI;

    //set_import_filepath_by_m(filepath, &config);

    //GBWModel::G_ip.generate_data(GBWModel::G_wrapper, &config, true);
    //GBWModel::G_ip.export_data(filepath);

    GBWModel::G_ip.import_data(filepath);

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
    
    std::mt19937 rng(1234);
    Nucleus nucleus(rng, 16);

    for (uint i=0, imax=1e2; i<imax; i++)
    {
        
    }
/*
    for (uint i=0, imax=1e8; i<imax; i++)
    {
        (void)SaturationModel::GeometryAverage::dsigma_d2b(1.0, 0.1, 0.2, 0.3, &nucleus);
    }
*/
    //Observables::calculate_dsigma_dt(true, false, Q_vec, Delta_vec, filepath_global);

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
    return 0;
}