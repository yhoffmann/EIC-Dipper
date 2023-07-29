#include "../include/include.h"
#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include "../../../Interpolation3D/include/Interpolator3D.h"
#include "../include/GBWModel.h"

double rng01() {
    return drand48();
}


Interpolator3D global_G_ip;


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/G_integral.txt";

    double counter = 0.0;
    while (false) {
        double bx = (5+counter)*2*(rng01()-0.5);
        double by = (5+counter)*2*(rng01()-0.5);
        double bbx = (5+counter)*2*(rng01()-0.5);
        double bby = (5+counter)*2*(rng01()-0.5);
        double rx = (5+counter)*2*(rng01()-0.5);
        double ry = (5+counter)*2*(rng01()-0.5);
        double rbx = (5+counter)*2*(rng01()-0.5);
        double rby = (5+counter)*2*(rng01()-0.5);

        double printthis = 9.0*SaturationModel::DDEigen(bx,by,bbx,bby,rx,ry,rbx,rby);

        //std::cout << bvalue << "\t" << bbvalue << "\t" << _a << "\t" << _b << "\t" << _c << "\t" << _d << "\t" << GBWModel::dsigma_dip_d2b(bx,by,rx,ry) << std::endl;
        std::cout << printthis << "\t" << bx << ", " << by << ", " << bbx << ", " << bby << "," << rx << ", " << ry << ", " << rbx << ", " << rby << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(150));
        counter += 0.1;
    }

    CubaConfig cuba_config;
    IntegrationConfig integration_config;
    AIntegrandParams A_integrand_params;

    integration_config.integrand_params = &A_integrand_params;

    cuba_config.progress_monitor = true;
    cuba_config.integrator = 'c';
    A_integrand_params.Delta = 0.0;


    global_G_ip.load_data(filepath);

    
    std::cout << GBWModel::G(1.5211,0.178946,1.77934,-0.847651) << " " << GBWModel::G_by_integration(1.5211,0.178946,1.77934,-0.847651) << std::endl;
    std::cout << std::stod("0") << std::endl;

    for (int i=0; i<100000000; i++)
    {
        (void)GBWModel::G_old(1,2,3,4);
    }

    //std::cout << Coherent::dsigma_dt(cuba_config,integration_config)[0] << std::endl;

    return 0;
}