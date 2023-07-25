#include "../include/include.h"
#include <random>
#include <chrono>
#include <thread>
#include <iostream>

double rng01() {
    return drand48();
}


int main (int argc, char** argv) {
    
    //std::cout << Coherent::A_integrand_function(1.20915, 0.199206, 2.82996, 1.30603,0.3,0.5,0.0,1) << std::endl;
    std::cout << SaturationModel::DD(0.299042, -1.3952, 4.59043, -4.67101, 4.9222, 4.78085, 1.82755, 1.69059) << std::endl;
    std::cout << 1.0/Nc/Nc*SaturationModel::DDEigen(0.299042, -1.3952, 4.59043, -4.67101, 4.9222, 4.78085, 1.82755, 1.69059) << std::endl;

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

    std::cout << GBWModel::G(-1,1,1,1) << std::endl;

    return 0;
}