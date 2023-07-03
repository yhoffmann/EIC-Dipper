#include "../include/include.h"
#include <random>
#include <chrono>
#include <thread>

double rng01() {
    return drand48();
}


int main (int argc, char** argv) {
    
    //std::cout << Coherent::A_integrand_function(1.20915, 0.199206, 2.82996, 1.30603,0.3,0.5,0.0,1) << std::endl;
    std::cout << 9.0*SaturationModel::DD(06.674, -25.3162, -79.9407, 03.12111, -19.0802, 79.47, -07.93541, -59.8289) << std::endl;
    std::cout << SaturationModel::DDEigen(06.674, -25.3162, -79.9407, 03.12111, -19.0802, 79.47, -07.93541, -59.8289) << std::endl;

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

    CubaConfig c_config;
    IntegrandParams i_params;

    c_config.progress_monitor = true;
    c_config.integrator = 'c';
    i_params.Delta = 0.0;

    std::cout << Incoherent::calculate_dsigma_dt(c_config, i_params)[0] << std::endl;

    return 0;
}