#include "../include/include.h"

int main (int argc, char** argv) {
    
    //std::cout << Coherent::A_integrand_function(1.20915, 0.199206, 2.82996, 1.30603,0.3,0.5,0.0,1) << std::endl;
    std::cout << 9.0*SaturationModel::DD(0.6674, -2.53162, -7.99407, 0.312111, -1.90802, 7.947, -0.793541, -5.98289) << std::endl;

    CubaConfig c_config;
    IntegrandParams i_params;

    c_config.progress_monitor = true;
    c_config.integrator = 'c';
    i_params.Delta = 0.0;

    std::cout << Incoherent::calculate_dsigma_dt(c_config, i_params)[0] << std::endl;
    
    return 0;
}