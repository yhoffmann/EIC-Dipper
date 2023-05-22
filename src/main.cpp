#include"include/include.h"


int main (int argc, char** argv) {
    
    CubaConfig c_config;
    IntegrandParams i_params;
    i_params.func = &Coherent::Trans::AIntegrandWrapper;

    std::cout << cuba_integrate(c_config, i_params) << std::endl;
    
    return 0;
}