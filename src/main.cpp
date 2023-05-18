//#include"include/include.h"
#include"include/integration_routines.cpp"

int main (int argc, char** argv) {
    IntegrationConfig ic;

    cuba_bessel_integrate(ic);
    
    return 0;
}