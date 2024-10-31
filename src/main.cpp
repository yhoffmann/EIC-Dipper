#include "../include/utilities.hpp"
#include "../include/Output.hpp"


int main (int argc, char** argv)
{
    set_parameters(argc, argv);

    // Output::dsdt_nucleus_avg_test(Q, std::sqrt(0.1));
    Output::dsigmadt_nucleus(1, 3, seed);
    // std::cout << NRPhoton::wave_function_factor_T << " " << NRPhoton::wave_function_factor_L(Q) << " " << NRPhoton::wave_function_factor(Q) << std::endl;
}