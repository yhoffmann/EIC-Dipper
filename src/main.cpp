#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>
#include <tuple>
#include <thread>
#include <future>
#include "../include/constants.hpp"
#include "../include/utilities.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Interpolation3D/include/Interpolator3D.hpp"
#include "../external/Interpolation3D/external/easy-progress-monitor/include/ProgressMonitor.hpp"
#include "../include/Output.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"
#include "../include/SaturationModel.hpp"
#include "../include/NRPhoton.hpp"


int main (int argc, char** argv)
{
    set_parameters(argc, argv);

    // Output::dsdt_nucleus_avg_test(Q, std::sqrt(0.1));
    Output::dsigmadt_nucleus(1, 3, seed);
    // std::cout << NRPhoton::wave_function_factor_T << " " << NRPhoton::wave_function_factor_L(Q) << " " << NRPhoton::wave_function_factor(Q) << std::endl;

    return 0;
}