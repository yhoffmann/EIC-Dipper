#include "../include/utilities.hpp"
#include "../include/Output.hpp"
#include "../include/Coherent.hpp"
#include "../include/NRPhoton.hpp"
#include "../include/Incoherent.hpp"

int main (int argc, char** argv)
{
    init(argc, argv);

    Output::dsigmadt_nucleus(1, 3, seed);
}