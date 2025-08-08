#include "../include/Output.hpp"
#include "../include/utilities.hpp"

int main(int argc, char** argv) {
  init(argc, argv);

  // Output::G(1000, g_filepath);
  Output::dsigmadt_nucleus(A, H, g_seed);
  // Output::hotspot_nucleus_thickness_avg(A, H, g_seed, 512, g_filepath);
}
