#include "../include/Output.hpp"
#include "../include/utilities.hpp"

int main(int argc, char** argv) {
  init(argc, argv);

  // Output::dsigmadt_nucleus(1, 3, g_seed);
  Output::hotspot_nucleus_thickness_avg(A, H, g_seed, 16, g_filepath);
}
