#include "../include/Output.hpp"
#include "../include/utilities.hpp"

int main(int argc, char** argv) {
  init(argc, argv);

  Output::dsigmadt_nucleus(1, 3, g_seed);
}
