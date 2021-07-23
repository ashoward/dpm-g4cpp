//
// Created by mbarbone on 7/22/21.
//

#include "testAPI.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <getopt.h>

static std::string gInputDataDir("../data");       // location of the pre-generated data
static double gPrimaryEnergy = 12.345;  // primary particle energy in [MeV] (> gamma-cut)
static int gNumPrimaries = 1.0E+5;  // generate 100 000 samples from the brem. interaction
static int gMaterialIndex = 0;  // material index (one of those data generated for)
//
static struct option options[] = {
  {"primary-energy      (in [MeV] units)                                - default: 12.345", required_argument, 0, 'e'},
  {"number-of-primaries (number of interactions to sample)              - default: 1.0E+5", required_argument, 0, 'n'},
  {"input-data-dir      (where the pre-generated data are located)      - default: ./data", required_argument, 0, 'd'},
  {"material-index      (one of the material index data generated for)  - default: 0", required_argument, 0, 'm'},
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};
// auxiliary functions for obtaining input arguments
void Help();
void GetOpt(int argc, char *argv[]);

int main(int argc, char *argv[]) {

  GetOpt(argc, argv);

  BremTest bremmTest(gInputDataDir,
                     gPrimaryEnergy,
                     gNumPrimaries,
                     gMaterialIndex);
  bremmTest.simulate();
  bremmTest.writeHists();

  MollerTest mollerTest(gInputDataDir,
                        gPrimaryEnergy,
                        gNumPrimaries,
                        gMaterialIndex);

  mollerTest.simulate();
  mollerTest.writeHists();

  MSCAngularDeflectionTest mscTest(gInputDataDir,
                                   gPrimaryEnergy,
                                   gNumPrimaries,
                                   gMaterialIndex);
  mscTest.simulate();
  mscTest.writeHists();
}

//
// Inplementation of the auxiliary functions for obtaining input ragumets
//
void Help() {
  std::cout << "\n " << std::setw(120) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
  std::cout << "  A bremsstrahlung interaction test."
            << std::endl;
  std::cout << "\n  Usage: test_brem [OPTIONS] \n" << std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout << "\n " << std::setw(120) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
}

void GetOpt(int argc, char *argv[]) {
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "h:e:n:d:m:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
      case 0:c = options[optidx].val;
        /* fall through */
      case 'e':gPrimaryEnergy = std::stod(optarg);
        break;
      case 'n':gNumPrimaries = std::stoi(optarg);
        break;
      case 'd':gInputDataDir = optarg;
        break;
      case 'm':gMaterialIndex = std::stoi(optarg);
        break;
      case 'h':Help();
        exit(-1);
        break;
      default:printf(" *** Unknown input argument: %c\n", c);
        Help();
        exit(-1);
    }
  }
  std::cout << "\n === The test_brem simulation confguration: \n"
            << "\n     - input data directory  : " << gInputDataDir
            << "\n     - primary energy        : " << gPrimaryEnergy << " [MeV]"
            << "\n     - number of interaction : " << gNumPrimaries
            << "\n     - material index        : " << gMaterialIndex
            << std::endl;
}
