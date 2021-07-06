#include "SimMaxScatStrength.hh"

#include <cstdio>
#include <iostream>

void  SimMaxScatStrength::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/el_scatStrength.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimMaxScatStrength::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 4 lines are comments
  for (int i=0; i<4; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the electron energy grid
  int numData;
  fscanf(f, "%d\n", &numData);
  if (verbose > 0) {
    std::cout << " == Loading Maximum Scattering Strength in MSC step data: "
              << numData << " discrete values for Spline interpolation. "
              << std::endl;
  }
  // first 4 lines are comments
  for (int i=0; i<4; ++i) {
    fgets(name, sizeof(name), f);
    if (i==2 && verbose>0) {
      std::cout << "    --- The K_1(E) data were computed for: " << name;
    }
  }
  // load the fNumData E, K_1(E) data and fill in the Spline interplator
  fData.SetSize(numData);
  for (int i=0; i<numData; ++i) {
    double ekin, val;
    fscanf(f, "%lg %lg", &ekin, &val);
    fData.FillData(i, ekin, val);
    if (i==0)         { fEmin = ekin; }
    if (i==numData-1) { fEmax = ekin; }
  }
  fclose(f);
}
