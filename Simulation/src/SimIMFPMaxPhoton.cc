#include "SimIMFPMaxPhoton.hh"

#include <cstdio>
#include <iostream>

void  SimIMFPMaxPhoton::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/imfp_globalMax.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimIMFPMaxPhoton::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 3 lines are comments
  for (int i=0; i<3; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the electron energy grid
  int numData;
  fscanf(f, "%d\n", &numData);
  if (verbose >0) {
    std::cout << " == Loading global max of total IMFP data: "
              << numData << " discrete values for Linear interpolation. "
              << std::endl;
  }
  // one additional line of comment
  fgets(name, sizeof(name), f);
  // load the fNumData E, tota-IMFP data and fill in the linear interplator
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
