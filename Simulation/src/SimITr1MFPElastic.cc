#include "SimITr1MFPElastic.hh"

#include <cstdio>
#include <iostream>

SimITr1MFPElastic::SimITr1MFPElastic() {
  fNumMaterial = -1;
  fEmin        = -1.;
  fEmax        = -1.;
}


void  SimITr1MFPElastic::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/el_itr1mfp.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimITr1MFPElastic::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 4 lines are comments
  for (int i=0; i<4; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the electron energy grid and #materials
  int numData;
  fscanf(f, "%d  %d\n", &numData, &fNumMaterial);
  if (verbose > 0) {
    std::cout << " == Loading Inverse Tr1-MFP (scalled) data: "
              << numData << " discrete values for Spline interpolation at each of the "
              << fNumMaterial << " different materials."
              << std::endl;
  }
  // skip one line
  fgets(name, sizeof(name), f);
  // allocate space for Spline-interpolation data and load data for each materials
  fDataPerMaterial.resize(fNumMaterial);
  for (int imat=0; imat<fNumMaterial; ++imat) {
    for (int i=0; i<3; ++i) {
      fgets(name, sizeof(name), f);
      if (i==1 && verbose>0) {
        std::cout << "    --- The Inverse Tr1-MFP data were computed for: " << name;
      }
    }
    // load the fNumData E, I-Tr1-MFP/density data and fill in the Spline interplator
    fDataPerMaterial[imat].SetSize(numData);
    for (int i=0; i<numData; ++i) {
      double ekin, val, ddum;
      fscanf(f, "%lg %lg %lg\n", &ekin, &ddum, &val);
      fDataPerMaterial[imat].FillData(i, ekin, val);
      if (imat==0 && i==0)         { fEmin = ekin; }
      if (imat==0 && i==numData-1) { fEmax = ekin; }
    }
  }
  fclose(f);
}
