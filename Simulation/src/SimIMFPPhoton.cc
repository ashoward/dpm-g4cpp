#include "SimIMFPPhoton.hh"

#include <cstdio>
#include <iostream>


SimIMFPPhoton::SimIMFPPhoton (int type) {
  // all will be set at LoadData()
  fNumMaterial = -1;
  fType        = type;
  fEmin        = -1.;
  fEmax        = -1.;
}

SimIMFPPhoton::~SimIMFPPhoton() {
  fDataPerMaterial.clear();
}

void  SimIMFPPhoton::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  const std::string strType = (fType == 0) ? "total" : ( fType == 1 ? "compton" : "pairp" );
  sprintf(name, "%s/imfp_%s.dat", dataDir.c_str(), strType.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimIMFPPhoton::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 3 lines are comments
  for (int i=0; i<3; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the electron energy grid and #materials
  int numData;
  fscanf(f, "%d  %d\n", &numData, &fNumMaterial);
  if (verbose > 0) {
    std::cout << " == Loading '" << strType << "' IMFP (scalled) data per-material: "
              << numData << " discrete values for Linear interpolation at each of the "
              << fNumMaterial << " different materials."
              << std::endl;
  }
  // skip one line
  fgets(name, sizeof(name), f);
  // allocate space for Linear-interpolation data and load data for each materials
  fDataPerMaterial.resize(fNumMaterial);
  for (int imat=0; imat<fNumMaterial; ++imat) {
    for (int i=0; i<3; ++i) {
      fgets(name, sizeof(name), f);
      if (i==1 && verbose>0) {
        std::cout << "    --- The IMFP data were computed for: " << name;
      }
    }
    // load the fNumData E, tota-IMFP/density data and fill in the Linear interplator
    fDataPerMaterial[imat].SetSize(numData);
    for (int i=0; i<numData; ++i) {
      double ekin, val;
      fscanf(f, "%lg %lg\n", &ekin, &val);
      fDataPerMaterial[imat].FillData(i, ekin, val);
      if (imat==0 && i==0)         { fEmin = ekin; }
      if (imat==0 && i==numData-1) { fEmax = ekin; }
    }
  }
  fclose(f);
}
