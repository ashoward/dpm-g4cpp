#include "SimMaterialData.hh"

#include <cstdio>
#include <iostream>

void  SimMaterialData::Load(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/mat.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimMaterialData::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 1 lines is comment followed by the secondary e- and gamma production
  // thresholds in [MeV]
  fgets(name, sizeof(name), f);
  fscanf(f, "%lg  %lg\n", &fElectronCut, &fGammaCut);
  // one line comment followed by the 3 MSC max-step lenght paraeters
  fgets(name, sizeof(name), f);
  fscanf(f, "%lg  %lg  %lg\n", &fMSCStepParSLow, &fMSCStepParSHigh, &fMSCStepParEcross);
  // one line comment followed by the number of materials
  fgets(name, sizeof(name), f);
  fscanf(f, "%d\n", &fNumMaterial);
  if (verbose > 0) {
    std::cout << " == Loading Material data for "
              << fNumMaterial << " different materials."
              << std::endl;
  }
  fMaterialName.resize(fNumMaterial);
  fMaterialDensity.resize(fNumMaterial);
  fMollerIMFPScaling.resize(fNumMaterial);
  for (int im=0; im<fNumMaterial; ++im) {
    // 2 lines comments followed by the materail index, name and denisty in [g/cm3]
    fgets(name, sizeof(name), f);
    fgets(name, sizeof(name), f);
    int idum;
    double ddum;
    fscanf(f, "%d %s %lg\n",&idum, name, &ddum);
    fMaterialName[im]   = std::string(name);
    fMaterialDensity[im] = ddum;
    // 3 lines not imprtant followed by 2 numbers: the first is the Moller IMFP scaling factor
    fgets(name, sizeof(name), f);
    fgets(name, sizeof(name), f);
    fgets(name, sizeof(name), f);
    fscanf(f, "%lg %lg\n",&(fMollerIMFPScaling[im]), &ddum);
    if (verbose > 0) {
      std::cout << "    --- Material data for: " << fMaterialName[im] << std::endl;
    }
  }
  fclose(f);
}
