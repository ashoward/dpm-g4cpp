#include "SimMollerTables.hh"

#include "SimLinAliasData.hh"

#include <iostream>
#include <cstdio>
#include <cmath>


SimMollerTables::SimMollerTables() {
  // all members will be set when loading the data from the file
  fSamplingTableSize        = -1;
  fNumPrimaryEnergies       = -1;
  fMinPrimaryEnergy         = -1.;
  fLogMinPrimaryEnergy      = -1.;
  fInvLogDeltaPrimaryEnergy = -1.;
}


void SimMollerTables::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/ioni_MollerDtrData.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimMollerTables::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 5 lines are comments
  for (int i=0; i<5; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the primary energy grid and the individual tables first
  fscanf(f, "%d %d", &fNumPrimaryEnergies, &fSamplingTableSize);
  if (verbose >0) {
    std::cout << " == Loading Moller tables: "
              << fNumPrimaryEnergies << " tables with a size of "
              << fSamplingTableSize << " each."
              << std::endl;
  }
  // clean tables if any
  CleanTables();
  fTheTables.resize(fNumPrimaryEnergies, nullptr);
  // load each primary energies and at each primary energy the corresponding table
  for (int ie=0; ie<fNumPrimaryEnergies; ++ie) {
    // load the primary particle kinetic energy value and use if it's needed
    double ddum;
    fscanf(f, "%lg", &ddum);
    if (ie==0) {
      // this is 2x electron-cut
      fMinPrimaryEnergy    = ddum;
      fLogMinPrimaryEnergy = std::log(ddum);
    }
    if (ie==1) {
      fInvLogDeltaPrimaryEnergy = 1./(std::log(ddum)-fLogMinPrimaryEnergy);
    }
    // construct a sampling table, load the data and fill in the sampling table
    fTheTables[ie] = new SimLinAliasData(fSamplingTableSize);
    for (int is=0; is<fSamplingTableSize; ++is) {
      double xdata, ydata, aliasw;
      int    aliasi;
      fscanf(f, "%lg %lg %lg %d", &xdata, &ydata, &aliasw, &aliasi);
      fTheTables[ie]->FillData(is, xdata, ydata, aliasw, aliasi);
    }
  }
  fclose(f);
}

// it is assumed that: 2 x electron-cut < eprim < E_max (also for e+)
double SimMollerTables::SampleEnergyTransfer(double eprim, double rndm1, double rndm2, double rndm3) {
  // determine the primary electron energy lower grid point and sample if that or one above is used now
  double lpenergy   = std::log(eprim);
  double phigher    = (lpenergy-fLogMinPrimaryEnergy)*fInvLogDeltaPrimaryEnergy;
  int penergyindx   = (int) phigher;
  phigher          -= penergyindx;
  if (rndm1<phigher) {
    ++penergyindx;
  }
  // should always be fine if 2 x electron-cut < eprim < E_max (also for e+) but make sure
//  penergyindx       = std::min(fNumPrimaryEnergies-1, penergyindx);
  // sample the transformed variable xi=[kappa-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)]
  // where kappa = ln(eps) with eps = T/T_0
  // so xi= [ln(T/T_0)-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)] that is in [0,1]
  const double   xi = fTheTables[penergyindx]->Sample(rndm2, rndm3);
  // fLogMinPrimaryEnergy is log(2*ecut) = log(ecut) - log(0.5)
  const double dum1 = lpenergy - fLogMinPrimaryEnergy;
  // return with the sampled kinetic energy transfered to the electron
  // (0.5*fMinPrimaryEnergy is the electron production cut)
  return std::exp(xi*dum1)*0.5*fMinPrimaryEnergy;
}


void SimMollerTables::CleanTables() {
  for (std::size_t i=0; i<fTheTables.size(); ++i) {
    if (fTheTables[i]) delete fTheTables[i];
  }
  fTheTables.clear();
}
