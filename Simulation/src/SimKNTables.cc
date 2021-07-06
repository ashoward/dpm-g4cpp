#include "SimKNTables.hh"

#include "SimLinAliasData.hh"

#include <iostream>
#include <cstdio>
#include <cmath>


SimKNTables::SimKNTables() {
  // all members will be set when loading the data from the file
  fSamplingTableSize        = -1;
  fNumPrimaryEnergies       = -1;
  fMinPrimaryEnergy         = -1.;
  fLogMinPrimaryEnergy      = -1.;
  fInvLogDeltaPrimaryEnergy = -1.;
}


void SimKNTables::LoadData(const std::string& dataDir, int verbose) {
  char name[512];
  sprintf(name, "%s/compton_KNDtrData.dat", dataDir.c_str());
  FILE* f = fopen(name, "r");
  if (!f) {
    std::cerr << " *** ERROR SimKNTables::LoadData: \n"
              << "     file = " << name << " not found! "
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // first 5 lines are comments
  for (int i=0; i<5; ++i) { fgets(name, sizeof(name), f); }
  // load the size of the primary energy grid and the individual tables first
  fscanf(f, "%d %d", &fNumPrimaryEnergies, &fSamplingTableSize);
  if (verbose > 0) {
    std::cout << " == Loading Klein-Nishina tables: "
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
      // this is the gamma-cut that is also the gamma absorption energy
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

// it is assumed that: gamma-cut < egamma < E_max
double SimKNTables::SampleEnergyTransfer(double egamma, double rndm1, double rndm2, double rndm3) {
  const double kEMC2 = 0.510991;
  // determine primary photon energy lower grid point and sample if that or one above is used now
  double lpenergy  = std::log(egamma);
  double phigher   = (lpenergy-fLogMinPrimaryEnergy)*fInvLogDeltaPrimaryEnergy;
  int penergyindx  = (int) phigher;
  phigher         -= penergyindx;
  if (rndm1<phigher) {
    ++penergyindx;
  }
  // should always be fine if gamma-cut < egamma < E_max but make sure
//  penergyindx      = std::min(fNumPrimaryEnergies-1, penergyindx);
  // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
  // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and momentum
  // conservation in case of scattering on free electron at rest).
  // where ep = E_1/E_0 and kappa = E_0/(mc^2)
  double xi = fTheTables[penergyindx]->Sample(rndm2, rndm3);
  // transform it back to eps = E_1/E_0
  // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
  double kappa = egamma/kEMC2;
  return std::exp(std::log(1.+2.*kappa)*(xi-1.)); // eps = E_1/E_0
}


void SimKNTables::CleanTables() {
  for (std::size_t i=0; i<fTheTables.size(); ++i) {
    if (fTheTables[i]) delete fTheTables[i];
  }
  fTheTables.clear();
}
