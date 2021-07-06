#include "SimPhotonData.hh"

#include "SimIMFPMaxPhoton.hh"
#include "SimIMFPPhoton.hh"

#include "SimKNTables.hh"

#include <iostream>

SimPhotonData::SimPhotonData() {
  fIMFPTotalMax      = nullptr;
  fIMFPTotal         = nullptr;
  fIMFPCompton       = nullptr;
  fIMFPPairProd      = nullptr;

  fTheKNTables       = nullptr;
}

SimPhotonData::~SimPhotonData() {
  if (fIMFPTotalMax)     delete fIMFPTotalMax;
  if (fIMFPTotal)        delete fIMFPTotal;
  if (fIMFPCompton)      delete fIMFPCompton;
  if (fIMFPPairProd)     delete fIMFPPairProd;

  if (fTheKNTables)      delete fTheKNTables;
}

void SimPhotonData::Load(const std::string& dataDir, int verbose) {
  //
  // construct and load the global max of the total IMFP data
  if (fIMFPTotalMax) delete fIMFPTotalMax;
  fIMFPTotalMax = new SimIMFPMaxPhoton();
  fIMFPTotalMax->LoadData(dataDir, verbose);
  // cross check
  if (verbose > 1) {
    std::cout << " Global Max of Total-IMFP(E=20.29795305) = "
              << fIMFPTotalMax->GetIMFP(1.5181)
              << std::endl;
  }
  //
  // construct and load the global max of the total IMFP data
  if (fIMFPTotal) delete fIMFPTotal;
  fIMFPTotal = new SimIMFPPhoton(0);
  fIMFPTotal->LoadData(dataDir, verbose);
  // cross check
  if (verbose > 1) {
    std::cout << "  Total-IMFP(E=20.29795305, imat=0) = "
              << fIMFPTotal->GetIMFPPerDensity(1.5181, 0)
              << std::endl;
  }
  //
  // construct and load the Compton IMFP data
  if (fIMFPCompton) delete fIMFPCompton;
  fIMFPCompton = new SimIMFPPhoton(1);
  fIMFPCompton->LoadData(dataDir, verbose);
  // cross check
  if (verbose > 1) {
    std::cout << " Compton-IMFP(E=20.29795305, imat=0) = "
              << fIMFPCompton->GetIMFPPerDensity(1.5181, 0)
              << std::endl;
  }
  // construct and load the Pair-production IMFP data
  if (fIMFPPairProd) delete fIMFPPairProd;
  fIMFPPairProd = new SimIMFPPhoton(2);
  fIMFPPairProd->LoadData(dataDir, verbose);
  // cross check
  if (verbose > 1) {
    std::cout << " PairProd-IMFP(E=20.29795305, imat=0) = "
              << fIMFPPairProd->GetIMFPPerDensity(1.5181, 0)
              << std::endl;
  }
  //
  // construct and load the Klein-Nishina samping tables
  if (fTheKNTables) delete fTheKNTables;
  fTheKNTables = new SimKNTables();
  fTheKNTables->LoadData(dataDir, verbose);
  //
}
