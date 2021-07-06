#include "SimElectronData.hh"

#include <iostream>

#include "SimITr1MFPElastic.hh"
#include "SimMaxScatStrength.hh"

#include "SimIMFPMoller.hh"
#include "SimIMFPBrem.hh"
#include "SimStoppingPower.hh"

#include "SimMollerTables.hh"
#include "SimSBTables.hh"
#include "SimGSTables.hh"

SimElectronData::SimElectronData() {
  fITr1MFPElastic    = nullptr;
  fMaxScatStrength   = nullptr;

  fIMFPMoller        = nullptr;
  fIMFPBrem          = nullptr;

  fDEDX              = nullptr;

  fTheMollerTables   = nullptr;
  fTheSBTables       = nullptr;
  fTheGSTables       = nullptr;
}

SimElectronData::~SimElectronData() {
  if (fITr1MFPElastic)   delete fITr1MFPElastic;
  if (fMaxScatStrength)  delete fMaxScatStrength;

  if (fIMFPMoller)       delete fIMFPMoller;
  if (fIMFPBrem)         delete fIMFPBrem;

  if (fDEDX)             delete fDEDX;

  if (fTheMollerTables)  delete fTheMollerTables;
  if (fTheSBTables)      delete fTheSBTables;
  if (fTheGSTables)      delete fTheGSTables;
}

void SimElectronData::Load(const std::string& dataDir, int verbose) {
  //
  // construct and load the elastic ITr1MFP data
  if (fITr1MFPElastic) delete fITr1MFPElastic;
  fITr1MFPElastic = new SimITr1MFPElastic();
  fITr1MFPElastic->LoadData(dataDir, verbose);
  // cross-check
  if (verbose>1) {
    std::cout << " ITr1-MFP(E=20.2427751, imat=0) = "
              << fITr1MFPElastic->GetITr1MFPPerDensity(20.2427751, 0)
              << std::endl;
  }
  //
  // construct and load the elastic ITr1MFP data
  if (fMaxScatStrength) delete fMaxScatStrength;
  fMaxScatStrength = new SimMaxScatStrength();
  fMaxScatStrength->LoadData(dataDir, verbose);
  // cross-check
  if (verbose>1) {
    std::cout << " K1(E=20.2427751) = "
              << fMaxScatStrength->GetMaxScatStrength(20.2427751)
              << std::endl;
  }
  //
  // construct and load the Moller IMFP data
  if (fIMFPMoller) delete fIMFPMoller;
  fIMFPMoller = new SimIMFPMoller();
  fIMFPMoller->LoadData(dataDir, verbose);
  // cross-check
  if (verbose>1) {
    std::cout << " IMFP-Moller(E=20.35517) = "
              << fIMFPMoller->GetIMFPPerDensity(20.35517)
              << std::endl;
  }
  //
  // construct and load the Brem. (Seltzer-Berger) IMFP data
  if (fIMFPBrem) delete fIMFPBrem;
  fIMFPBrem = new SimIMFPBrem();
  fIMFPBrem->LoadData(dataDir, verbose);
  // cross-check
  if (verbose>1) {
    std::cout << " IMFP-Brem(E=20.244377, imat=0) = "
              << fIMFPBrem->GetIMFPPerDensity(20.244377, 0)
              << std::endl;
  }
  //
  // construct and load the restricted stopping power data
  if (fDEDX) delete fDEDX;
  fDEDX = new SimStoppingPower();
  fDEDX->LoadData(dataDir, verbose);
  // cross-check
  if (verbose>1) {
    std::cout << " dE/dx(E=20.244377, imat=0) = "
              << fDEDX->GetDEDXPerDensity(20.244377, 0)
              << std::endl;
  }
  //
  // construct and load the Moller samping tables
  if (fTheMollerTables) delete fTheMollerTables;
  fTheMollerTables = new SimMollerTables();
  fTheMollerTables->LoadData(dataDir, verbose);

  // construct and load the Seltzer-Berger samping tables
  if (fTheSBTables) delete fTheSBTables;
  fTheSBTables = new SimSBTables();
  fTheSBTables->LoadData(dataDir, verbose);

  // construct and load the GS samping tables
  if (fTheGSTables) delete fTheGSTables;
  fTheGSTables = new SimGSTables();
  fTheGSTables->LoadData(dataDir, verbose);
}
