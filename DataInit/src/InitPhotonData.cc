
#include "InitPhotonData.hh"


#include "PhotonData.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4ProductionCutsTable.hh"
#include "G4DataVector.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"

#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

void InitPhotonData(PhotonData& phData, double electronCutInEnergy, double gammaCutInEnergy) {
  phData.fElectronCut = electronCutInEnergy;
  phData.fGammaCut    = gammaCutInEnergy;
  InitMXsecData(phData);
  phData.fKNTables->BuildTables(phData.fGammaCut, phData.fPhEGrid[phData.fNumPhEData-1]);
}

void InitMXsecData(PhotonData& phData) {
  // NOTE: we will use the Livermore models for photons to get the cross sections
  //       that are created and initialised here
  G4LivermoreComptonModel         theCompton;
  G4LivermoreGammaConversionModel thePair;
  G4LivermorePhotoElectricModel   thePhoelec;

  theCompton.SetHighEnergyLimit(100.0*CLHEP::MeV);
  thePair.SetHighEnergyLimit(100.0*CLHEP::MeV);
  thePhoelec.SetHighEnergyLimit(100.0*CLHEP::MeV);

  G4ParticleDefinition* part = G4Gamma::Gamma();

  // get cuts for secondary e-
  const G4DataVector* theElCuts = static_cast<const G4DataVector*>(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1));
  theCompton.Initialise(part, *theElCuts);
  thePair.Initialise(part, *theElCuts);
  thePhoelec.Initialise(part, *theElCuts);
  // get cuts for secondary gamma
//  const G4DataVector* theGamCuts = static_cast<const G4DataVector*>(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(0));
//  theBrem.Initialise(part, *theGamCuts);
  //
  // compute the macroscopic cross sections of the 3 interactions for the individual
  // materials and also find the overall maximum value of their sum that will be used
  // to play the Woodcoc trick (delta scattering)
  int numMat = phData.fNumMaterial;
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
    G4Material* mat = G4NistManager::Instance()->GetMaterial(im+1);
    phData.fDataPerMaterial[im]->fMaterialName    = mat->GetName();
    phData.fDataPerMaterial[im]->fMaterialDensity = mat->GetDensity()/(CLHEP::g/CLHEP::cm3);
    // compute the macroscopic cross sections over the photon energy grid
    int numPhE = phData.fNumPhEData;
    for (int ie=0; ie<numPhE; ++ie) {
      double phE = phData.fPhEGrid[ie];
      double mxCompton = std::max(1.0E-20, theCompton.CrossSectionPerVolume(mat, part, phE));
      double mxPair    = std::max(1.0E-20, thePair.CrossSectionPerVolume(mat, part, phE));
      double mxPhoElec = std::max(1.0E-20, thePhoelec.CrossSectionPerVolume(mat, part, phE));
      double mxTotal   = mxCompton + mxPair + mxPhoElec;
      if (im==0 || phData.fGlobalMaxTotalMXsec[ie] < mxTotal) {
        phData.fGlobalMaxTotalMXsec[ie] = mxTotal;
      };
      phData.fDataPerMaterial[im]->fComptonMXsec[ie] = mxCompton;
      phData.fDataPerMaterial[im]->fPairMXsec[ie]    = mxPair;
      phData.fDataPerMaterial[im]->fPhoElecMXsec[ie] = mxPhoElec;
      phData.fDataPerMaterial[im]->fTotalMXsec[ie]   = mxTotal;
    }
  }
  //
  //
}
