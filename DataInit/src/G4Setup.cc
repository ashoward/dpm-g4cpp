
#include "G4Setup.hh"

// Geant4 includes
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4String.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"

#include "G4DataVector.hh"
#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"

#include "Spline.hh"

// builds a fake Geant4 geometry with all materials given just to be able to produce material-cuts couple
void FakeG4Setup (const std::vector<std::string>& g4Materials, double electronCutInEnergy, double gammaCutInEnergy, int verbose) {
  //
  // --- Geometry definition: create the word
  G4double wDimX      = 0.6*mm;
  G4double wDimY      = 0.6*mm;
  G4double wDimZ      = 0.6*mm;
  G4Material* wMat    = new G4Material("Universe", CLHEP::universe_mean_density, 1);
  wMat->AddElement(G4NistManager::Instance()->FindOrBuildElement(1), 1.0);
  G4Box*           sW = new G4Box ("Box",wDimX, wDimY, wDimZ);
  G4LogicalVolume* lW = new G4LogicalVolume(sW,wMat,"Box",0,0,0);
  G4PVPlacement*   pW = new G4PVPlacement(0,G4ThreeVector(),"Box",lW,0,false,0);
  //
  // --- Build all NIST materials and set a logical volume for each
//  const std::vector<G4String>& namesMat = G4NistManager::Instance()->GetNistMaterialNames();
  const G4int     numMat = g4Materials.size();
  std::vector<double> theCutInRange(numMat);
  std::vector<double> theCutInAbsLenght(numMat);
  std::vector<G4LogicalVolume*> theLogicVolumes(numMat);
  const G4double  halfX  =  0.5/numMat;  // half width of one material-box
  const G4double     x0  = -0.5+halfX;   // start x-position of the first material-box
  for (int im=0; im<numMat; ++im) {
    G4Material*      mat  = G4NistManager::Instance()->FindOrBuildMaterial(g4Materials[im]);
    if (mat==nullptr) {
      std::cerr << "  *** G4Setup::FakeG4Setup(): unknown Geant4 NIST material = " << g4Materials[im]
                << std::endl;
      exit(-1);
    }
    theCutInRange[im]     = ConvertCutToRange(electronCutInEnergy, mat);
    theCutInAbsLenght[im] = ConvertCutToAbsLength(gammaCutInEnergy, mat);
//    std::cout << mat->GetName() << " : " << theCutInAbsLenght[im] << " " << theCutInRange[im] << std::endl;
    G4Box*             ss = new G4Box ("Box", halfX, 0.5, 0.5);
    G4LogicalVolume*   ll = new G4LogicalVolume(ss, mat, "Box", 0, 0, 0);
    theLogicVolumes[im] = ll;
    new G4PVPlacement(0, G4ThreeVector(x0+im*halfX , 0, 0), "Box", ll, pW, false, 0);
  }
  //
  // --- Create particles that has secondary production threshold
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4Proton::Proton();
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();
  //
  for (G4int im=0; im<numMat; ++im) {
    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
//    std::cout << mat->GetName() << std::endl;
    G4ProductionCuts* iPcut = new G4ProductionCuts();
    iPcut->SetProductionCut(theCutInAbsLenght[im], 0); // set cut for gamma
    iPcut->SetProductionCut(theCutInRange[im], 1); // set cut for e-
    iPcut->SetProductionCut(1.0, 2); // set cut for e+
    iPcut->SetProductionCut(1.0, 3); // set cut for p+
    G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(mat, iPcut);
    couple->SetIndex(im+1);
    G4Region* iReg = new G4Region(mat->GetName()+"-region");
    iReg->AddRootLogicalVolume(theLogicVolumes[im]);
    iReg->UsedInMassGeometry(true);
    iReg->SetProductionCuts(iPcut);
    iReg->RegisterMaterialCouplePair(mat, couple);
  }
  // --- Update the couple tables
  G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  theCoupleTable->UpdateCoupleTable(pW);
  if (verbose>0) {
    theCoupleTable->DumpCouples();
  }
  //
}


// convert the secondary e- production cut in energy to lenght
double ConvertCutToRange(double electronCutInEnergy, const G4Material* theMaterial) {
  const int    numEbins     = 301;
  const double minCutEnergy = 100.0*CLHEP::eV;
  const double maxCutEnergy =  10.0*CLHEP::GeV;

  // compute the energy grid
  double delta = std::log(maxCutEnergy/minCutEnergy)/(numEbins-1);
  double base  = std::log(minCutEnergy)/delta;
  std::vector<double> theEnergyGrid(numEbins);
  theEnergyGrid[0]          = minCutEnergy;
  theEnergyGrid[numEbins-1] = maxCutEnergy;
  for (int ie = 1; ie<numEbins-1; ++ie) {
    theEnergyGrid[ie] = std::exp((base+ie)*delta);
  }

  // compute and approximated dedx for this material
  std::vector < double > theRange(numEbins);
  const double* elemDensity       = theMaterial->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* elemVect = theMaterial->GetElementVector();
  const int numElems = elemVect->size();
  double sum = 0.0;
  for (int ie = 0; ie<numEbins; ++ie) {
    const double ekin = theEnergyGrid[ie];
    double dedx       = 0.0;
    for (int iel=0; iel<numElems; ++iel) {
      const G4Element *elem = (*elemVect)[iel];
      double zet  = elem->GetZ();
      dedx += elemDensity[iel]*ComputeApproximatedDEDX(ekin, zet);
    }
    // integrate the 1/dEdx for getting the range
    const double dum0 = theEnergyGrid[ie]/dedx;
    if (ie == 0) {
      theRange[0] = dum0*delta;
      sum = 0.5*dum0;
    } else {
      sum += dum0;
      theRange[ie] = (sum-0.5*dum0)*delta;
    }
  }
  //
  if (electronCutInEnergy <= minCutEnergy) {
    return theRange[0];
  } else if (electronCutInEnergy >= maxCutEnergy) {
    return theRange[numEbins-1];
  }
  // set a Spline on the energy-lenght function then interpolate to get the cut
  // in range
  Spline* sp = new Spline(theEnergyGrid.data(), theRange.data(), numEbins, true);
  double cutInRange = sp->GetValueAt(electronCutInEnergy);
  delete sp;
  // cut in range [mm] that correspondes to the given cut in energy [MeV] in this material
  return cutInRange;
}

double ComputeApproximatedDEDX(double ekin, double zet) {
  const double cbr1  = 0.02;
  const double cbr2  = -5.7e-5;
  const double cbr3  = 1.;
  const double cbr4  = 0.072;
  const double tlow  = 10.0*CLHEP::keV;
  const double thigh =  1.0*CLHEP::GeV;
  const double mass  = CLHEP::electron_mass_c2;
  const double taul  = tlow/mass;
  const double cpot  = 1.6e-5*CLHEP::MeV;
  const double fact  = CLHEP::twopi*CLHEP::electron_mass_c2*CLHEP::classic_electr_radius*CLHEP::classic_electr_radius;

  double ionpot     = cpot*std::exp(0.9*std::log(zet))/mass;
  double ionpotlog  = std::log(ionpot);

  // calculate approximated dE/dx for electrons
  double tau  = ekin/mass;
  double dEdx = 0.0;
  if (tau<taul) {
    double t1    = taul+1.;
    double t2    = taul+2.;
    double tsq   = taul*taul;
    double beta2 = taul*t2/(t1*t1);
    double f     = 1.-beta2+std::log(tsq/2.)+(0.5+0.25*tsq+(1.+2.*taul)*std::log(0.5))/(t1*t1);
    dEdx         = (std::log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx         = fact*zet*dEdx;
    double clow  = dEdx*std::sqrt(taul);
    dEdx         = clow/std::sqrt(tau);
  } else {
    double t1    = tau+1.;
    double t2    = tau+2.;
    double tsq   = tau*tau;
    double beta2 = tau*t2/(t1*t1);
    double f     = 1.-beta2+std::log(tsq/2.)+(0.5+0.25*tsq+(1.+2.*tau)*std::log(0.5))/(t1*t1);
    dEdx         = (std::log(2.*tau+4.)-2.*ionpotlog+f)/beta2;
    dEdx         = fact*zet*dEdx;
    // loss from bremsstrahlung follows
    double cbrem = (cbr1+cbr2*zet)*(cbr3+cbr4*std::log(ekin/thigh));
    cbrem        = 0.1*zet*(zet+1.)*cbrem*tau/beta2;
    dEdx        += fact*cbrem;
  }
  return dEdx;
}


// convert the secondary e- production cut in energy to lenght
double ConvertCutToAbsLength(double gammaCutInEnergy, const G4Material* theMaterial) {
  const int    numEbins     = 301;
  const double minCutEnergy = 100.0*CLHEP::eV;
  const double maxCutEnergy =  10.0*CLHEP::GeV;

  // compute the energy grid
  double delta = std::log(maxCutEnergy/minCutEnergy)/(numEbins-1);
  double base  = std::log(minCutEnergy)/delta;
  std::vector<double> theEnergyGrid(numEbins);
  theEnergyGrid[0]          = minCutEnergy;
  theEnergyGrid[numEbins-1] = maxCutEnergy;
  for (int ie = 1; ie<numEbins-1; ++ie) {
    theEnergyGrid[ie] = std::exp((base+ie)*delta);
  }

  // compute and approximated absorption length table for thismaterial
  std::vector < double > theAbsLength(numEbins);
  const double* elemDensity       = theMaterial->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* elemVect = theMaterial->GetElementVector();
  const int numElems   = elemVect->size();
  double maxAbsLenght  = -1.0;
  int    indxMaxLength = -1;
  for (int ie = 0; ie<numEbins; ++ie) {
    const double ekin = theEnergyGrid[ie];
    double macXsec    = 0.0;
    for (int iel=0; iel<numElems; ++iel) {
      const G4Element *elem = (*elemVect)[iel];
      double zet  = elem->GetZ();
      macXsec += elemDensity[iel]*ComputeApproximatedAbsXsecPerAtom(ekin, zet);
    }
    double absLength = 5.0/macXsec;
    theAbsLength[ie] = absLength;
    if (maxAbsLenght < absLength) {
      maxAbsLenght   = absLength;
      indxMaxLength  = ie;
    }
  }
  //
  if (gammaCutInEnergy <= minCutEnergy) {
    return theAbsLength[0];
  } else if (gammaCutInEnergy >= maxCutEnergy) {
    return theAbsLength[numEbins-1];
  }
  int istart  = 0;
  int numData = numEbins;
  if (gammaCutInEnergy < theEnergyGrid[indxMaxLength]) {
    numData = indxMaxLength+1;
  } else {
    istart  = indxMaxLength;
    numData = numEbins-indxMaxLength;
  }
  // set a Spline on the energy-lenght function then interpolate to get the cut
  // in range
  Spline* sp = new Spline(&(theEnergyGrid[istart]), &(theAbsLength[istart]), numData, true);
  double cutInAbsLenght = sp->GetValueAt(gammaCutInEnergy);
  delete sp;
  // cut in absorption length [mm] that correspondes to the given cut in energy [MeV] in this material
  return cutInAbsLenght;
}


// Compute the photon "absorption" cross section: sum of destructive (approximated) cross sections like
// pair production, Compton scattering and photoelectric effect
double ComputeApproximatedAbsXsecPerAtom(double ekin, double zet) {
  const double t1keV      =   1.0*CLHEP::keV;
  const double t200keV    = 200.0*CLHEP::keV;
  const double t100MeV    = 100.0*CLHEP::MeV;
  const double Zsquare    = zet*zet;
  const double Zlog       = std::log(zet);
  const double Zlogsquare = Zlog*Zlog;
  // set some Z dependent variables
  double s200keV = (0.2651-0.1501*Zlog+0.02283*Zlogsquare)*Zsquare;
  double tmin    = (0.552+218.5/zet+557.17/Zsquare)*CLHEP::MeV;
  double smin    = (0.01239+0.005585*Zlog-0.000923*Zlogsquare)*std::exp(1.5*Zlog);
  double cmin    = std::log(s200keV/smin)/(std::log(tmin/t200keV)*std::log(tmin/t200keV));
  double tlow    = 0.2*std::exp(-7.355/std::sqrt(zet))*CLHEP::MeV;
  double slow    = s200keV*std::exp(0.042*zet*std::log(t200keV/tlow)*std::log(t200keV/tlow));
  double s1kev   = 300.0*Zsquare;
  double clow    = std::log(s1kev/slow)/std::log(tlow/t1keV);
  double chigh   = (7.55e-5-0.0542e-5*zet)*Zsquare*zet/std::log(t100MeV/tmin);
  // calculate the absorption cross section (using an approximate empirical formula)
  double xs = 0.0;
  if (ekin<tlow) {
    if (ekin<t1keV)
      xs = slow*std::exp(clow*std::log(tlow/t1keV));
    else
      xs = slow*std::exp(clow*std::log(tlow/ekin));
  } else if (ekin<t200keV) {
    xs = s200keV*std::exp(0.042*zet*std::log(t200keV/ekin)*std::log(t200keV/ekin));
  } else if (ekin<tmin) {
    double dum = std::log(tmin/ekin);
    xs = smin*std::exp(cmin*dum*dum);
  } else {
    xs = smin+chigh*std::log(ekin/tmin);
  }
  return xs*CLHEP::barn;
}
