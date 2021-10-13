//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     DPMMaterialBuilder
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
// 31-10-05 Add chemical effect and gas properties (V.Ivanchenko)
// 27.02.06 V.Ivanchenko add ConstructNewGasMaterial
// 11.05.06 V.Ivanchenko add warning flag to FindMaterial method
// 27.06.06 V.Ivanchenko fix graphite description
// 27.07.07 V.Ivanchenko remove dependence on NistManager
// 30.10.09 V.Ivanchenko update density of G4_GRAFITE from PDG'2008
//                       added G4_GRAPHITE_POROUS
// 03.11.09 A.Lechner changed following material names:
//                    From G4_NYLON-6/6 to G4_NYLON-6-6
//                    From G4_NYLON-6/10 to G4_NYLON-6-10
// 12.12.10 A.Ivantchenko added following materials methodes:
//                    BioChemicalMaterials() and SpaceMaterials(),
//                    where new materials are introduced
// 14.06.11 A.Ivantchenko updated body materials (G4_....ICRP)
//                    according ICRU Report 46 (1992) instead of 1975 
//                    data from ICRU Report 37 used previously
// 26.10.11 new scheme for G4Exception  (mma)
// 09.02.12 P.Gumplinger add ConstructNewIdealGasMaterial
// 30.04.13 M.Trocme & S.Seltzer:
//        - Replace AddElementByWeightFraction() by AddElementByAtomCount() 
//          as much as possible
//        - Comment out ill-defined material GLUCOSE 
//        - Fixed density and atom composition of  POLYCHLOROSTYRENE, 
//          POLYVINYL_BUTYRAL, TERPHENYL
// -------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DPMMaterialBuilder.hh"
#include "DPMElementBuilder.hh"
#include "G4Element.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

#ifdef G4MULTITHREADED
G4Mutex DPMMaterialBuilder::nistMaterialMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DPMMaterialBuilder::DPMMaterialBuilder(DPMElementBuilder* eb, G4int vb)
: elmBuilder(eb),
  verbose(vb),
  nMaterials(0),
  nComponents(0),
  nCurrent(0)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DPMMaterialBuilder::~DPMMaterialBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::FindOrBuildMaterial(const G4String& matname,
                                                       G4bool, G4bool warning)
{
  if(verbose > 1) {
    G4cout << "DPMMaterialBuilder::FindOrBuildMaterial " 
	   << matname << G4endl;
  }
  G4Material* mat = FindMaterial(matname);
  if(mat != nullptr) { return mat; }
  G4String name = matname;
  if(name == "G4_NYLON-6/6" || name == "G4_NYLON-6/10") {
    if("G4_NYLON-6/6" == matname)  { name = "G4_NYLON-6-6"; }
    else { name = "G4_NYLON-6-10";}
    mat = FindMaterial(name);
  }
  return (mat == nullptr) ? BuildNistMaterial(name, warning) : mat; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::BuildNistMaterial(const G4String& name, 
						     G4bool warning)
{
  G4Material* mat = nullptr;
  // Check if name inside DB
  for (G4int i=0; i<nMaterials; ++i) {

    if (name == names[i]) {
#ifdef G4MULTITHREADED
      G4MUTEXLOCK(&nistMaterialMutex);
#endif
      if(matIndex[i] == -1) { 
	// Build new Nist material 
	mat = BuildMaterial(i); 
      } else { 
	// Nist material was already built
	const G4MaterialTable* theMaterialTable = 
	  G4Material::GetMaterialTable();
	mat = (*theMaterialTable)[matIndex[i]]; 
      }
#ifdef G4MULTITHREADED
      G4MUTEXUNLOCK(&nistMaterialMutex);
#endif
      return mat;
    }
  }

  if( (verbose == 1 && warning) || verbose > 1) {
    G4cout << "DPMMaterialBuilder::FindOrBuildMaterial WARNING:"
	   << " material <" << name
	   << "> is not found." << G4endl;
  }
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* 
DPMMaterialBuilder::FindOrBuildSimpleMaterial(G4int Z, G4bool warn)
{
  G4Material* mat = FindSimpleMaterial(Z);
  if(mat == nullptr) {
    mat = BuildNistMaterial(names[Z], warn); 
  }
  return mat;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::BuildMaterial(G4int i)
{
  if (verbose > 1) {
    G4cout << "DPMMaterialBuilder: BuildMaterial #" << i
	   << G4endl;
  }
  G4Material* mat = 0;
  if(i >= nMaterials) { return mat; }

  G4int nc = components[i];

  // Check gas parameters:
  // defaults may be changed via AddGas() method
  G4double t = NTP_Temperature;
  G4double p = STP_Pressure;
  if(kStateGas == states[i]) {
    size_t nn = idxGas.size();
    if(nn > 0) {
      for(size_t j=0; j<nn; ++j) {
        if(i == idxGas[j]) {
	  t = gasTemperature[j];
          p = gasPressure[j];
          break;
	}
      }
    }
  }
    // liquids
  //  } else if( !STP[i] ) { t = 0.0; }

  mat = new G4Material(names[i],densities[i],nc,states[i],t,p);

  if (verbose>1) { G4cout << "New material nComponents= " << nc << G4endl; }
  if (nc > 0) {
    G4int idx = indexes[i];
    for (G4int j=0; j<nc; ++j) {
      G4int Z = elements[idx+j];
      G4Element* el = elmBuilder->FindOrBuildElement(Z);
      if(!el) {
	G4cout << "DPMMaterialBuilder::BuildMaterial:"
	       << "  ERROR: elements Z= " << Z << " is not found"
	       << " for material " << names[i]
	       << G4endl;
	G4Exception("DPMMaterialBuilder::BuildMaterial()", "mat103",
	             FatalException, "Failed to construct material");
	return 0;
      }
      if(atomCount[i]) {
	mat->AddElement(el,G4lrint(fractions[idx+j]));
      } else {
	mat->AddElement(el,fractions[idx+j]);
      }
    }
  }

  // Ionisation potential can be defined via NIST DB or 
  // Chemical Formula (ICRU37 Report data)
  G4IonisParamMat* ion = mat->GetIonisation();
  G4double exc0 = ion->GetMeanExcitationEnergy();
  G4double exc1 = exc0;
  if(chFormulas[i] != "") {
    mat->SetChemicalFormula(chFormulas[i]);
    exc1 = ion->FindMeanExcitationEnergy(mat);
  }
  // If exists, NIST DB data always overwrites other data 
  if(ionPotentials[i] > 0.0) { exc1 = ionPotentials[i]; }
  if(exc0 != exc1) { ion->SetMeanExcitationEnergy(exc1); }

  // Index in Material Table
  matIndex[i] = mat->GetIndex();
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens, 
				      G4bool,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) { 
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return 0;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && 
     (temp != NTP_Temperature || pres != CLHEP::STP_Pressure))
    { stp = false; }

  AddMaterial(name,dens*CLHEP::cm3/CLHEP::g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
  }

  return BuildMaterial(nMaterials-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& w,
				      G4double dens, 
				      G4bool,
				      G4State state,     
				      G4double temp,  
				      G4double pres)
{
  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) { 
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
	   << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return 0;
  } 

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(state == kStateGas && 
     (temp != NTP_Temperature || pres != CLHEP::STP_Pressure))
    { stp = false; }
  AddMaterial(name,dens*CLHEP::cm3/CLHEP::g,0,0.,els,state,stp);
  if(!stp) { AddGas(name,temp,pres); }

  for (G4int i=0; i<els; ++i) {
    AddElementByWeightFraction(elmBuilder->GetZ(elm[i]), w[i]);
  }
  
  return BuildMaterial(nMaterials-1); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::ConstructNewGasMaterial(
				      const G4String& name,
				      const G4String& nameDB,
				      G4double temp, 
				      G4double pres, 
				      G4bool)
{
  // Material name is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) { 
    G4cout << "DPMMaterialBuilder::ConstructNewGasMaterial:"
           << "  WARNING: the material <" << name
	   << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
	   << G4endl;
    return mat; 
  }

  G4Material* bmat = FindOrBuildMaterial(nameDB);
  if(!bmat) {
    G4cout << "DPMMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING: the Name <" << nameDB 
	   << "> is NOT in the database: no new gas will be constructed."
	   << G4endl;
    return 0;
  }
  if(bmat->GetState() != kStateGas) {
    G4cout << "DPMMaterialBuilder::ConstructNewGasMaterial:"
	   << "  WARNING:  <" << nameDB 
	   << "> is NOT a gas -  no new gas will be constructed."
	   << G4endl;
    return 0;
  }

  G4double dens = bmat->GetDensity()*pres*bmat->GetTemperature()
    /(temp*bmat->GetPressure());
  mat = new G4Material(name,dens,bmat,kStateGas,temp,pres);

  if (verbose>1) {
    G4cout << "DPMMaterialBuilder::ConstructNewGasMaterial: done" << G4endl;
    G4cout << &mat << G4endl; 
  }	
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DPMMaterialBuilder::ConstructNewIdealGasMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
                                      G4bool,
                                      G4double temp,
                                      G4double pres)
{
  G4State state = kStateGas;

  // Material is in DB
  G4Material* mat = FindOrBuildMaterial(name);
  if(mat) {
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: the material <" << name
           << "> already exists." << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return mat;
  }

  // Material not in DB
  G4int els = elm.size();
  if(els == 0) {
    G4cout << "DPMMaterialBuilder::ConstructNewMaterial:"
           << "  WARNING: empty list of elements for " << name
           << G4endl;
    G4cout << "      New material will NOT be built!"
           << G4endl;
    return 0;
  }

  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential is not defined
  G4bool stp = true;
  if(temp != NTP_Temperature || pres != CLHEP::STP_Pressure)
    { stp = false; }

  G4double massPerMole = 0.;

  G4int Z = 0;
  for (G4int i=0; i<els; ++i) {
    Z = elmBuilder->GetZ(elm[i]);
    massPerMole += nbAtoms[i]*elmBuilder->GetAtomicMassAmu(Z)*CLHEP::amu_c2;
  }

  G4double dens = massPerMole/(CLHEP::Avogadro*CLHEP::k_Boltzmann*temp/pres);

  if (els == 1) { AddMaterial(name,dens,Z,0.,els,state,stp); }
  else {
    AddMaterial(name,dens,0,0.,els,state,stp);
    for (G4int i=0; i<els; ++i) {
      AddElementByAtomCount(elmBuilder->GetZ(elm[i]), nbAtoms[i]);
    }
  }

  if(!stp) { AddGas(name,temp,pres); }

  return BuildMaterial(nMaterials-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::AddMaterial(const G4String& nameMat, G4double dens,
					G4int Z, G4double pot, 
					G4int ncomp, G4State state, 
					G4bool stp)
{
  // add parameters of material into internal vectors
  // density in g/cm3, mean ionisation potential in eV

  // if ncomp == 1 then Z should be defined and 
  // AddElement should not be applied

  if (nCurrent != 0) {
    G4cout << "DPMMaterialBuilder::AddMaterial WARNING: previous "
	   << "mixture " << nMaterials << " " << names[nMaterials]
	   << " is not yet complete!"
	   << G4endl;
    G4cout << "         New material " << nameMat << " will not be added."
	   << G4endl;
    return;
  }

  // density in g/cm3, mean ionisation potential in eV

  names.push_back(nameMat);
  chFormulas.push_back("");
  densities.push_back(dens*CLHEP::g/CLHEP::cm3);
  ionPotentials.push_back(pot*CLHEP::eV);
  states.push_back(state);
  components.push_back(ncomp);
  indexes.push_back(nComponents);
  STP.push_back(stp);
  matIndex.push_back(-1);
  atomCount.push_back(false);

  if (1 == ncomp && Z > 0) {
    elements.push_back(Z);
    fractions.push_back(1.0);
    atomCount[nMaterials] = true;
    ++nComponents;
    nCurrent = 0;
  } else {
    nCurrent = ncomp;
  }

  ++nMaterials;

  if(verbose > 1) {
    G4cout << "New material " << nameMat << " is prepared; "
           << " nMaterials= " << nMaterials
           << " nComponents= " << nComponents
           << " nCurrent= " << nCurrent
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::SetVerbose(G4int val)
{
  verbose = val;
  elmBuilder->SetVerbose(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::ListMaterials(const G4String& mnam) const
{
  if (mnam == "simple")           { ListNistSimpleMaterials(); }
  else if (mnam == "compound")    { ListNistCompoundMaterials(); }

  else if (mnam == "all") {
    ListNistSimpleMaterials();
    ListNistCompoundMaterials();

  } else {
    G4cout << "### DPMMaterialBuilder::ListMaterials: Warning " 
	   << mnam << " list is not known." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::ListNistSimpleMaterials() const
{
  G4cout << "=======================================================" << G4endl;
  G4cout << "###   Simple Materials from the NIST Data Base      ###" << G4endl;
  G4cout << "=======================================================" << G4endl;
  G4cout << " Z   Name   density(g/cm^3)  I(eV)                     " << G4endl;
  G4cout << "=======================================================" << G4endl;
  for (G4int i=1; i<nElementary; ++i) {DumpElm(i);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::ListNistCompoundMaterials() const
{
  G4cout << "=============================================================" << G4endl;
  G4cout << "###    Compound Materials from the NIST Data Base          ##" << G4endl;
  G4cout << "=============================================================" << G4endl;
  G4cout << " Ncomp             Name      density(g/cm^3)  I(eV) ChFormula" << G4endl;
  G4cout << "=============================================================" << G4endl;
  for (G4int i=nElementary; i<nNIST; ++i) {DumpMix(i);}
  DumpMix(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::DumpElm(G4int i) const
{
  G4cout << std::setw(2)  << i << " " 
         << std::setw(6)  << names[i] 
         << std::setw(14) << densities[i]*cm3/g 
         << std::setw(11) << ionPotentials[i]/eV
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::DumpMix(G4int i) const
{
  G4int nc = components[i];
  G4cout << std::setw(2)  << nc << " "
	 << std::setw(26) << names[i] << " " 
         << std::setw(10) << densities[i]*cm3/g 
	 << std::setw(10) << ionPotentials[i]/eV
	 << "   " << chFormulas[i]
	 << G4endl;
  if (nc > 1) {
    G4int imin = indexes[i];
    G4int imax = imin + nc;
    for (G4int j=imin; j<imax; ++j) {
      G4cout << std::setw(10) << elements[j] << std::setw(14) << fractions[j] 
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
DPMMaterialBuilder::AddGas(const G4String& nameMat, G4double t, G4double p)
{
  for(G4int i=0; i<nMaterials; ++i) {
    if(nameMat == names[i]) {
      idxGas.push_back(i);
      gasTemperature.push_back(t);
      gasPressure.push_back(p);
      return;
    }
  }
  G4cout << "WARNING: DPMMaterialBuilder::AddGas problem: there is no "
	 << nameMat << " in the list of materials."
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::AddElementByWeightFraction(G4int Z, G4double w)
{
  elements.push_back(Z);
  fractions.push_back(w);
  --nCurrent;
  ++nComponents;
  if (nCurrent == 0) {
    G4int n = nMaterials - 1;
    G4double sum = 0.0;
    G4int imin = indexes[n];
    G4int imax = imin + components[n];

    if(!atomCount[n]) {
      for(G4int i=imin; i<imax; ++i) {sum += fractions[i];}
      if (sum > 0.0) for (G4int i=imin; i<imax; ++i) {fractions[i] /= sum;}
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::AddElementByWeightFraction(const G4String& name,
                                                       G4double w)
{
  G4int Z = elmBuilder->GetZ(name);
  AddElementByWeightFraction(Z, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::AddElementByAtomCount(G4int Z, G4int nb)
{
  atomCount[nMaterials-1] = true;
  G4double w = (G4double)nb;
  AddElementByWeightFraction(Z, w);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::AddElementByAtomCount(const G4String& name,
						  G4int nb)
{
  atomCount[nMaterials-1] = true;
  G4int Z = elmBuilder->GetZ(name);
  G4double w = (G4double)nb;
  AddElementByWeightFraction(Z, w);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::Initialise()
{
  if (verbose > 0) {
    G4cout << "### DPMMaterialBuilder::Initialise()" << G4endl;
  }
  NistSimpleMaterials();
  NistCompoundMaterials();
  NistCompoundMaterials2();
  NistCompoundMaterials3();

  if (verbose > 1) { ListMaterials("all"); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::NistSimpleMaterials()
{
  // density in g/cm3, mean ionisation potential in eV

  AddMaterial("G4_WATER", 1.0,0, 78., 2);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("O" ,  1);
  chFormulas[nMaterials-1] = "H_2O";

  AddMaterial("G4_H" ,  8.37480e-5,  1,  19.2, 1, kStateGas);
  AddMaterial("G4_He",  1.66322e-4,  2,  41.8, 1, kStateGas);
  AddMaterial("G4_Li",  0.534     ,  3,  40. );
  AddMaterial("G4_Be",  1.848     ,  4,  63.7);
  AddMaterial("G4_B" ,  2.37      ,  5,  76. );
  AddMaterial("G4_C" ,  2.        ,  6,  81. );
  AddMaterial("G4_N" ,  1.16520e-3,  7,  82. , 1, kStateGas);
  AddMaterial("G4_O" ,  1.33151e-3,  8,  95. , 1, kStateGas);
  AddMaterial("G4_F" ,  1.58029e-3,  9, 115. , 1, kStateGas);
  AddMaterial("G4_Ne",  8.38505e-4, 10, 137. , 1, kStateGas);
  AddMaterial("G4_Na",  0.971     , 11, 149. );
  AddMaterial("G4_Mg",  1.74      , 12, 156. );
  AddMaterial("G4_Al",  2.699     , 13, 166. );
  AddMaterial("G4_Si",  2.33      , 14, 173. );
  AddMaterial("G4_P" ,  2.2       , 15, 173. );
  AddMaterial("G4_S" ,  2.0       , 16, 180. );
  AddMaterial("G4_Cl",  2.99473e-3, 17, 174. , 1, kStateGas);
  AddMaterial("G4_Ar",  1.66201e-3, 18, 188.0, 1, kStateGas);
  AddMaterial("G4_K" ,  0.862     , 19, 190. );
  AddMaterial("G4_Ca",  1.55      , 20, 191. );
  AddMaterial("G4_Sc",  2.989     , 21, 216. );
  AddMaterial("G4_Ti",  4.54      , 22, 233. );
  AddMaterial("G4_V" ,  6.11      , 23, 245. );
  AddMaterial("G4_Cr",  7.18      , 24, 257. );
  AddMaterial("G4_Mn",  7.44      , 25, 272. );
  AddMaterial("G4_Fe",  7.874     , 26, 286. );
  AddMaterial("G4_Co",  8.9       , 27, 297. );
  AddMaterial("G4_Ni",  8.902     , 28, 311. );
  AddMaterial("G4_Cu",  8.96      , 29, 322. );
  AddMaterial("G4_Zn",  7.133     , 30, 330. );
  AddMaterial("G4_Ga",  5.904     , 31, 334. );
  AddMaterial("G4_Ge",  5.323     , 32, 350. );
  AddMaterial("G4_As",  5.73      , 33, 347. );
  AddMaterial("G4_Se",  4.5       , 34, 348. );
  AddMaterial("G4_Br",  7.07210e-3, 35, 343. , 1, kStateGas);
  AddMaterial("G4_Kr",  3.47832e-3, 36, 352. , 1, kStateGas);
  AddMaterial("G4_Rb",  1.532     , 37, 363. );
  AddMaterial("G4_Sr",  2.54      , 38, 366. );
  AddMaterial("G4_Y" ,  4.469     , 39, 379. );
  AddMaterial("G4_Zr",  6.506     , 40, 393. );
  AddMaterial("G4_Nb",  8.57      , 41, 417. );
  AddMaterial("G4_Mo", 10.22      , 42, 424. );
  AddMaterial("G4_Tc", 11.50      , 43, 428. );
  AddMaterial("G4_Ru", 12.41      , 44, 441. );
  AddMaterial("G4_Rh", 12.41      , 45, 449. );
  AddMaterial("G4_Pd", 12.02      , 46, 470. );
  AddMaterial("G4_Ag", 10.5       , 47, 470. );
  AddMaterial("G4_Cd",  8.65      , 48, 469. );
  AddMaterial("G4_In",  7.31      , 49, 488. );
  AddMaterial("G4_Sn",  7.31      , 50, 488. );
  AddMaterial("G4_Sb",  6.691     , 51, 487. );
  AddMaterial("G4_Te",  6.24      , 52, 485. );
  AddMaterial("G4_I" ,  4.93      , 53, 491. );
  AddMaterial("G4_Xe",  5.48536e-3, 54, 482. , 1, kStateGas);
  AddMaterial("G4_Cs",  1.873     , 55, 488. );
  AddMaterial("G4_Ba",  3.5       , 56, 491. );
  AddMaterial("G4_La",  6.154     , 57, 501. );
  AddMaterial("G4_Ce",  6.657     , 58, 523. );
  AddMaterial("G4_Pr",  6.71      , 59, 535. );
  AddMaterial("G4_Nd",  6.9       , 60, 546. );
  AddMaterial("G4_Pm",  7.22      , 61, 560. );
  AddMaterial("G4_Sm",  7.46      , 62, 574. );
  AddMaterial("G4_Eu",  5.243     , 63, 580. );
  AddMaterial("G4_Gd",  7.9004    , 64, 591. );
  AddMaterial("G4_Tb",  8.229     , 65, 614. );
  AddMaterial("G4_Dy",  8.55      , 66, 628. );
  AddMaterial("G4_Ho",  8.795     , 67, 650. );
  AddMaterial("G4_Er",  9.066     , 68, 658. );
  AddMaterial("G4_Tm",  9.321     , 69, 674. );
  AddMaterial("G4_Yb",  6.73      , 70, 684. );
  AddMaterial("G4_Lu",  9.84      , 71, 694. );
  AddMaterial("G4_Hf", 13.31      , 72, 705. );
  AddMaterial("G4_Ta", 16.654     , 73, 718. );
  AddMaterial("G4_W" , 19.30      , 74, 727. );
  AddMaterial("G4_Re", 21.02      , 75, 736. );
  AddMaterial("G4_Os", 22.57      , 76, 746. );
  AddMaterial("G4_Ir", 22.42      , 77, 757. );
  AddMaterial("G4_Pt", 21.45      , 78, 790. );
  AddMaterial("G4_Au", 19.32      , 79, 790. );
  AddMaterial("G4_Hg", 13.546     , 80, 800. );
  AddMaterial("G4_Tl", 11.72      , 81, 810. );
  AddMaterial("G4_Pb", 11.35      , 82, 823. );
  AddMaterial("G4_Bi",  9.747     , 83, 823. );
  AddMaterial("G4_Po",  9.32      , 84, 830. );
  AddMaterial("G4_At",  9.32      , 85, 825. );
  AddMaterial("G4_Rn",  9.00662e-3, 86, 794. , 1, kStateGas);
  AddMaterial("G4_Fr",  1.00      , 87, 827. );
  AddMaterial("G4_Ra",  5.00      , 88, 826. );
  AddMaterial("G4_Ac", 10.07      , 89, 841. );
  AddMaterial("G4_Th", 11.72      , 90, 847. );
  AddMaterial("G4_Pa", 15.37      , 91, 878. );
  AddMaterial("G4_U" , 18.95      , 92, 890. );
  AddMaterial("G4_Np", 20.25      , 93, 902. );
  AddMaterial("G4_Pu", 19.84      , 94, 921. );
  AddMaterial("G4_Am", 13.67      , 95, 934. );
  AddMaterial("G4_Cm", 13.51      , 96, 939. );
  AddMaterial("G4_Bk", 14.00      , 97, 952. );
  AddMaterial("G4_Cf", 10.00      , 98, 966. );

  nElementary = nMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DPMMaterialBuilder::NistCompoundMaterials()
{
  AddMaterial("G4_A-150_TISSUE", 1.127, 0, 65.1, 6);
  AddElementByWeightFraction( 1, 0.101327);
  AddElementByWeightFraction( 6, 0.775501);
  AddElementByWeightFraction( 7, 0.035057);
  AddElementByWeightFraction( 8, 0.052316);
  AddElementByWeightFraction( 9, 0.017422);
  AddElementByWeightFraction(20, 0.018378);

  AddMaterial("G4_ADIPOSE_TISSUE_ICRP", 0.95, 0, 63.2, 7);
  AddElementByWeightFraction( 1, 0.114);
  AddElementByWeightFraction( 6, 0.598);
  AddElementByWeightFraction( 7, 0.007);
  AddElementByWeightFraction( 8, 0.278);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(16, 0.001);
  AddElementByWeightFraction(17, 0.001);

  AddMaterial("G4_AIR", 0.00120479, 0, 85.7, 4, kStateGas);
  AddElementByWeightFraction( 6, 0.000124);
  AddElementByWeightFraction( 7, 0.755267);
  AddElementByWeightFraction( 8, 0.231781);
  AddElementByWeightFraction(18, 0.012827);

  AddMaterial("G4_B-100_BONE", 1.45, 0, 85.9, 6);
  AddElementByWeightFraction( 1, 0.065471);
  AddElementByWeightFraction( 6, 0.536945);
  AddElementByWeightFraction( 7, 0.0215  );
  AddElementByWeightFraction( 8, 0.032085);
  AddElementByWeightFraction( 9, 0.167411);
  AddElementByWeightFraction(20, 0.176589);

  AddMaterial("G4_BLOOD_ICRP", 1.06, 0, 75.2, 10);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.110);
  AddElementByWeightFraction( 7, 0.033);
  AddElementByWeightFraction( 8, 0.745);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.002);
  AddElementByWeightFraction(26, 0.001);

  AddMaterial("G4_BONE_COMPACT_ICRU", 1.85, 0, 91.9, 8);
  AddElementByWeightFraction( 1, 0.064);
  AddElementByWeightFraction( 6, 0.278);
  AddElementByWeightFraction( 7, 0.027);
  AddElementByWeightFraction( 8, 0.410);
  AddElementByWeightFraction(12, 0.002);
  AddElementByWeightFraction(15, 0.07 );
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(20, 0.147);

  // Sceleton Cortical bone for Adult ICRU 46
  AddMaterial("G4_BONE_CORTICAL_ICRP", 1.92, 0, 110, 9);
  AddElementByWeightFraction( 1, 0.034);
  AddElementByWeightFraction( 6, 0.155);
  AddElementByWeightFraction( 7, 0.042);
  AddElementByWeightFraction( 8, 0.435);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(12, 0.002);
  AddElementByWeightFraction(15, 0.103);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(20, 0.225);

  AddMaterial("G4_BRAIN_ICRP", 1.04, 0, 73.3, 9);
  AddElementByWeightFraction( 1, 0.107);
  AddElementByWeightFraction( 6, 0.145);
  AddElementByWeightFraction( 7, 0.022);
  AddElementByWeightFraction( 8, 0.712);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.004);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.003);

  AddMaterial("G4_EYE_LENS_ICRP", 1.07, 0, 73.3, 8);
  AddElementByWeightFraction( 1, 0.096);
  AddElementByWeightFraction( 6, 0.195);
  AddElementByWeightFraction( 7, 0.057);
  AddElementByWeightFraction( 8, 0.646);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.001);

}

void DPMMaterialBuilder::NistCompoundMaterials2()
{
  //Adult Lung congested
  AddMaterial("G4_LUNG_ICRP", 1.04, 0, 75.3, 9);
  AddElementByWeightFraction( 1, 0.105);
  AddElementByWeightFraction( 6, 0.083);
  AddElementByWeightFraction( 7, 0.023);
  AddElementByWeightFraction( 8, 0.779);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.002);

  AddMaterial("G4_MS20_TISSUE", 1.0, 0, 75.1, 6);
  AddElementByWeightFraction( 1, 0.081192);
  AddElementByWeightFraction( 6, 0.583442);
  AddElementByWeightFraction( 7, 0.017798);
  AddElementByWeightFraction( 8, 0.186381);
  AddElementByWeightFraction(12, 0.130287);
  AddElementByWeightFraction(17, 0.0009  );

  AddMaterial("G4_MUSCLE_SKELETAL_ICRP", 1.05, 0, 75.3, 9);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.143);
  AddElementByWeightFraction( 7, 0.034);
  AddElementByWeightFraction( 8, 0.710);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.001);
  AddElementByWeightFraction(19, 0.004);

  // from old ICRU report
  AddMaterial("G4_MUSCLE_STRIATED_ICRU", 1.04, 0, 74.7, 8);
  AddElementByWeightFraction( 1, 0.102);
  AddElementByWeightFraction( 6, 0.123);
  AddElementByWeightFraction( 7, 0.035);
  AddElementByWeightFraction( 8, 0.729);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.004);
  AddElementByWeightFraction(19, 0.003);

  AddMaterial("G4_MUSCLE_WITH_SUCROSE", 1.11, 0, 74.3, 4);
  AddElementByWeightFraction( 1, 0.098234);
  AddElementByWeightFraction( 6, 0.156214);
  AddElementByWeightFraction( 7, 0.035451);
  AddElementByWeightFraction( 8, 0.7101  );
  
  AddMaterial("G4_MUSCLE_WITHOUT_SUCROSE", 1.07, 0, 74.2, 4);
  AddElementByWeightFraction( 1, 0.101969);
  AddElementByWeightFraction( 6, 0.120058);
  AddElementByWeightFraction( 7, 0.035451);
  AddElementByWeightFraction( 8, 0.742522);

  AddMaterial("G4_POLYETHYLENE", 0.94, 0, 57.4, 2);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  2);
  chFormulas[nMaterials-1] = "(C_2H_4)_N-Polyethylene";

  AddMaterial("G4_SKIN_ICRP", 1.09, 0, 72.7, 9);
  AddElementByWeightFraction( 1, 0.100);
  AddElementByWeightFraction( 6, 0.204);
  AddElementByWeightFraction( 7, 0.042);
  AddElementByWeightFraction( 8, 0.645);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.003);
  AddElementByWeightFraction(19, 0.001);

  AddMaterial("G4_TESTIS_ICRP", 1.04, 0, 75., 9);
  AddElementByWeightFraction( 1, 0.106);
  AddElementByWeightFraction( 6, 0.099);
  AddElementByWeightFraction( 7, 0.020);
  AddElementByWeightFraction( 8, 0.766);
  AddElementByWeightFraction(11, 0.002);
  AddElementByWeightFraction(15, 0.001);
  AddElementByWeightFraction(16, 0.002);
  AddElementByWeightFraction(17, 0.002);
  AddElementByWeightFraction(19, 0.002);

  // TISSUE_SOFT_MALE ICRU-44/46 (1989)
  AddMaterial("G4_TISSUE_SOFT_ICRP", 1.03, 0, 72.3, 9);
  AddElementByWeightFraction( 1, 0.105);
  AddElementByWeightFraction( 6, 0.256);
  AddElementByWeightFraction( 7, 0.027);
  AddElementByWeightFraction( 8, 0.602);
  AddElementByWeightFraction(11, 0.001);
  AddElementByWeightFraction(15, 0.002);
  AddElementByWeightFraction(16, 0.003);
  AddElementByWeightFraction(17, 0.002);
  AddElementByWeightFraction(19, 0.002);

  // Tissue soft adult ICRU-33 (1980)
  AddMaterial("G4_TISSUE_SOFT_ICRU-4", 1.0, 0, 74.9, 4);
  AddElementByWeightFraction( 1, 0.101);
  AddElementByWeightFraction( 6, 0.111);
  AddElementByWeightFraction( 7, 0.026);
  AddElementByWeightFraction( 8, 0.762);

  AddMaterial("G4_TISSUE-METHANE", 0.00106409, 0, 61.2, 4, kStateGas);
  AddElementByWeightFraction( 1, 0.101869);
  AddElementByWeightFraction( 6, 0.456179);
  AddElementByWeightFraction( 7, 0.035172);
  AddElementByWeightFraction( 8, 0.40678 );

  AddMaterial("G4_TISSUE-PROPANE", 0.00182628, 0, 59.5, 4, kStateGas);
  AddElementByWeightFraction( 1, 0.102672);
  AddElementByWeightFraction( 6, 0.56894 );
  AddElementByWeightFraction( 7, 0.035022);
  AddElementByWeightFraction( 8, 0.293366);

  AddMaterial("G4_UREA", 1.323, 0, 72.8, 4);
  AddElementByAtomCount("C" ,  1);
  AddElementByAtomCount("H" ,  4);
  AddElementByAtomCount("N" ,  2);
  AddElementByAtomCount("O" ,  1);

  AddMaterial("G4_WATER_VAPOR", 0.000756182, 0, 71.6, 2, kStateGas);
  AddElementByAtomCount("H" ,  2);
  AddElementByAtomCount("O" ,  1);
  chFormulas[nMaterials-1] = "H_2O-Gas";

  // nNIST = nMaterials;
}

void DPMMaterialBuilder::NistCompoundMaterials3()
{
   // Soft tissue (ICRP - NIST)
    AddMaterial("SoftTissue", 1.00, 0, 69.02, 13);    
    AddElementByWeightFraction(1, 10.4472*perCent);
    AddElementByWeightFraction(6, 23.219*perCent);
    AddElementByWeightFraction(7, 2.488*perCent);
    AddElementByWeightFraction(8, 63.0238*perCent);
    AddElementByWeightFraction(11, 0.113*perCent);
    AddElementByWeightFraction(12, 0.0113*perCent);
    AddElementByWeightFraction(15, 0.113*perCent);
    AddElementByWeightFraction(16, 0.199*perCent);
    AddElementByWeightFraction(17, 0.134*perCent);
    AddElementByWeightFraction(19, 0.199*perCent);
    AddElementByWeightFraction(20, 0.023*perCent);
    AddElementByWeightFraction(26, 0.005*perCent);
    AddElementByWeightFraction(30, 0.003*perCent);

    //  Lung Inhale
    AddMaterial("LungInhale", 0.217, 0, 69.44, 9);    
    AddElementByWeightFraction(1,0.103);
    AddElementByWeightFraction(6,0.105);
    AddElementByWeightFraction(7,0.031);
    AddElementByWeightFraction(8,0.749);
    AddElementByWeightFraction(11,0.002);
    AddElementByWeightFraction(15,0.002);
    AddElementByWeightFraction(16,0.003);
    AddElementByWeightFraction(17,0.002);
    AddElementByWeightFraction(19,0.003);

    // Lung exhale
    AddMaterial("LungExhale", 0.508, 0, 69.44, 9);    
    AddElementByWeightFraction(1,0.103);
    AddElementByWeightFraction(6,0.105);
    AddElementByWeightFraction(7,0.031);
    AddElementByWeightFraction(8,0.749);
    AddElementByWeightFraction(11,0.002);
    AddElementByWeightFraction(15,0.002);
    AddElementByWeightFraction(16,0.003);
    AddElementByWeightFraction(17,0.002);
    AddElementByWeightFraction(19,0.003);

    // Adipose tissue
    AddMaterial("AdiposeTissue", 0.967, 0, 62.23, 7);    
    AddElementByWeightFraction(1,0.114);
    AddElementByWeightFraction(6,0.598);
    AddElementByWeightFraction(7,0.007);
    AddElementByWeightFraction(8,0.278);
    AddElementByWeightFraction(11,0.001);
    AddElementByWeightFraction(16,0.001);
    AddElementByWeightFraction(17,0.001);

    // Brain (ICRP - NIST)
    AddMaterial("BrainTissue", 1.03, 0, 73.3, 13);    
    AddElementByWeightFraction(1, 11.0667*perCent);
    AddElementByWeightFraction(6, 12.542*perCent);
    AddElementByWeightFraction(7, 1.328*perCent);
    AddElementByWeightFraction(8, 73.7723*perCent);
    AddElementByWeightFraction(11, 0.1840*perCent);
    AddElementByWeightFraction(12, 0.015*perCent);
    AddElementByWeightFraction(15, 0.356*perCent);
    AddElementByWeightFraction(16, 0.177*perCent);
    AddElementByWeightFraction(17, 0.236*perCent);
    AddElementByWeightFraction(19, 0.31*perCent);
    AddElementByWeightFraction(20, 0.009*perCent);
    AddElementByWeightFraction(26, 0.005*perCent);
    AddElementByWeightFraction(30, 0.001*perCent); 

   // Breast
    AddMaterial("Breast", 0.990, 0, 70.3, 8);    
    AddElementByWeightFraction(1,0.109);
    AddElementByWeightFraction(6,0.506);
    AddElementByWeightFraction(7,0.023);
    AddElementByWeightFraction(8,0.358);
    AddElementByWeightFraction(11,0.001);
    AddElementByWeightFraction(15,0.001);
    AddElementByWeightFraction(16,0.001);
    AddElementByWeightFraction(17,0.001);

    // Spinal Disc
    AddMaterial("SpinalDisc", 1.10, 0, 74.75, 8);    
    AddElementByWeightFraction(1, 9.60*perCent);
    AddElementByWeightFraction(6, 9.90*perCent);
    AddElementByWeightFraction(7, 2.20*perCent);
    AddElementByWeightFraction(8, 74.40*perCent);
    AddElementByWeightFraction(11, 0.50*perCent);
    AddElementByWeightFraction(15, 2.20*perCent);
    AddElementByWeightFraction(16, 0.90*perCent);
    AddElementByWeightFraction(17, 0.30*perCent);

    // Muscle
    AddMaterial("Muscle", 1.061, 0, 75.3, 9);    
    AddElementByWeightFraction(1,0.102);
    AddElementByWeightFraction(6,0.143);
    AddElementByWeightFraction(7,0.034);
    AddElementByWeightFraction(8,0.710);
    AddElementByWeightFraction(11,0.001);
    AddElementByWeightFraction(15,0.002);
    AddElementByWeightFraction(16,0.003);
    AddElementByWeightFraction(17,0.001);
    AddElementByWeightFraction(19,0.004);

    // Liver
    AddMaterial("Liver", 1.071, 0, 69.02, 9);    
    AddElementByWeightFraction(1,0.102);
    AddElementByWeightFraction(6,0.139);
    AddElementByWeightFraction(7,0.030);
    AddElementByWeightFraction(8,0.716);
    AddElementByWeightFraction(11,0.002);
    AddElementByWeightFraction(15,0.003);
    AddElementByWeightFraction(16,0.003);
    AddElementByWeightFraction(17,0.002);
    AddElementByWeightFraction(19,0.003);

    // Tooth Dentin
    AddMaterial("ToothDentin", 2.14, 0, 80, 10);    
    AddElementByWeightFraction(1, 2.67*perCent);
    AddElementByWeightFraction(6, 12.77*perCent);
    AddElementByWeightFraction(7, 4.27*perCent);
    AddElementByWeightFraction(8, 40.40*perCent);
    AddElementByWeightFraction(11, 0.65*perCent);
    AddElementByWeightFraction(12, 0.59*perCent);
    AddElementByWeightFraction(15, 11.86*perCent);
    AddElementByWeightFraction(17, 0.04*perCent);
    AddElementByWeightFraction(20, 26.74*perCent);
    AddElementByWeightFraction(30, 0.01*perCent);

    // Trabecular Bone
    AddMaterial("TrabecularBone", 1.159, 0, 63.97, 12);    
    AddElementByWeightFraction(1,0.085);
    AddElementByWeightFraction(6,0.404);
    AddElementByWeightFraction(7,0.058);
    AddElementByWeightFraction(8,0.367);
    AddElementByWeightFraction(11,0.001);
    AddElementByWeightFraction(12,0.001);
    AddElementByWeightFraction(15,0.034);
    AddElementByWeightFraction(16,0.002);
    AddElementByWeightFraction(17,0.002);
    AddElementByWeightFraction(19,0.001);
    AddElementByWeightFraction(20,0.044);
    AddElementByWeightFraction(26,0.001);

    // Trabecular bone used in the DICOM Head   
    AddMaterial("TrabecularBone_HEAD", 1.18, 0, 63.97, 12);    
    AddElementByWeightFraction(1, 8.50*perCent);
    AddElementByWeightFraction(6, 40.40*perCent);
    AddElementByWeightFraction(7, 2.80*perCent);
    AddElementByWeightFraction(8, 36.70*perCent);
    AddElementByWeightFraction(11, 0.10*perCent);
    AddElementByWeightFraction(12, 0.10*perCent);
    AddElementByWeightFraction(15, 3.40*perCent);
    AddElementByWeightFraction(16, 0.20*perCent);
    AddElementByWeightFraction(17, 0.20*perCent);
    AddElementByWeightFraction(19, 0.10*perCent);
    AddElementByWeightFraction(20, 7.40*perCent);
    AddElementByWeightFraction(26, 0.10*perCent);

    // Dense Bone
    AddMaterial("DenseBone", 1.575, 0, 91.9, 11);    
    AddElementByWeightFraction(1,0.056);
    AddElementByWeightFraction(6,0.235);
    AddElementByWeightFraction(7,0.050);
    AddElementByWeightFraction(8,0.434);
    AddElementByWeightFraction(11,0.001);
    AddElementByWeightFraction(12,0.001);
    AddElementByWeightFraction(15,0.072);
    AddElementByWeightFraction(16,0.003);
    AddElementByWeightFraction(17,0.001);
    AddElementByWeightFraction(19,0.001);
    AddElementByWeightFraction(20,0.146);

    // Cortical Bone (ICRP - NIST)
    AddMaterial("CorticalBone", 1.85, 0, 106.4, 9);    
    AddElementByWeightFraction(1, 4.7234*perCent);
    AddElementByWeightFraction(6, 14.4330*perCent);
    AddElementByWeightFraction(7, 4.199*perCent);
    AddElementByWeightFraction(8, 44.6096*perCent);
    AddElementByWeightFraction(12, 0.22*perCent);
    AddElementByWeightFraction(15, 10.497*perCent);
    AddElementByWeightFraction(16, 0.315*perCent);
    AddElementByWeightFraction(20, 20.993*perCent);
    AddElementByWeightFraction(30, 0.01*perCent);

    // Tooth enamel 
    AddMaterial("ToothEnamel", 2.89, 0, 80, 10);    
    AddElementByWeightFraction(1, 0.95*perCent);
    AddElementByWeightFraction(6, 1.11*perCent);
    AddElementByWeightFraction(7, 0.23*perCent);
    AddElementByWeightFraction(8,41.66*perCent);
    AddElementByWeightFraction(11, 0.79*perCent);
    AddElementByWeightFraction(12, 0.23*perCent);
    AddElementByWeightFraction(15, 18.71*perCent);
    AddElementByWeightFraction(17, 0.34*perCent);
    AddElementByWeightFraction(20, 35.97*perCent);
    AddElementByWeightFraction(30, 0.02*perCent);

    nNIST = nMaterials;
}
