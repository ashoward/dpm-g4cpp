
#include "SBTableBuilder.hh"

#include "Spline.hh"
#include "AliasTable.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include <iostream>
#include <cmath>


SBTableBuilder::SBTableBuilder(const std::string& datadir) {
  fDataFileLocation                  = datadir;
  fLoadedDCSNumElectronEnergies      = 113;
  fLoadedDCSNumReducedPhotonEnergies = 32;
  fLoadedDCSElectronEnergyGrid.resize(fLoadedDCSNumElectronEnergies);
  fLoadedDCSReducedPhotonEnergyGrid.resize(fLoadedDCSNumReducedPhotonEnergies);
  fLoadedDCSForElements.resize(100, std::vector<double>( fLoadedDCSNumElectronEnergies*fLoadedDCSNumReducedPhotonEnergies) );
  //
  fMinElecEnergy                     = 0.001; // will be set to the gamma cut in BuildTables()
  fMaxElecEnergy                     =  21.0; // will be set to the max-energy in BuildTables()
  fElEnLMin                          =  -1.0;
  fElEnILDelta                       =  -1.0;
  fNumSamplingElecEnergies           =   128;
  fSamplingElecEnergies.resize(fNumSamplingElecEnergies);
  fLSamplingElecEnergies.resize(fNumSamplingElecEnergies);
  //
  fNumSamplingPhotEnergies           = 64;
  //
  fAliasSampler = new AliasTable();
  //
}

SBTableBuilder::~SBTableBuilder() {
  for (std::size_t imat=0; imat<fTablesPerMaterial.size(); ++imat) {
    for (int ie=0; ie<fNumSamplingElecEnergies; ++ie) {
      fTablesPerMaterial[imat].fTheTables[ie]->Free();
    }
    fTablesPerMaterial[imat].fTheTables.clear();
  }
  fTablesPerMaterial.clear();
  delete fAliasSampler;
}

// here it must be sure that eekin > k_c !!!
double SBTableBuilder::SampleEnergyTransfer(int imat, double eekin, double gcut, double r1, double r2, double r3) {
  // determine primary electron energy lower grid point
  double leenergy  = std::log(eekin);
  int eenergyindx  = (int) ((leenergy-fElEnLMin)*fElEnILDelta);
  eenergyindx      = std::min(fNumSamplingElecEnergies-2, eenergyindx);
  double ploweenr  = (fLSamplingElecEnergies[eenergyindx+1]-leenergy)*fElEnILDelta;
  if (r1>ploweenr) {
    ++eenergyindx;
  }
  // sample the transformed variable
  LinAlias* theTable =  fTablesPerMaterial[imat].fTheTables[eenergyindx];
  double xi = fAliasSampler->SampleLinear(theTable->fXdata, theTable->fYdata,
                                          theTable->fAliasW, theTable->fAliasIndx,
                                          fNumSamplingPhotEnergies,r2,r3);
  // transform it back to kappa then to gamma energy
  double kappac = gcut/eekin;
  double kappa  = kappac*std::exp(-xi*std::log(kappac));
  return kappa*eekin;
}

void SBTableBuilder::BuildTables(int numMaterials, double gcut, double maxEnergy) {
  // load numerical sclalled DCS data
  LoadDCSData();
  //
  int theNumMaterials       = numMaterials;
  double theGammaCutE       = gcut;
  double theMaxEnergy       = maxEnergy;
  // build the discrete kinetic energy grid for the smapling tables between k_c and E_max
  fMinElecEnergy            = theGammaCutE;
  fMaxElecEnergy            = theMaxEnergy;
  double theElEnLMin        = std::log(fMinElecEnergy);
  fElEnLMin                 = theElEnLMin;
  double delta              = std::log(fMaxElecEnergy/fMinElecEnergy)/(fNumSamplingElecEnergies-1.0);
  double theElEnILDelta     = 1.0/delta;
  fElEnILDelta              = theElEnILDelta;
  fSamplingElecEnergies[0]  = fMinElecEnergy;
  fLSamplingElecEnergies[0] = theElEnLMin;
  fSamplingElecEnergies[fNumSamplingElecEnergies-1]  = fMaxElecEnergy;
  fLSamplingElecEnergies[fNumSamplingElecEnergies-1] = std::log(fMaxElecEnergy);
  for (int i=1; i<fNumSamplingElecEnergies-1; ++i) {
    fLSamplingElecEnergies[i] = theElEnLMin+i*delta;
    fSamplingElecEnergies[i]  = std::exp(theElEnLMin+i*delta);
  }
  //
  //
  // allocate container to store a set of tables for the individual materials
  fTablesPerMaterial.resize(theNumMaterials);
  for (int imat=0; imat<theNumMaterials; ++imat) {
    // the first material is the 'Universe' so skipp that one
    const G4Material* theMaterial = G4NistManager::Instance()->GetMaterial(imat+1);
    // allocate container to store a set of tables for this material at the individual
    // discrete electron energies
    fTablesPerMaterial[imat].fTheTables.resize(fNumSamplingElecEnergies);
    fTablesPerMaterial[imat].fMaterialName = theMaterial->GetName();
    for (int ie=0; ie<fNumSamplingElecEnergies; ++ie) {
      // allocate data structure to store a single table at this electron energy
      LinAlias* aTable = new LinAlias(fNumSamplingPhotEnergies);
      double elEnergy = fSamplingElecEnergies[ie];
      // add 1 eV to the very first kinetic energy that is equal to k_c
      if (ie==0) {
        elEnergy += 1.0E-6;
      }
      // lower the last kinetic energy with 1 eV
      if (ie==fNumSamplingElecEnergies-1) {
        elEnergy -= 1.0E-6;
      }
      BuildOneLinAlias(theMaterial, theGammaCutE, aTable, elEnergy);
      fTablesPerMaterial[imat].fTheTables[ie] = aTable;
    }
  }
}

// builds one sampling tables for the given material at the given discrete e-
// kinetic energy `E` ( k_c <= E <= E_max)
void SBTableBuilder::BuildOneLinAlias(const G4Material *mat, double gcut, LinAlias* theTable, double eener) {
//  std::cout << " Building Table for material: " << mat->GetName() << " E = " << eener << std::endl;
  // we will need the element composition of this material
  const double* elemDensity       = mat->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* elemVect = mat->GetElementVector();
  const int numElems              = elemVect->size();
  // find the discrete e- energy index `i` in the available (loaded) e- energy grid
  // grid such that E_i <= E < E_i+1
  // NOTE: E is such that E_0 < E < E_max for sure with E_0 = 1 keV and E_max = 1 GeV
  int eenerindx = 0;
  for (; eenerindx<fLoadedDCSNumElectronEnergies && eener>=fLoadedDCSElectronEnergyGrid[eenerindx]; ++eenerindx) {}
  --eenerindx; // to make it lower index
  // using the scalled DCS-s at E_i and E_i+1, interpolate the DCS to E (using linear on log scale)
  const int numKappas = fLoadedDCSNumReducedPhotonEnergies;
  std::vector<double> theDCS(numKappas);
  std::vector<double> logReducedPhotonEnergyGrid(numKappas);
  // ln(x)-ln(x1)
  double dum0 = std::log(eener/fLoadedDCSElectronEnergyGrid[eenerindx]);
  // ln(x2)-ln(x1)
  double dum1 = std::log(fLoadedDCSElectronEnergyGrid[eenerindx+1]/fLoadedDCSElectronEnergyGrid[eenerindx]);
  for (int irpener=0; irpener<numKappas; ++irpener) {
     // transform the \kappa values to \phi(\kappa) = \ln [ 1.0 - \kappa]
     // NOTE: \kappa_max = 1.0 so add 1.0E-12 here to avoid probelms
     logReducedPhotonEnergyGrid[irpener] = std::log(1.0-fLoadedDCSReducedPhotonEnergyGrid[irpener]+1.0e-12);
     int indxdcsl = eenerindx*numKappas + irpener;
     int indxdcsh = indxdcsl + numKappas; //(eenerindx+1)*numKappas + irpener;
     theDCS[irpener] = 0.0;
     for (int ielem=0; ielem<numElems; ++ielem) {
       const G4Element *elem = (*elemVect)[ielem];
       const double      zet = elem->GetZ();
       const int        izet = std::lrint(zet);
       // ln(y2) -ln(y1)
       double dum2 = std::log(fLoadedDCSForElements[izet][indxdcsh]/fLoadedDCSForElements[izet][indxdcsl]);
       double dcs  = dum2/dum1*dum0+std::log(fLoadedDCSForElements[izet][indxdcsl]); //this is ln(dcs)
       dcs = std::exp(dcs);
       dcs *= elemDensity[ielem]*zet*zet;
       theDCS[irpener] += dcs;
     }
   }
   //
   // set up a spline on this interpolated scalled DCS
   Spline     *sp = new Spline();
   sp->SetUpSpline(logReducedPhotonEnergyGrid.data(), theDCS.data(), fLoadedDCSNumReducedPhotonEnergies);//true, true);
   //
   // compute \kappa_c = k_c/E (that is for sure < 1.0) and find \kappa_i <= \kappa_c < \kappa_i+1
   const double theKappaC    = gcut/eener;
   const double theLogKappaC = std::log(theKappaC);
   int kappaindx = 0;
   for (; kappaindx<fLoadedDCSNumReducedPhotonEnergies-1 && theKappaC>=fLoadedDCSReducedPhotonEnergyGrid[kappaindx]; ++kappaindx) {}
   --kappaindx;
   //
   // NOTE: that the sampling table is built in \xi(\kappa) = 1 - \ln(\kappa)/ln(\kappac)
   //       and the discrete loaded kappa grid was transformed to \phi(\kappa) = \ln [1.0-\kappa+1.0E-12]
   //       a. \kappa_max = 1 so add 1.0E-12 here to avoid problem
   //       b. \kappa_c = k_c/E < 1 is guarantied since we add 1 eV to the very first E = k_c
   // 1. fill in the very first value that corresponds to \kappa_c :
   //     - this is \xi(\kappa = \kappa_c) = 0 and \phi(\kappa = \kappa_c) = \ln [1.0E-6+1.0E-12]
   int numdata = 1;
   double xi  = 0.0;
   double phi = std::log(1.0-theKappaC+1.0e-12);
   theTable->fXdata[0] = 0.0;
   theTable->fYdata[0] = sp->GetValueAt(phi);
   // 2. then add all \xi(\kappa_i) values such tha \kappa_c < \kappa_i < 1
   for (int k=kappaindx+1; k<fLoadedDCSNumReducedPhotonEnergies-1; ++k) {
     double theKappa = fLoadedDCSReducedPhotonEnergyGrid[k];
     xi = 1.0 - std::log(theKappa)/theLogKappaC;
     theTable->fXdata[numdata] = xi;
     theTable->fYdata[numdata] = theDCS[k];
     ++numdata;
   }
   // 3. then insert the very last xi(\kappa=1) = 1 value.
   theTable->fXdata[numdata] = 1.0;
   theTable->fYdata[numdata] = theDCS[numKappas-1];
   ++numdata;
   //
   // Then expand the muber of \xi values up to fNumSamplingPhotEnergies by
   // inserting a new \xi value between \xi_i and \xi_i+1 that gives the largest
   // linear interpolation error
   while (numdata<fNumSamplingPhotEnergies) {
     // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
     double maxerr     = 0.0; // value of the current maximum error
     double thexval    = 0.0;
     double theyval    = 0.0;
     int    maxerrindx = 0;   // the lower index of the corresponding bin

     for (int i=0; i<numdata-1; ++i) {
       double xx = 0.5*(theTable->fXdata[i]+theTable->fXdata[i+1]); // mid point
       double yy = 0.5*(theTable->fYdata[i]+theTable->fYdata[i+1]); // lin func val at the mid point
       // compute \kappa(\xi=xx) = \kappa_c \exp{-xi*\ln{\kappa_c}} and \phi(\kappa)
       double thekappa = theKappaC*std::exp(-xx*theLogKappaC);
       double   thePhi = std::log(1.0-thekappa+1.0e-12);
       double    spval = sp->GetValueAt(thePhi); // spline intp. val. at mid point
       double    err   = std::fabs(yy-spval);
       if (err>maxerr) {
         maxerr     = err;
         maxerrindx = i;
         thexval    = xx;
         theyval    = spval;
       }
     }
     // extend x,y data by puting a spline interp.ted value at the mid point of
     // the highest error bin: first shift all values to the right
     for (int j=numdata; j>maxerrindx+1; --j) {
       theTable->fXdata[j] = theTable->fXdata[j-1];
       theTable->fYdata[j] = theTable->fYdata[j-1];
     }
     // fill x mid point
     theTable->fXdata[maxerrindx+1] = thexval;
     theTable->fYdata[maxerrindx+1] = theyval;
     // increase number of data
     ++numdata;
   }
   // set up a linear alias smapler on this data
   fAliasSampler->PreparLinearTable(theTable->fXdata, theTable->fYdata,
                                    theTable->fAliasW, theTable->fAliasIndx,
                                    fNumSamplingPhotEnergies);
   //
   delete sp;
}


void SBTableBuilder::LoadDCSData() {
  const double kMBarnToMM2 = 1.0E-25;
  char baseFilename[512];
  sprintf(baseFilename,"%s/NIST_BREM/nist_brems_",fDataFileLocation.c_str());
  // load the discrete kinetic energy and reduced photon energy grids
  FILE *f = nullptr;
  char filename[512];
  sprintf(filename,"%sgrid",baseFilename);
  f = fopen(filename,"r");
  if (f==nullptr) {
    std::cerr << "******   ERROR in SBTableBuilder::LoadDCSData() \n"
              << "         " << filename << " could not be found!\n"
              << std::endl;
    exit(1);
  }
  // skipp the first line 3 number
  int idum;
  fscanf(f,"%d%d%d",&idum,&idum,&idum);
  // read the fLoadedDCSNumElectronEnergies discrete kinetic energy values [MeV]
  for (int i=0; i<fLoadedDCSNumElectronEnergies; ++i) {
    fscanf(f,"%lg",&(fLoadedDCSElectronEnergyGrid[i]));
  }
  // read the fLoadedDCSNumReducedPhotonEnergies reduced photon energy values
  for (int i=0; i<fLoadedDCSNumReducedPhotonEnergies; ++i) {
    fscanf(f,"%lg",&(fLoadedDCSReducedPhotonEnergyGrid[i]));
  }
  fclose(f);
  //
  // now go and load the scalled SB brem. DCS for each of the 99 element
  const int numDataPerElement = fLoadedDCSNumElectronEnergies*fLoadedDCSNumReducedPhotonEnergies;
  const int numElem = fLoadedDCSForElements.size();
  for (int ie=1; ie<numElem; ++ie) {
    // std::cout << " loading data for element Z = " << ie << std::endl;
    sprintf(filename,"%s%d",baseFilename,ie);
    f = fopen(filename,"r");
    if (f==nullptr) {
      std::cerr << "******   ERROR in SBTableBuilder::LoadDCSData() \n"
                << "         " << filename << " could not be found!\n"
                << std::endl;
      exit(1);
    }
    for (int i=0; i<numDataPerElement; ++i) {
      fscanf(f,"%lg",&(fLoadedDCSForElements[ie][i]));
      // change to [mm2] units because they are stored in [mbarn] units
      fLoadedDCSForElements[ie][i] *= kMBarnToMM2;

    }
    fclose(f);
  }
}


void SBTableBuilder::Write(const std::string& dirname) {
  char name[512];
  FILE* f = nullptr;
  sprintf(name, "%s/brem_SBDtrData.dat", dirname.c_str());
  f = fopen(name, "w");
  fprintf(f, "# Tables for rejection free sampling of the emitted photon energy\n");
  fprintf(f, "# in electron bremsstrahlung interaction according to the Seltzer-Berger DCS.\n");
  fprintf(f, "# The size of the discrete primary (np), secondary (ns) energies and #materials first \n");
  fprintf(f, "# followed by the discrete primary energy grid and for each material (np) tables with \n");
  fprintf(f, "# (4 x ns) sampling table data each. \n");
  // size of the discrete primary and secondary energy grids
  fprintf(f, "%d %d %d\n", fNumSamplingElecEnergies, fNumSamplingPhotEnergies, (int)(fTablesPerMaterial.size()));
  // first we write out the fNumSamplingElecEnergies size discrete primary eletron
  // energy grid (this is common for each material, with a minimum of gamma-cut)
  for (int ie=0; ie<fNumSamplingElecEnergies; ++ie) {
    fprintf(f,"%22.14E ", fSamplingElecEnergies[ie]);
    if ( (ie+1)%6==0 || ie==fNumSamplingElecEnergies-1) {
      fprintf(f, "\n");
    }
  }
  // then for each material:
  //  - first the material name in the format:
  //      #
  //      #  MATERIAL NAME: name
  //      #
  //  - then the `fNumSamplingElecEnergies` sampling tables (Alias table with
  //    linear pdf approximation) with a size of `4 x fNumSamplingPhotEnergies`
  //    each
  for (std::size_t im=0; im<fTablesPerMaterial.size(); ++im) {
    // write the material name into comment
    fprintf(f, "#\n# Material name: %s\n#\n", fTablesPerMaterial[im].fMaterialName.c_str());
    // then the `fNumSamplingElecEnergies` tables
    for (int ie=0; ie<fNumSamplingElecEnergies; ++ie) {
      // each having the 4 x fNumSamplingPhotEnergies data
      for (int is=0; is<fNumSamplingPhotEnergies; ++is) {
        fprintf(f, "%22.14E  %22.14E  %22.14E  %d ",
                   fTablesPerMaterial[im].fTheTables[ie]->fXdata[is],
                   fTablesPerMaterial[im].fTheTables[ie]->fYdata[is],
                   fTablesPerMaterial[im].fTheTables[ie]->fAliasW[is],
                   fTablesPerMaterial[im].fTheTables[ie]->fAliasIndx[is]
                 );
        if ( (is+1)%2==0 || is==fNumSamplingPhotEnergies-1) {
          fprintf(f, "\n");
        }
      }
    }
  }
  fclose(f);
}
