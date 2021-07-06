#include "MollerTableBuilder.hh"

#include "AliasTable.hh"

#include <cmath>
#include <cstdio>

MollerTableBuilder::MollerTableBuilder() {
  fMinPrimaryEnergy                   = 0.001;  // will be set to 2x e- cut in BuildTables()
  fMaxPrimaryEnergy                   = 21.0;   // will be set in BuildTables()
  fLogMinPrimaryEnergy                = -1.0;
  fInvLogDeltaPrimaryEnergy           = -1.0;
  //
  fLogTwoElectronProdCut             = -1.0;
  //
  fNumSamplingPrimaryEnergies = 128;
  fSamplingPrimaryEnergies.resize(fNumSamplingPrimaryEnergies);
  fLSamplingPrimaryEnergies.resize(fNumSamplingPrimaryEnergies);
  //
  fNumSamplingSecondaryEnergies = 64;
  fTables.resize(fNumSamplingPrimaryEnergies, nullptr);
  //
  fAliasSampler = new AliasTable();
}


MollerTableBuilder::~MollerTableBuilder() {
  for (int ie=0; ie<fNumSamplingPrimaryEnergies; ++ie) {
    if (fTables[ie]) {
      fTables[ie]->Free();
    }
    delete fTables[ie];
  }
  fTables.clear();
}

// here it must be sure that eekin > 2*ecut
double MollerTableBuilder::SampleEnergyTransfer(double eekin, double ecut, double r1, double r2, double r3) {
  // determine primary electron energy lower grid point
  double leenergy  = std::log(eekin);
  int eenergyindx  = (int) ((leenergy-fLogMinPrimaryEnergy)*fInvLogDeltaPrimaryEnergy);
  eenergyindx      = std::min(fNumSamplingPrimaryEnergies-2, eenergyindx);
  double ploweenr  = (fLSamplingPrimaryEnergies[eenergyindx+1]-leenergy)*fInvLogDeltaPrimaryEnergy;
  if (r1>ploweenr) {
    ++eenergyindx;
  }
  // sample the transformed variable
  LinAlias* theTable =  fTables[eenergyindx];
  // sample the transformed variable xi=[kappa-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)]
  // where kappa = ln(eps) with eps = T/T_0
  // so xi= [ln(T/T_0)-ln(T_cut/T_0)]/[ln(T_cut/T_0)-ln(T_max/T_0)] that is in [0,1]
  double xi = fAliasSampler->SampleLinear(theTable->fXdata, theTable->fYdata,
                                          theTable->fAliasW, theTable->fAliasIndx,
                                          fNumSamplingSecondaryEnergies,r2,r3);
  // transform it back
//  const double kLogHalf = -0.69314718055994;
//  double dum0 = eekin/ecut;
//  double dum1 = std::log(dum0);
//  if (fIsElectron)
//    dum1 += kLogHalf;
  const double dum1 = leenergy - fLogTwoElectronProdCut;
  // return with the sampled kinetic energy transfered to the electron
  return std::exp(xi*dum1)*ecut;
}

void MollerTableBuilder::BuildTables(double ecut, double maxEnergy) {
  // build the discrete primary energy grid
  fLogTwoElectronProdCut = std::log(2.0*ecut);
  fMinPrimaryEnergy       = 2.0*ecut;
  fMaxPrimaryEnergy       = maxEnergy;
  fLogMinPrimaryEnergy    = std::log(fMinPrimaryEnergy);
  double delta   = std::log(fMaxPrimaryEnergy/fMinPrimaryEnergy)/(fNumSamplingPrimaryEnergies-1.0);
  fInvLogDeltaPrimaryEnergy = 1.0/delta;
  fSamplingPrimaryEnergies[0]  = fMinPrimaryEnergy;
  fLSamplingPrimaryEnergies[0] = fLogMinPrimaryEnergy;
  fSamplingPrimaryEnergies[fNumSamplingPrimaryEnergies-1]  = fMaxPrimaryEnergy;
  fLSamplingPrimaryEnergies[fNumSamplingPrimaryEnergies-1] = std::log(fMaxPrimaryEnergy);
  for (int i=1; i<fNumSamplingPrimaryEnergies-1; ++i){
    fLSamplingPrimaryEnergies[i] = fLogMinPrimaryEnergy+i*delta;
    fSamplingPrimaryEnergies[i]  = std::exp(fLogMinPrimaryEnergy+i*delta);
  }
  //
  // go over the discrete primary electron energy grid
  for (int ie=0; ie<fNumSamplingPrimaryEnergies; ++ie) {
    double primener = fSamplingPrimaryEnergies[ie];
    // add 1 eV to the very first kinetic energy that is equal to 2xecut
    if (ie==0) {
      primener += 1.0E-6;
    }
    // lower the last kinetic energy with 1 eV
    if (ie==fNumSamplingPrimaryEnergies-1) {
      primener -= 1.0E-6;
    }
    // minimum kinetic energy transferable to the scattered e-
    // double tmax = 0.5*primener;
    // allocate data structure to store a single table at this electron energy
    LinAlias* aTable = new LinAlias(fNumSamplingSecondaryEnergies);
    // fill the x-scale i.e. the transformed xi variable values that are \in [0,1]
    // and the corresponding scattered e- distribution
    double adum = 1.0/(fNumSamplingSecondaryEnergies-1.0);
    for (int i=0; i<fNumSamplingSecondaryEnergies; ++i) {
      double xi = i*adum;
      if (i==0) {
        xi = 0.0;
      } else if (i==fNumSamplingSecondaryEnergies-1) {
        xi = 1.0;
      }
      aTable->fXdata[i] = xi;
      aTable->fYdata[i] = ComputeMollerPDF(xi, ecut, primener);
    }
    // init the alias data structure for this table
    fAliasSampler->PreparLinearTable(aTable->fXdata, aTable->fYdata,
                                     aTable->fAliasW, aTable->fAliasIndx,
                                     fNumSamplingSecondaryEnergies);//fAliasData[ialias]->fNumdata);
    fTables[ie] = aTable;
  }
}


void MollerTableBuilder::Write(const std::string& dirname) {
  char name[512];
  FILE* f = nullptr;
  sprintf(name, "%s/ioni_MollerDtrData.dat", dirname.c_str());
  f = fopen(name, "w");
  fprintf(f, "# Tables for rejection free sampling of energy transferred\n");
  fprintf(f, "# to the secondary electron in Moller interaction.\n");
  fprintf(f, "# The size of the discrete primary (np) and secondary (ns) energies first \n");
  fprintf(f, "# then the primary energy and the corresponding (4 x ns) sampling table data \n");
  fprintf(f, "# each of the (np) discrete primary energies. \n");
  // size of the discrete primary and secondary energy grids
  fprintf(f, "%d %d\n", fNumSamplingPrimaryEnergies, fNumSamplingSecondaryEnergies);
  // there is a sampling table (Alias table with linear pdf approximation) built
  // with fNumSamplingSecondaryEnergies discrete data at each of the
  // fNumSamplingPrimaryEnergies discrete primary kinetic energy points.
  // At each discrete primary energies, first the kinetic energy value then the
  // sampling table is witten
  for (int ie=0; ie<fNumSamplingPrimaryEnergies; ++ie) {
    fprintf(f,"%22.14E\n", fSamplingPrimaryEnergies[ie]);
    for (int is=0; is<fNumSamplingSecondaryEnergies; ++is) {
      fprintf(f, "%22.14E  %22.14E  %22.14E  %d ",
                 fTables[ie]->fXdata[is],  fTables[ie]->fYdata[is],
                 fTables[ie]->fAliasW[is], fTables[ie]->fAliasIndx[is]);
      if ( (is+1)%2==0 || is==fNumSamplingSecondaryEnergies-1) {
        fprintf(f, "\n");
      }
    }
  }
  fclose(f);
}


double MollerTableBuilder::ComputeMollerPDF(double xi, double prodcutenergy, double particleekin) {
  const double kEMC2 = 0.510991;
  double tau         = particleekin/kEMC2;      // i.e. E_kin in (mc^2) units
  double gamma       = tau + 1.0;              // = E_t/(mc^2)  i.e. E_t in (mc^2) units
  double gamma2      = gamma*gamma;            // \gamma^2 = [E_t/(mc^2)]^2
  //double beta2     = tau*(tau+2.0)/gamma2; // \beta2 i.e. [P_t/E_t]^2
  double C1          = (gamma-1.0)/gamma;
  C1 *=C1;
  double C2          = (2.0*gamma-1.0)/gamma2;

  double dum0        = prodcutenergy/particleekin;
  double dum1        = std::log(0.5/dum0);
  double a           = std::exp(xi*dum1)*dum0; // this is eps =  exp(xi*ln(0.5*T_0/T_cut))*T_cut/T_0
  double b           = 1.0-a;                  // eps'

  return ((1.0/a-C2)+a*C1+a/b*(1.0/b-C2)) *dum0; // xdum0 is just scaling; this is the shape
}
