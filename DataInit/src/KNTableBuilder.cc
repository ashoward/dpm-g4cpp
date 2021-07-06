#include "KNTableBuilder.hh"

#include "AliasTable.hh"

#include <cmath>
#include <cstdio>

KNTableBuilder::KNTableBuilder() {
  fMinPrimaryEnergy                   = 0.001;  // will be set to gamma cut in BuildTables()
  fMaxPrimaryEnergy                   = 21.0;   // will be set in BuildTables()
  fLogMinPrimaryEnergy                = -1.0;
  fInvLogDeltaPrimaryEnergy           = -1.0;
  //
  fNumSamplingPrimaryEnergies = 32;
  fSamplingPrimaryEnergies.resize(fNumSamplingPrimaryEnergies);
  fLSamplingPrimaryEnergies.resize(fNumSamplingPrimaryEnergies);
  //
  fNumSamplingSecondaryEnergies = 64;
  fTables.resize(fNumSamplingPrimaryEnergies, nullptr);
  //
  fAliasSampler = new AliasTable();
}


KNTableBuilder::~KNTableBuilder() {
  for (int ie=0; ie<fNumSamplingPrimaryEnergies; ++ie) {
    if (fTables[ie]) {
      fTables[ie]->Free();
    }
    delete fTables[ie];
  }
  fTables.clear();
}

// here it must be sure that pekin > gcut where gcut is the photon absorption as well
double KNTableBuilder::SampleEnergyTransfer(double pekin, double r1, double r2, double r3) {
  const double kEMC2 = 0.510991;
  // determine primary photon energy lower grid point
  double lpenergy  = std::log(pekin);
  int penergyindx  = (int) ((lpenergy-fLogMinPrimaryEnergy)*fInvLogDeltaPrimaryEnergy);
  penergyindx      = std::min(fNumSamplingPrimaryEnergies-2, penergyindx);
  double ploweenr  = (fLSamplingPrimaryEnergies[penergyindx+1]-lpenergy)*fInvLogDeltaPrimaryEnergy;
  if (r1>ploweenr) {
    ++penergyindx;
  }
  // sample the transformed variable
  LinAlias* theTable =  fTables[penergyindx];
  // sample the transformed variable xi=[\alpha-ln(ep)]/\alpha (where \alpha=ln(1/(1+2\kappa)))
  // that is in [0,1] when ep is in [ep_min=1/(1+2\kappa),ep_max=1] (that limits comes from energy and momentum
  // conservation in case of scattering on free electron at rest).
  // where ep = E_1/E_0 and kappa = E_0/(mc^2)
  double xi = fAliasSampler->SampleLinear(theTable->fXdata, theTable->fYdata,
                                          theTable->fAliasW, theTable->fAliasIndx,
                                          fNumSamplingSecondaryEnergies,r2,r3);
  // transform it back to eps = E_1/E_0
  // \epsion(\xi) = \exp[ \alpha(1-\xi) ] = \exp [\ln(1+2\kappa)(\xi-1)]
  double kappa = pekin/kEMC2;
  return std::exp(std::log(1.+2.*kappa)*(xi-1.)); // eps = E_1/E_0
}

void KNTableBuilder::BuildTables(double gcut, double maxEnergy) {
  const double kEMC2 = 0.510991;
  // build the discrete primary energy grid
  fMinPrimaryEnergy       = gcut;
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
    // add 1 eV to the very first kinetic energy that is equal to gcut
    if (ie==0) {
      primener += 1.0E-6;
    }
    // lower the last kinetic energy with 1 eV
    if (ie==fNumSamplingPrimaryEnergies-1) {
      primener -= 1.0E-6;
    }
    // primary gamma energy in electron rest mas energy units
    double kappa  = primener/kEMC2;
    // allocate data structure to store a single table at this electron energy
    LinAlias* aTable = new LinAlias(fNumSamplingSecondaryEnergies);
    //
    // Build the table
    // note: the transformd variable (xi) is in [0,1] when eps=E_1/E_0 in [eps_min, eps_max]
    // so fill 3 initial values of xi:
    //  -  xi_0 = x_min = 0
    //  -  xi_1 = (x_max-x_min)/2 = 0.5
    //  -  xi_2 = x_max = 1
    // and the corresponding y(i.e.~PDF) values
    aTable->fXdata[0] = 0.0;
    aTable->fXdata[1] = 0.5;
    aTable->fXdata[2] = 1.0;
    aTable->fYdata[0] = ComputeDXSection(aTable->fXdata[0],kappa);
    aTable->fYdata[1] = ComputeDXSection(aTable->fXdata[1],kappa);
    aTable->fYdata[2] = ComputeDXSection(aTable->fXdata[2],kappa);
    int curNumData = 3;
    // expand the data up to numdata points by inserting new points where the
    // linear interpolation error (that will be used at sampling) is the highest
    while (curNumData<aTable->fNumdata) {
      // find the lower index of the bin, where we have the biggest linear interp. error compared to spline
      double maxerr     = 0.0; // value of the current maximum error
      double thexval    = 0.0;
      double theyval    = 0.0;
      int    maxerrindx = 0;   // the lower index of the corresponding bin
      for (int i=0; i<curNumData-1; ++i) {
        double xx    = 0.5*(aTable->fXdata[i]+aTable->fXdata[i+1]);    // mid x point
        double yy    = 0.5*(aTable->fYdata[i]+aTable->fYdata[i+1]);    // lin. interpolated pdf value at the mid point
        double val   = ComputeDXSection(xx,kappa); // real pdf value at the mid point
        double err   = std::abs(1.-(yy/val));
        if (err>maxerr) {
          maxerr     = err;
          maxerrindx = i;
          thexval    = xx;
          theyval    = val;
        }
      }
      // extend x,y data by puting a new real value at the mid point of the highest error bin
      // first shift all values to the right
      for (int j=curNumData; j>maxerrindx+1; --j) {
        aTable->fXdata[j] = aTable->fXdata[j-1];
        aTable->fYdata[j] = aTable->fYdata[j-1];
      }
      // fill x mid point
      aTable->fXdata[maxerrindx+1] = thexval;
      aTable->fYdata[maxerrindx+1] = theyval;
      // increase number of data
      ++curNumData;
    } // end while
    // init the alias data structure for this table
    fAliasSampler->PreparLinearTable(aTable->fXdata, aTable->fYdata,
                                     aTable->fAliasW, aTable->fAliasIndx,
                                     fNumSamplingSecondaryEnergies);//fAliasData[ialias]->fNumdata);
    fTables[ie] = aTable;
  }
}


void KNTableBuilder::Write(const std::string& dirname) {
  char name[512];
  FILE* f = nullptr;
  sprintf(name, "%s/compton_KNDtrData.dat", dirname.c_str());
  f = fopen(name, "w");
  fprintf(f, "# Tables for rejection free sampling of reduced energy transfer\n");
  fprintf(f, "# in Compton scattering according to the Klein-Nishina DCS.\n");
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


// pdf in the transformed variable
double KNTableBuilder::ComputeDXSection(double xi, double kappa) {
  double inv2Kappa  = 1./(1.+2.*kappa);
  double linv2Kappa = std::log(inv2Kappa);
  double eps        = std::exp(linv2Kappa*(1.-xi));
  double invEps     = 1./eps;
  double beta       = (1.-invEps)/kappa;

  return eps*(eps+invEps)*(1.+eps*beta*(2.+beta)/(1.+eps*eps));
}
