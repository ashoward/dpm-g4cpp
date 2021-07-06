#ifndef MollerTableBuilder_HH
#define MollerTableBuilder_HH

// M. Novak: 2021
//
// A utility class to build sampling tables that can be used at run time to
// obtain values of energy transferred to the secondary electron in Moller
// interaction without rejection.


#include <vector>
#include <string>


class AliasTable;

class MollerTableBuilder {
public:
  MollerTableBuilder();
 ~MollerTableBuilder();

 // tables are built between 2x the secondary electron production energy and the
 // maximum energy used in the simulation
 void   BuildTables(double ecut, double maxEnergy);
 double SampleEnergyTransfer(double eekin, double ecut, double r1, double r2, double r3);

 // write the data into file
 void   Write(const std::string& dirname);

private:
 double ComputeMollerPDF(double xi, double prodcutenergy, double particleekin);

private:
  // Sampling table: at one discrete electron energy
  struct LinAlias{
    LinAlias(int size) : fNumdata(size) {
      fXdata     = new double[size]();
      fYdata     = new double[size]();
      fAliasW    = new double[size]();
      fAliasIndx = new int[size]();
    }
    void Free () {
      delete [] fXdata;
      delete [] fYdata;
      delete [] fAliasW;
      delete [] fAliasIndx;
    }

    // Number of data points i.e. fNumSamplingPhotEnergies
    int     fNumdata;
    // Reduced photon energy related transformed variable values.
    double *fXdata;
    // The probability density function values (not necessarily normalised) over fXData
    double *fYdata;
    // The alias probabilities (not necessarily normalised)  over fXData
    double *fAliasW;
    // The alias indices  over fXData
    int    *fAliasIndx;
  };

  // logarithm of twice e- production threshold
  double fLogTwoElectronProdCut;

  // The minimum discrete primary electron energy at which Moller is possible (will be set to 2 x e- cut)
  double                fMinPrimaryEnergy;
  // The maximum discrete primary electron energy used in the simulation
  double                fMaxPrimaryEnergy;
  // The logarithm of the minimum electron energy
  double                fLogMinPrimaryEnergy;
  // Inverse delta log kinetic energy
  double                fInvLogDeltaPrimaryEnergy;
  // The number of discrete primary energy points at which sampling tables will be built between the abve two.
  int                   fNumSamplingPrimaryEnergies;
  // The discrete energy grid for the sampling tables
  std::vector<double>   fSamplingPrimaryEnergies; // [fNumSamplingPrimaryElectronEnergies]
  // The logarithm of the above grid values
  std::vector<double>   fLSamplingPrimaryEnergies; // [fNumSamplingPrimaryElectronEnergies]

  // The number of secondary electron energy related transormed variables over tables are built.
  int                   fNumSamplingSecondaryEnergies;

  //
  // Set of sampling tables
  std::vector<LinAlias*> fTables; // [fNumSamplingPrimaryElectronEnergie]

  AliasTable*            fAliasSampler;
};

#endif
