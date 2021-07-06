#ifndef KNTableBuilder_HH
#define KNTableBuilder_HH

// M. Novak: 2021
//
// A utility class to build sampling tables that can be used at run time to
// obtain values of energy transferred to the secondary electron in case of
// compton scattering according to the Klein-Nishina DCS without rejection.


#include <vector>
#include <string>

class AliasTable;

class KNTableBuilder {
public:
  KNTableBuilder();
 ~KNTableBuilder();

 // tables are built between the photon absorption energy (gamma-cut) and the
 // maximum energy used in the simulation
 void   BuildTables(double gcut, double maxEnergy);
 double SampleEnergyTransfer(double pekin, double r1, double r2, double r3);

 // write the data into file
 void   Write(const std::string& dirname);

private:
 double ComputeDXSection(double xi, double kappa);

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

  // The minimum discrete primary photon energy at which Compton is possible (will be set to gamma cut)
  double                fMinPrimaryEnergy;
  // The maximum discrete primary photon energy used in the simulation
  double                fMaxPrimaryEnergy;
  // The logarithm of the minimum primary energy
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
