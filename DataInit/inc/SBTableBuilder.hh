
#ifndef SBTableBuilder_HH
#define SBTableBuilder_HH

// M. Novak: 2021
//
// A utility class to build sampling tables that can be used at run time to
// obtain energy values emitted in a form of a bremsstrahlung photon according
// to the Seltzer-Berger (SB) scalled numerical DCS without rejection.

#include <vector>
#include <string>

class G4Material;
class AliasTable;

class SBTableBuilder {
public:
  SBTableBuilder(const std::string& datadir);
 ~SBTableBuilder();

  void   BuildTables(int numMaterials, double gcut, double maxEnergy);
  double SampleEnergyTransfer(int imat, double eekin, double gcut, double r1, double r2, double r3);

  // writing out the sampling table data per material into a file
  void   Write(const std::string& dirname);

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

  struct TableForAMaterial {
    std::string            fMaterialName;
    std::vector<LinAlias*> fTheTables; // [fNumSamplingElecEnergies]
  };

  // Loads the numerical SB DCS data for all the elements
  void LoadDCSData();
  void BuildOneLinAlias(const G4Material *mat, double gcut, LinAlias* theTable, double eener);

// data memebrs
private:
  std::string          fDataFileLocation;

  // Number of electron energies that numerical DCS are available in the file.
  int                  fLoadedDCSNumElectronEnergies;
  // Number of reduced photon energies at each of the above electron energies.
  int                  fLoadedDCSNumReducedPhotonEnergies;
  // The electron energy grid.
  std::vector<double>  fLoadedDCSElectronEnergyGrid;
  // The reduced photon energy grid.
  std::vector<double>  fLoadedDCSReducedPhotonEnergyGrid;
  // The numerical DCS data per elements. Each of the sub-vectors are a 113*32 size vector
  // holding the DCS data for a given element. There are 99 elements
  // pointer referes to a double array with size of 113*32
  std::vector< std::vector<double> > fLoadedDCSForElements;

  // The minimum discrete electron energy at which brem. is possible (will be set to the gamma cut)
  double                fMinElecEnergy;
  // The maximum discrete electron energy used in the simulation
  double                fMaxElecEnergy;
  // The logarithm of the minimum electron energy
  double                fElEnLMin;
  // Inverse delta log kinetic energy
  double                fElEnILDelta;
  // The number of discrete energy points at which sampling tables will be built between the abve two.
  int                   fNumSamplingElecEnergies;
  // The discrete energy grid for the sampling tables
  std::vector<double>   fSamplingElecEnergies;
  // The logarithm of the above grid values
  std::vector<double>   fLSamplingElecEnergies;

  // The number of reduced photon energy related transormed variables over tables are built.
  int                   fNumSamplingPhotEnergies;

  //
  // Set of sampling tables for each material
  std::vector<TableForAMaterial> fTablesPerMaterial; // [#materials]

  AliasTable*           fAliasSampler;
};

#endif
