#ifndef SimKNTables_HH
#define SimKNTables_HH

//
// M. Novak: 2021
//
// Sampling tables for generating (reduced) energy transfer values in Compton
// scattering according to the Klein-Nishina DCS without rejection.
//
// It is assumed, that the required data has been generated and written into
// the `compton_KNDtrData.dat` file previously by using the `KNTableBuilder`.
// This file contains a discrete primary energy grid and a sampling table
// (`SimLinAliasData`) at each of the discrete primary energies. All the data
// can be loaded by using the provided `LoadData(dataDir)` method then the
// `SampleEnergyTransfer()` can be used to generate (reduced) energy transfer
// values in Compton scattering according to the Klein-Nishina DCS at the given
// photon energy without rejection.


#include <vector>
#include <string>

class SimLinAliasData;


class SimKNTables {

public:

  // CTR and DTR
  SimKNTables();
 ~SimKNTables() { CleanTables(); }

  // loads the previously prepared data from the appropriate file locating in
  // the directory specified by the `dataDir` input argument
  void   LoadData(const std::string& dataDir, int verbose);

  // samples (reduced) energy transfer in Compton scattering according to the
  // Klein-Nishina DCS at the given `egamma` primary gamma energy. The three
  // additional input arguments are uniformly random values on [0,1].
  //
  // NOTE: it is assumed that: gamma-cut < egamma < E_max
  double SampleEnergyTransfer(double egamma, double rndm1, double rndm2, double rndm3);



private:

  void   CleanTables();


private:

  // The size of the individual sampling tables
  int                   fSamplingTableSize;
  // The number of discrete primary energies at which sampling tables are built.
  int                   fNumPrimaryEnergies;
  // Minimum primary energy from which sampling tables are built
  double                fMinPrimaryEnergy;
  // The logarithm of the minimum primary energy.
  double                fLogMinPrimaryEnergy;
  // Inverse delta log kinetic energy (i.e. bin size in log-scale).
  double                fInvLogDeltaPrimaryEnergy;

  // The collection of sampling tables over the discrete primary energy grid.
  std::vector<SimLinAliasData*>  fTheTables;

};

#endif // SimKNTables_HH
