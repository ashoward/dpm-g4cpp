#ifndef SimMollerTables_HH
#define SimMollerTables_HH

//
// M. Novak: 2021
//
// Sampling tables for generating energy transfer values in Moller interaction
// without rejection.
//
// It is assumed, that the required data has been generated and written into
// the `ioni_MollerDtrData.dat` file previously by using the `MollerTableBuilder`.
// This file contains a discrete primary energy grid and a sampling table
// (`SimLinAliasData`) at each of the discrete primary energies. All the data
// can be loaded by using the provided `LoadData(dataDir)` method then the
// `SampleEnergyTransfer()` can be used to generate values of energy transferred
// to the secondary electron in Moller interaction without rejection.


#include <vector>
#include <string>

class SimLinAliasData;


class SimMollerTables {

public:

  // CTR and DTR
  SimMollerTables();
 ~SimMollerTables() { CleanTables(); }

  // loads the previously prepared data from the appropriate file locating in
  // the directory specified by the `dataDir` input argument
  void   LoadData(const std::string& dataDir, int verbose);

  // samples value of energy transferred to the secondary electron in Moller
  // interaction at the given `eprim` primary electron energy. The three
  // additional input arguments are uniformly random values on [0,1].
  //
  // NOTE: it is assumed that: 2 x electron-cut < eprim < E_max (also for e+)
  double SampleEnergyTransfer(double eprim, double rndm1, double rndm2, double rndm3);



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

#endif // SimMollerTables_HH
