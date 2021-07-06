#ifndef SimGSTables_HH
#define SimGSTables_HH

//
// M. Novak: 2021
//
// Sampling tables for generating (cosine of) angular deflection of electrons
// due to multiple Coulomb scattering along a given (fixed) path length in a
// given (reference) material.
//
// It is assumed, that the required data has been generated and written into
// the `el_GSDtrData.dat` file previously by using `InitElectronData::InitGSData`.
// This file contains GS sampling tables for all the materials that are in the
// target but the DATA ARE TAKEN ONLY FOR THE VERY FIRST MATERIAL, that is
// considered to be THE REFERENCE MATERIAL. The file contains infomation regarding
// the discrete electron energy grid used to build tables as well as the number
// of discrete cumulative values. The former determines the number of tables for
// a given material while the latter gives the size of the individual tables.
// For each materials, the data building up the tables are listed at each
// individual discrete electron energies as follows:
//  - table-index, electron-energy, transformation-paramater values
//  - 3 x #cumulatives values that compose the sampling tables
// Note, that these sampling tabels are different than those used for SB, KN or
// the Moller interaction since this sampling is based on the numerical inversion
// of the cumulative of a transformed (as close to uniform as possible) pdf.
//
// All the data can be loaded by using the provided `LoadData(dataDir)` method
// then `SampleAngularDeflection()` can be used to generate cosine of the angular
// deflection along the predefined, fixed path length in the given reference
// material without rejection. The path length is predefined by the `sLow`, `sHigh`
// and `eCross` parameters that must have been given during the data preparation.


#include <vector>
#include <string>


class SimGSTables {

public:

  // CTR and DTR
  SimGSTables();
 ~SimGSTables() { CleanTables(); }

  // loads the previously prepared data from the appropriate file locating in
  // the directory specified by the `dataDir` input argument
  void   LoadData(const std::string& dataDir, int verbose);

  // samples cosine of the angular deflection at the given `ekin` electron energy
  // (travelling the predefined path length in the predefined referecne material).
  // The two additional input arguments are uniformly random values on [0,1].
  //
  // NOTE: it is assumed that the `eprim` electron energy: electron-cut < eprim <E_max
  double SampleAngularDeflection(double ekin, double rndm1, double rndm2);


private:

  void   CleanTables();


private:

  // The size of the individual sampling tables
  int                   fSamplingTableSize;
  // The number of discrete primary energies at which sampling tables are built.
  int                   fNumPrimaryEnergies;
  // Minimum primary energy from which sampling tables are built
  double                fMinPrimaryEnergy;
  // The electron kinetic energy grid
  std::vector<double>   fPrimaryEnergyGrid;  // [#fNumPrimaryEnergies]
  // The logarithm of the minimum primary energy.
  double                fLogMinPrimaryEnergy;
  // Inverse delta log kinetic energy (i.e. bin size in log-scale).
  double                fInvLogDeltaPrimaryEnergy;

  // Delta cumulative value: 1/(fSamplingTableSize-1)
  double                fDeltaCum;

  struct OnePoint {
    double fVarU;     // the transfomred variable `u` that corresponds this cumulative value
    double fParmA;    // interpolation parameters for the approximate inversion of the cumulative
    double fParmB;
  };
  struct OneGSTable {
    double                 fTransformParam; // transformation paraneter `a` used to generate the corresponding smooth pdf
    std::vector<OnePoint>  fGSTable;        // the representation of the inverse cumulative of this pdf
  };

  // The collection of sampling tables over the discrete primary energy grid for this material.
  std::vector<OneGSTable*> fTheTables;

};

#endif // SimGSTables_HH
