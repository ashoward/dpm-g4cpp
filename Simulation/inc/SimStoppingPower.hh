#ifndef SimStoppingPower_HH
#define SimStoppingPower_HH

//
// M. Novak: 2021
//
// Utility class to load (restricted) total (radiative and collisonal) Stopping
// Power and provide (Spline) interpolated values at run time for average energy
// loss due to sub-threshold ionisation and bremsstrahlung interactions.
//
// The dE/dx is loaded FOR ALL THE INDIVIDUAL MATERIALS and SCALLED (divided) by
// the corresponding material DENSITY IN [g/cm3].
// So eventually it will make sence to use the values only after scaling back
// (multiplying) by the actual (voxel) material density in [g/cm3]. After this
// scaling, the dE/dx will have the correct [MeV/mm] units.
//
// It is assumed that the data has already been generated (by using the
// `InitElectronData::InitElossData` method) and written into the `eloss_rdedx.dat`
// file already divided by the material denisty in [g/cm3].

#include <vector>
#include <string>

#include "SimDataSpline.hh"

class SimStoppingPower {

public:

  SimStoppingPower();
 ~SimStoppingPower();

  void  LoadData(const std::string& dataDir, int verbose);

  // the stopping power in [MeV/mm] [cm3/g] scalled units
  double GetDEDXPerDensity(double ekin, int imat) {
    // check vacuum case i.e. imat = -1
    if (imat<0) return 1.0E-20;
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValue(e));
  }
  double GetDEDXPerDensity(double ekin, double logekin, int imat) {
    // check vacuum case i.e. imat = -1
    if (imat<0) return 1.0E-20;
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValue(e, logekin));
  }
  double GetDEDXPerDensity(double ekin, int ilow, int imat) {
    // check vacuum case i.e. imat = -1
    if (imat<0) return 1.0E-20;
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValueAt(ekin, ilow));
  }


private:

  // number of materials (dE/dx data are used for the given material)
  int     fNumMaterial;

  // store the min/max kinetic enrgy values and the correspnding dE/dx values
  double  fEmin;
  double  fEmax;

  // the dE/dx data per material, divided by the material density in [g/cm3],
  // ready for the run-time spline interpolation
  std::vector<SimDataSpline>   fDataPerMaterial;
};

#endif // SimStoppingPower_HH
