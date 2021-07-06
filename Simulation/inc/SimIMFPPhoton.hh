#ifndef SimIMFPPhoton_HH
#define SimIMFPPhoton_HH

//
// M. Novak: 2021
//
// Utility class to load IMFP for photon interactions including Compton, Pair -
// production and total (sum of Compton, Pair-prodcution, Photoelectric) and
// provide (linearly) interpolated values at run time. The type of the IMFP is
// eventually determined by the `type` integer received at construction:
//  - `type = 0` - total,
//  - `type = 1` - Compton,
//  - `type = 2` - Pair-production,
//
// All IMFP are loaded FOR ALL THE INDIVIDUAL MATERIALS and SCALLED (divided)
// by the corresponding material DENSITY IN [g/cm3].
// So eventually it will make sence to use the values only after scaling back
// (multiplying) by the actual (voxel) material density in [g/cm3]. After this
// scaling, the IMFP will have the correct [1/mm] units.
//
// It is assumed that the data has already been generated (by using the
// `InitPhotonData::InitMXsecData` method) and written into the `imfp_total.dat`
// file already divided by the material denisty in [g/cm3].

#include <vector>
#include <string>

#include "SimDataLinear.hh"

class SimIMFPPhoton {

public:

  SimIMFPPhoton(int type);
 ~SimIMFPPhoton();

  void  LoadData(const std::string& dataDir, int verbose);

  // the inverse IMFP in [1/mm] [cm3/g] scalled units
  double GetIMFPPerDensity(double ekin, int imat) {
    // check vacuum case i.e. imat = -1
    if (imat<0) return 1.0E-20;
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValue(e));
  }
  double GetIMFPPerDensity(double ekin, int ilow, int imat) {
    // check vacuum case i.e. imat = -1
    if (imat<0) return 1.0E-20;
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValueAt(ekin, ilow));
  }


private:
  // number of materials (photon data are taken to the actual material during the sim.)
  int     fNumMaterial;

  // type of the IMFP (total, Compton or Pair-production)
  int     fType;

  // store the min/max kinetic enrgy values and the correspnding IMFP values
  double  fEmin;
  double  fEmax;

  // the IMFP data per material, divided by the material density in [g/cm3],
  // ready for the run-time linear interpolation
  std::vector<SimDataLinear>   fDataPerMaterial;
};

#endif // SimIMFPPhoton_HH
