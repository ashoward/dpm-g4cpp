#ifndef SimMaxScatStrength_HH
#define SimMaxScatStrength_HH

//
// M. Novak: 2021
//
// Utility class to load the maximum scattering strength (K_1(E)~S_max(E)/tr1-mfp(E)))
// that is allowed to consume in a single (MSC) step in the reference material
// and provide (Spline) interpolated values at run time for (MSC) elastic steps.
//
// K_1(E) is loaded ONLY FOR THE REFERENCE MATERIAL (i.e. material with
// index 0). Note, that K_1(E) is not scalled (not divided) by the material
// density and it has no units [].
//
// It is assumed that the data has already been generated (by using the
// `InitElectronData::InitElasticData` method) and written into the
// `el_scatStrength.dat` file.

#include <vector>
#include <string>

#include "SimDataSpline.hh"

class SimMaxScatStrength {

public:

  SimMaxScatStrength() {}
 ~SimMaxScatStrength() {}

  void  LoadData(const std::string& dataDir, int verbose);

  // the maximum allowed scattering strength (K_1(E)) in a single MSC step in the
  // reference material approximated as `S_max(E)/tr1-mfp(E)` (no units) where
  // S_max(E) is a sigmoid-like function of the maximum allowed MSC step length
  // determined by the `slow` and `shigh` parameters with the smooth transition
  // around `ecross`.
  //
  double GetMaxScatStrength(double ekin) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fData.GetValue(e));
  }
  double GetMaxScatStrength(double ekin, double logekin) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fData.GetValue(e, logekin));
  }
  double GetMaxScatStrength(double ekin, int ilow) {
    return std::max(1.0E-20, fData.GetValueAt(ekin, ilow));
  }


private:

  // store the min/max kinetic enrgy values and the correspnding IMFP values
  double          fEmin;
  double          fEmax;

  // the IMFP data, divided by the ref. material density in [g/cm3], ready for
  // the run-time spline interpolation
  SimDataSpline   fData;
};

#endif // SimMaxScatStrength_HH
