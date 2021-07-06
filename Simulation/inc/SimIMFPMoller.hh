#ifndef SimIMFPMoller_HH
#define SimIMFPMoller_HH

//
// M. Novak: 2021
//
// Utility class to load (restricted) IMFP (macroscopic cross section) and provide
// (Spline) interpolated values at run time for (above threshold) Moller interactions.
//
// Moller IMFP is loaded ONLY FOR THE REFERENCE MATERIAL (i.e. material with
// index 0) and SCALLED (divided) by the reference material DENSITY IN [g/cm3].
// So eventually it will make sence to use the values only after scaling back
// (multiplying) by the actual material density in [g/cm3]. After this scaling,
// the IMFP will have the correct [1/mm] units.
//
// It is assumed that the data has already been generated (by using the
// `InitElectronData::InitElossData` method) and written into the `imfp_moller.dat`
// file already divided by the (reference) material denisty in [g/cm3].

#include <vector>
#include <string>

#include "SimDataSpline.hh"

class SimIMFPMoller {

public:

  SimIMFPMoller() {}
 ~SimIMFPMoller() {}

  void  LoadData(const std::string& dataDir, int verbose);

  // the inverse MFP in [1/mm] [cm3/g] scalled units
  double GetIMFPPerDensity(double ekin) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fData.GetValue(e));
  }
  double GetIMFPPerDensity(double ekin, double logekin) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fData.GetValue(e, logekin));
  }
  double GetIMFPPerDensity(double ekin, int ilow) {
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

#endif // SimIMFPMoller_HH
