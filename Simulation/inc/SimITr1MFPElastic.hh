#ifndef SimITr1MFPElastic_HH
#define SimITr1MFPElastic_HH

//
// M. Novak: 2021
//
// Utility class to load inverse first transport mean free path (ITr1MFP) and
// provide (Spline) interpolated values at run time for elatsic interactions.
//
// The ITr1MFP is loaded FOR ALL THE INDIVIDUAL MATERIALS and SCALLED (divided)
// by the corresponding material DENSITY IN [g/cm3].
// So eventually it will make sence to use the values only after scaling back
// (multiplying) by the actual material density in [g/cm3]. After this scaling,
// the ITr1MFP will have the correct [1/mm] units.
//
// It is assumed that the data has already been generated (by using the
// `InitElectronData::InitElasticData` method) and written into the `el_itr1mfp.dat`
// file already divided by the (reference) material denisty in [g/cm3].

#include <vector>
#include <string>

#include "SimDataSpline.hh"

class SimITr1MFPElastic {

public:

  SimITr1MFPElastic();
 ~SimITr1MFPElastic() {}

  void  LoadData(const std::string& dataDir, int verbose);

  // the inverse Tr1-MFP in [1/mm] [cm3/g] scalled units
  double GetITr1MFPPerDensity(double ekin, int imat) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValue(e));
  }
  double GetITr1MFPPerDensity(double ekin, double logekin, int imat) {
    // make sure that E_min <= ekin < E_max
    const double e = std::min(fEmax-1.0E-6, std::max(fEmin, ekin));
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValue(e, logekin));
  }
  double GetITr1MFPPerDensity(double ekin, int ilow, int imat) {
    return std::max(1.0E-20, fDataPerMaterial[imat].GetValueAt(ekin, ilow));
  }


private:

  int             fNumMaterial;

  // store the min/max kinetic enrgy values and the correspnding IMFP values
  double          fEmin;
  double          fEmax;

  // the inverse Tr1-MFP data per materail, divided by the  material density in
  // [g/cm3], ready for the run-time spline interpolation
  std::vector<SimDataSpline>   fDataPerMaterial;
};

#endif // SimITr1MFPElastic_HH
