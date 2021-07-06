#ifndef SimMaterialData_HH
#define SimMaterialData_HH

//
// M. Novak: 2021
//
// A simple object to store a very minimal set of data per material
// needed during the simulation

#include <string>
#include <vector>

class SimMaterialData {

public:
  SimMaterialData() {}
 ~SimMaterialData() {}

  void Load(const std::string& dataDir, int verbose=0);

public:

  //
  // additional infomation stored in the material file:
  //
  // secondary e- and gamma production threshold in [MeV]
  double                     fElectronCut;
  double                     fGammaCut;

  // MSC max-step lenght parameters
  double                     fMSCStepParSLow;
  double                     fMSCStepParSHigh;
  double                     fMSCStepParEcross;

  //
  // data per material:
  //
  int                        fNumMaterial;
  // name and density in [g/cm3]
  std::vector<std::string>   fMaterialName;
  std::vector<double>        fMaterialDensity; // in [g/cm3]

  // the material scaling factor for the Moller inverse-mf: [A/(Z\rho/)]_ref [(Z\rho)/A]_actual
  // or more exactly its [A/Z)]_ref [(Z)/A]_actual part
  std::vector<double>        fMollerIMFPScaling;

};

#endif // SimMaterialData_HH
