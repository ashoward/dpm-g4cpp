#ifndef Configuration_HH
#define Configuration_HH

//
// M. Novak: 2021
//
// A configuration class that encapsulates the configuration parameters of the
// geometry and simulation. While the data generation part heavily relies on
// this infomation only the `fGeomIndex` pre-defined geometry index is used
// at run time since all other infomation is read form the generated data files.
//

#include <iostream>
#include <vector>
#include <string>

class Configuration {

public:

  // default number kinetic energy grid points for e-, photon data and the
  // min/max kinetic energy of the simualtion in [MeV]
  int    fNumEkinForElectron;
  int    fNumEkinForPhoton;
  double fMinEkin;
  double fMaxEkin;
  // secondary e-, gamma prodcution threshold energies [MeV] that also serve as
  // absorption energies or tracking cuts
  double fElectronCut;
  double fGammaCut;
  // the MSC maximum step length parameters [mm] [mm] and [MeV] units
  double fMscStepSlow;
  double fMscStepShigh;
  double fMscEcross;
  //
  // index of the pre-defined geometry: see more at `Geom.hh`
  // NOTE: this infomation is used only during the simulation phase but the
  //       the following list of materials is strongly connected to the pre-
  //       defined geometry since trhough the material indices. So we store this
  //       here where the index of the materials will be determined (by their
  //       order)
  int    fGeomIndex;
  //
  // list of G4-NIST material names that are in the geometry.
  // NOTE: their order will determine their index and the one with index 0 will
  //       be the reference material for the elastic interaction
  std::vector<std::string> fTheG4MaterialNames;
  //
  // setting the defaults:
  Configuration() {
    fNumEkinForElectron   =   128;
    fNumEkinForPhoton     =  1024;
    fMinEkin              = 0.001; //  1 [keV]
    fMaxEkin              =  21.0;  // 21 [MeV]
    //
    fElectronCut          =  0.20;  // 200 [keV]
    fGammaCut             =  0.05;  //  50 [keV]
    //
    fMscStepSlow          =   5.0;  //   5 [mm]
    fMscStepShigh         =  10.0;  //   1 [cm]
    fMscEcross            =  12.0;  //  12 [MeV]
    //
    fGeomIndex            =     0;  // homogeneous geometry (see `Geom.hh`)
    fTheG4MaterialNames.push_back("G4_WATER");
  }


  // a set of pre-defined configurations are available for specific simulations
  void PreDefinedConfigurations(int indx) {
    switch (indx) {
      // indx = 0 : homogeneous Water target all the default configuration
      case 0: break;
      // indx = 1 : 0-1 cm Water, 1-3 cm Air, 3 - ... cm Water with all default
      //            but the material list and the geometry configuration
      case 1: fGeomIndex = 1;
              fTheG4MaterialNames.resize(2);
              fTheG4MaterialNames[0] = "G4_WATER";
              fTheG4MaterialNames[1] = "G4_AIR";
              break;
      // indx = 2 : 0-2 cm Water, 2-4 cm Bone, 4 - ... cm Water with all default
      //            but the material list and the geometry configuration
      case 2: fGeomIndex = 2;
              fTheG4MaterialNames.resize(2);
              fTheG4MaterialNames[0] = "G4_WATER";
              fTheG4MaterialNames[1] = "G4_BONE_CORTICAL_ICRP";
              break;
      // indx = 3 : 0-1 cm Water, 1-1.5 cm Titatnium, 1.5 - 2.5 cm Bone and
      //            2.5 - ... cm Water with all default but the material list and
      //            the geometry configuration
      case 3: fGeomIndex = 3;
              fTheG4MaterialNames.resize(3);
              fTheG4MaterialNames[0] = "G4_WATER";
              fTheG4MaterialNames[1] = "G4_Ti";
              fTheG4MaterialNames[2] = "G4_BONE_COMPACT_ICRU";
              break;
      // indx = 4 : homogeneous Titatnium with all default but the material list
      case 4: fTheG4MaterialNames.resize(1);
              fTheG4MaterialNames[0] = "G4_Ti";
              break;
      // indx = 5 : homogeneous Tungsten with all default but the material list
      //            and the MSC maximum step lenght parameters
      case 5: fMscStepSlow  = 1.0; // [mm]
              fMscStepShigh = 1.0; // [mm]
              fTheG4MaterialNames.resize(1);
              fTheG4MaterialNames[0] = "G4_W";
              break;
      // default: homogeneous Water target all the default configuration
    }
  }


  //
  // Write the configuration:
  void Write() {
    std::cout << "\n === Pre-init/data generation configuration: "
              << "\n      - Min-Ekin          : " << fMinEkin << "  [MeV]"
              << "\n      - Max-Ekin          : " << fMaxEkin << "  [MeV]"
              << "\n      - #Ekins for e-     : " << fNumEkinForElectron
              << "\n      - #Ekins for gamma  : " << fNumEkinForPhoton
              << "\n      - cut for e-        : " << fElectronCut << "  [MeV]"
              << "\n      - cut for gamma     : " << fGammaCut    << "  [MeV]"
              << "\n      - MSC step param    : "
              << "\n             s-low    = " << fMscStepSlow  << "  [mm]"
              << "\n             s-high   = " << fMscStepShigh << "  [mm]"
              << "\n             e-cross  = " << fMscEcross    << " [MeV]"
              << "\n      - pre-defined geom. : " << fGeomIndex
              << "\n      - materials         : ";
    for (std::size_t i=0; i<fTheG4MaterialNames.size(); ++i) {
      std::cout<<"\n             [" << i << "] = " << fTheG4MaterialNames[i];
    }
    std::cout << "\n ............................................... \n"
              << std::endl;
  }
};

#endif
