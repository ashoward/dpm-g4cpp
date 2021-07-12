
#include "SimMaterialData.hh"
#include "SimElectronData.hh"
#include "SimPhotonData.hh"

#include "SimDPMLike.hh"

#include "Configuration.hh"

#include <iostream>
#include <iomanip>
#include <string>

#include <getopt.h>

//
// M. Novak: 2021
//
// The simulation phase of the `dpm++` prototype. This can be executed only after
// the pre-init/data generation part (i.e. `dpm_GenerateData`) since this relies
// on all the data files generated during the pre-init phase while performs the
// simulation of with the specified configuration using the given input arguments.
//

//
// Input arguments for the simulation phase with their default values.
static std::string   gPrimaryParticle("e-");        // primary particle is electron
static std::string   gInputDataDir("./data");       // location of the pre-generated data
static double        gPrimaryEnergy     = 15.0;     // primary particle energy in [MeV]
static double        gVoxelSize         =  1.0;     // geometry voxel/box size in [mm]
static int           gNumPrimaries      =  1.0E+5;  // simulate 100 000 primary events
static int           gConfigIndex       =  0;       // 0 that corresponds to a homogeneous Water geometry
//
static struct option options[] = {
  {"primary-particle    (possible particle names: e-, gamma)            - default: e-"     , required_argument, 0, 'p'},
  {"primary-energy      (in [MeV] units)                                - default: 15"     , required_argument, 0, 'e'},
  {"number-of-histories (number of primary events to simulate)          - default: 1.0E+5" , required_argument, 0, 'n'},
  {"input-data-dir      (where the pre-generated data are located)      - default: ./data" , required_argument, 0, 'd'},
  {"voxel-size          (size of the voxel/box in [mm])                 - default: 1.0"    , required_argument, 0, 'b'},
  {"configuration-index (one of the pre-defined configuration index)    - default: 0"      , required_argument, 0, 'c'},
  {"help"                                                                                  , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
// auxiliary functions for obtaining input arguments
void Help();
void GetOpt(int argc, char *argv[]);


//
//
// The main: obtaines the input arguments, creates the selected configuration
//           settings (mainly for infomation only since only its `fGeomIndex`
//           member is need now), loads all the data needed for the simulation
//           and executes the simulation with the given configuration and input
//           arguments. At termination, the `hist.sim` file contains the simulated
//           depth dose distribution.
int main(int argc, char *argv[]) {
  //
  // get the input arguments
  GetOpt(argc, argv);
  //
  // get the predefined configuration
  Configuration  theConfig;
  theConfig.PreDefinedConfigurations(gConfigIndex);
  theConfig.Write();
  //
  // Load data for simulation:
  // - electron related data
  SimElectronData theSimElectronData;
  theSimElectronData.Load(gInputDataDir);
  // - photon related data
  SimPhotonData theSimPhotonData;
  theSimPhotonData.Load(gInputDataDir);
  // - configuration and material related data
  SimMaterialData theSimMaterialData;
  theSimMaterialData.Load(gInputDataDir);
  //
  // Execute the simulation according to the iput arguments
  bool isElectron = (gPrimaryParticle=="e-");
  Simulate(gNumPrimaries, gPrimaryEnergy, isElectron, gVoxelSize, theSimMaterialData, theSimElectronData, theSimPhotonData, theConfig.fGeomIndex);
  //
  return 0;
}


//
// Inplementation of the auxiliary functions for obtaining input ragumets
//
void Help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  The dpm++ simulation phase."
           << std::endl;
  std::cout<<"\n  Usage: dpm_Simulate [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}


void GetOpt(int argc, char *argv[]) {
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "hp:e:n:d:b:c:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
       c = options[optidx].val;
       /* fall through */
    case 'p':
       gPrimaryParticle = optarg;
       if ( !(gPrimaryParticle=="e-" || gPrimaryParticle=="gamma") ) {
         std::cout << " *** Unknown primary particle name -p: " << optarg << std::endl;
         Help();
         exit(-1);
       }
       break;
    case 'e':
       gPrimaryEnergy = std::stod(optarg);
       break;
    case 'n':
       gNumPrimaries  = std::stoi(optarg);
       break;
    case 'd':
       gInputDataDir  = optarg;
       break;
    case 'b':
       gVoxelSize     = std::stod(optarg);
       break;
    case 'c':
       gConfigIndex   = std::stoi(optarg);
       break;
    case 'h':
       Help();
       exit(-1);
       break;
    default:
      printf(" *** Unknown input argument: %c\n",c);
      Help();
      exit(-1);
    }
   }
   std::cout << "\n === The dpm++ simulation confguration: \n"
            << "\n     - input data directory  : " << gInputDataDir
            << "\n     - primary particle      : " << gPrimaryParticle
            << "\n     - primary energy        : " << gPrimaryEnergy << " [MeV]"
            << "\n     - number of histories   : " << gNumPrimaries
            << "\n     - geometry voxel size   : " << gVoxelSize  << " [mm]"
            << "\n     - confoguration index   : " << gConfigIndex
            << std::endl;
}
