
#include "G4Setup.hh"
#include "ElectronData.hh"
#include "InitElectronData.hh"
#include "PhotonData.hh"
#include "InitPhotonData.hh"

#include "Configuration.hh"

#include <iostream>
#include <iomanip>
#include <string>


#include <getopt.h>

//
// M. Novak: 2021
//
// The pre-init/data generationpart of the `dpm++` prototype. This needs to be
// exeuted in order to generate all the data files that are required by the
// simulation phase. This executable depends on `Geant4` while the second, that
// performs the simulation based on the data generated here, is independent of
// anything.
//

//
// Input argument options for the pre-init/data generation phase with their
// default values.
//
// location where the generated data should be written (the dir. must exist)
static std::string gOutputDataDir("./data");
// one of the pre-defined configuration settings index (see more in `Configuration.hh`)
static int         gConfigIndex = 0;
//
static struct option options[] = {
  {"input-data-dir      (where the pre-generated data will be written)  - default: ./data" , required_argument, 0, 'd'},
  {"configuration-index (one of the pre-defined configuration settings) - default: 0"      , required_argument, 0, 'c'},
  {"help"                                                                                  , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
// auxiliary functions fot obtaining input arguments
void Help();
void GetOpt(int argc, char *argv[]);


//
//
// The main: obtaines the input arguments, creates the selected configuration
//           settings, generates all the data needed for the simulation and
//           writes them out into files under the specified location.
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
  // Start the data generation:
  //
  // Builds the materials, initialises a fake G4 geometry with having each of
  // the materials in a different detector regions. This later is required since
  // G4 requires the secondary production threshold in length and not in energy.
  // Therefore, a fixed energy threshold means different range and absorption
  // lengths in the different materials: these lengths are also computed and
  // the appropriate lengths, leading the given production threshold in energy,
  // is set for e- (range) and gamma (absorption length) in the given detector
  // region.
  // Note: this step is required only for being able to initialise the G4 models
    FakeG4Setup ( theConfig.fTheG4MaterialNames,
                  theConfig.fElectronCut,
                  theConfig.fGammaCut
                );
  //
  // message
  std::cout << "\n === Data are generated now (don't care about the "
            << "\n     Geant4 warning: 'WARNING from G4PenelopeIonisationModel...')\n"
            << std::endl;
  //
  // Create an electron data object to store the related elastic and energy loss
  // related data per material:
  // - elastic data:     elastic and first transport MFP, screening parameter
  // - energy loss data: (restricted) stopping power (radiative plus collisonal) and range
  // - sampling tables
  // - etc.
  // NOTE: the common kinetic energy grid is from the secondary e- production threshold
  //       till the maximum kinetic energy specified above
  ElectronData theElectronData ( "../SBTables/data",
                                 theConfig.fTheG4MaterialNames.size(),
                                 theConfig.fElectronCut*0.99,
                                 theConfig.fMaxEkin,
                                 theConfig.fNumEkinForElectron
                               );
  // Fill  this above electron data object
  InitElectronData ( theElectronData,
                     theConfig.fElectronCut,
                     theConfig.fGammaCut,
                     theConfig.fMscStepShigh,
                     theConfig.fMscStepSlow,
                     theConfig.fMscEcross
                   );

  // Create the photon data object and fill it
  PhotonData thePhotonData( theConfig.fTheG4MaterialNames.size(),
                            theConfig.fGammaCut*0.99,
                            theConfig.fMaxEkin,
                            theConfig.fNumEkinForPhoton
                          );
  InitPhotonData ( thePhotonData,
                   theConfig.fElectronCut,
                   theConfig.fGammaCut
                 );

  //
  // message
  std::cout << "\n === All data have been generated and will be written now "
            << "\n     into files under the directory: " << gOutputDataDir << "\n"
           << std::endl;
  //
  // Write out all the data for simulation
  theElectronData.WriteData(gOutputDataDir);
  thePhotonData.WriteData(gOutputDataDir);

return 0;
}



//
// Inplementation of the auxiliary functions for obtaining input ragumets
//
void Help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  The dpm++ pre-init/data generation phase."
           << std::endl;
  std::cout<<"\n  Usage: dpm_GenerateData [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}


void GetOpt(int argc, char *argv[]) {
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "hd:c:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
       c = options[optidx].val;
       /* fall through */
    case 'd':
       gOutputDataDir = optarg;
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
  std::cout << "\n === The dpm++ data generation configuration: \n"
            << "\n     - output data directory : " << gOutputDataDir
            << "\n     - configuration index   : " << gConfigIndex
            << std::endl;
}
