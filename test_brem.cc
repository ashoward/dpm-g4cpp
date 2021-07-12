//
// == M. Novak: 2021
//
// This model-level test can be used to test the rejection-free sampling of the
// energy distribution, transferred to the emitted photon in case of e-
// bremsstrahlung, according to the Seltzer-Berger DCS.
//
// NOTE: `dpm_GenerateData` needs to be executed before since this test rely on
//       the generated data. Furthermore, exatly the same function ( the
//       `SimDPMLike::PerformBrem()`) that is used during teh simualtion, is
//       used by the test to generate samples according to the DCS determined
//       by the primary e- primary energy and the material. The materials, as
//       well as the secondary gamma productio threshold, are those used during
//       the data generation and can be selected by their index(that corresponds
//       to their order of definition).
//
// NOTE: if the matrial index is valid or the primary energy is indeed higher
//       than the gamma cut, ect. are not tested (since developers should be
//       well aware of these constraints of a restricted brem. description).


#include "SimDPMLike.hh"
#include "SimSBTables.hh"
#include "SimMaterialData.hh"
#include "Random.hh"
#include "Track.hh"
#include "TrackStack.hh"


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdio>

#include <getopt.h>

static std::string   gInputDataDir("./data");       // location of the pre-generated data
static double        gPrimaryEnergy     =  12.345;  // primary particle energy in [MeV] (> gamma-cut)
static int           gNumPrimaries      =  1.0E+5;  // generate 100 000 samples from the brem. interaction
static int           gMaterialIndex     =       0;  // material index (one of those data generated for)
//
static struct option options[] = {
  {"primary-energy      (in [MeV] units)                                - default: 12.345" , required_argument, 0, 'e'},
  {"number-of-primaries (number of interactions to sample)              - default: 1.0E+5" , required_argument, 0, 'n'},
  {"input-data-dir      (where the pre-generated data are located)      - default: ./data" , required_argument, 0, 'd'},
  {"material-index      (one of the material index data generated for)  - default: 0"      , required_argument, 0, 'm'},
  {"help"                                                                                  , no_argument      , 0, 'h'},
  {0, 0, 0, 0}
};
// auxiliary functions for obtaining input arguments
void Help();
void GetOpt(int argc, char *argv[]);


int main(int argc, char *argv[]) {
  //
  // get the input arguments
  GetOpt(argc, argv);
  //
  // Load Material data
  SimMaterialData theSimMaterialData;
  theSimMaterialData.Load(gInputDataDir);
  // Load SB sampling tables
  SimSBTables theSBTables;
  theSBTables.LoadData(gInputDataDir, 1);
  //
  const double gammaCut = theSimMaterialData.fGammaCut;
  //
  // create a histogram of the x =  log10(k/E_0)
  int    hbins = 101;
  double  xmin = std::log10(gammaCut/gPrimaryEnergy);
  double  xmax = 0.1;
  double  hdel = (xmax-xmin)/(hbins-1);
  double ihdel = 1./hdel;
  std::vector<double> theHist(hbins, 0.);
  // create a primary and secondary track
  Track  primaryTrack;
  Track  secondaryTrack;
  for (int is=0; is<gNumPrimaries; ++is) {
    // set initial properties of the primary electron (that are updated in the
    // interaction) according to the input arguments
    // NOTE: material indes do not changed but keep them together for clarity
    primaryTrack.fEkin         = gPrimaryEnergy;
    primaryTrack.fMatIndx      = gMaterialIndex;
    primaryTrack.fDirection[0] = 0.;
    primaryTrack.fDirection[1] = 0.;
    primaryTrack.fDirection[2] = 1.;
    //
    // use the fuction from SimDPMLike to perform the brem intercation
    PerformBrem(primaryTrack, &theSBTables);
    // NOTE: in each interactions, secondary tracks are pushed to the global
    //       TrackStack so get the secondary track
    TrackStack::Instance().PopIntoThisTrack(secondaryTrack);
    // get the secondary (emitted photon) energy
    const double theK  =  secondaryTrack.fEkin;
    // compute the reduced photon energy k/E_0
    double redPhEnergy = theK/gPrimaryEnergy;
    if (redPhEnergy>0.0) {
      redPhEnergy = std::log10(redPhEnergy);
      theHist[(int)((redPhEnergy-xmin)*ihdel)] += 1.0;
    }
  }
  // normalisation
  double norm = 0.0;
  for (int ih=0; ih<hbins-1; ++ih) {
    norm += 0.5*hdel*(theHist[ih]+theHist[ih+1]);
  }
  norm = 1./norm;
  FILE* f = fopen("res_brem_test_ph_energy.dat","w");
  for (int ih=0; ih<hbins-1; ++ih) {
    fprintf(f, "%d %lg %lg\n", ih, xmin+(ih+0.5)*hdel, theHist[ih]*norm);
  }
  fclose(f);

 return 0;
}


//
// Inplementation of the auxiliary functions for obtaining input ragumets
//
void Help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  A bremsstrahlung interaction test."
           << std::endl;
  std::cout<<"\n  Usage: test_brem [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}


void GetOpt(int argc, char *argv[]) {
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "h:e:n:d:m:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
       c = options[optidx].val;
       /* fall through */
    case 'e':
       gPrimaryEnergy = std::stod(optarg);
       break;
    case 'n':
       gNumPrimaries  = std::stoi(optarg);
       break;
    case 'd':
       gInputDataDir  = optarg;
       break;
    case 'm':
       gMaterialIndex = std::stoi(optarg);
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
   std::cout << "\n === The test_brem simulation confguration: \n"
            << "\n     - input data directory  : " << gInputDataDir
            << "\n     - primary energy        : " << gPrimaryEnergy << " [MeV]"
            << "\n     - number of interaction : " << gNumPrimaries
            << "\n     - material index        : " << gMaterialIndex
            << std::endl;
}
