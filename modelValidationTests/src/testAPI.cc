//
// Created by mbarbone on 7/22/21.
//

#include "testAPI.h"
#include "SimSBTables.hh"

#include <cmath>
#include <numeric>
#include <Track.hh>
#include <SimDPMLike.hh>
#include <TrackStack.hh>
#include <utility>
#include <SimMollerTables.hh>
#include <SimGSTables.hh>

AbstractTest::AbstractTest(std::string g_input_data_dir,
                           double gPrimaryEnergy,
                           int gNumPrimaries,
                           int gMaterialIndex)
  : gInputDataDir(std::move(g_input_data_dir)),
    gPrimaryEnergy(gPrimaryEnergy),
    gNumPrimaries(gNumPrimaries),
    gMaterialIndex(gMaterialIndex) {
  theSimMaterialData.Load(gInputDataDir);
  // create a histogram of the x =  log10(k/E_0)
  const double gammaCut = theSimMaterialData.fGammaCut;
  hbins = 101;
  xmin = std::log10(gammaCut / gPrimaryEnergy);
  xmax = 0.1;
  hdel = (xmax - xmin) / (hbins - 1);
  ihdel = 1. / hdel;
  theHist = std::vector<double>(hbins, 0.);
  theCosineHist = std::vector<double>(hbins, 0.);
}

void AbstractTest::writeHists(const std::string &filenamePrefix) {
  double norm = 0.0;
  for (int ih = 0; ih < hbins - 1; ++ih) {
    norm += 0.5 * hdel * (theHist[ih] + theHist[ih + 1]);
  }
  norm = 1. / norm;

  FILE *f = fopen((filenamePrefix + "sec_energy.dat").c_str(), "w");
  for (int ih = 0; ih < hbins - 1; ++ih) {
    fprintf(f, "%d %lg %lg\n", ih, xmin + (ih + 0.5) * hdel, theHist[ih] * norm);
  }
  fclose(f);
  const double cos_norm =
    1. / std::accumulate(theCosineHist.begin(), theCosineHist.end(), 0.);

  FILE *f_cos = fopen((filenamePrefix + "_cost.dat").c_str(), "w");
  const auto step  = 2. / (hbins-1);
  auto position = -1.;
  for (int ih = 0; ih < hbins; ++ih) {
    fprintf(f_cos, "%d %lg %lg\n", ih, position, theCosineHist[ih] * cos_norm);
    position+=step;
  }
  fclose(f_cos);
}

__attribute__((unused)) std::vector<std::tuple<double, double>> AbstractTest::getEnergyHist() {
  std::vector<std::tuple<double, double>> normalized;
  double norm = 0.0;
  for (int ih = 0; ih < hbins - 1; ++ih) {
    norm += 0.5 * hdel * (theHist[ih] + theHist[ih + 1]);
  }
  norm = 1. / norm;
  normalized.reserve(hbins - 1);
  for (int ih = 0; ih < hbins - 1; ++ih) {
    normalized.emplace_back(xmin + (ih + 0.5) * hdel, theHist[ih] * norm);
  }
  return normalized;
}

__attribute__((unused)) std::vector<std::tuple<double, double>> AbstractTest::getCosHist() {
  std::vector<std::tuple<double, double>> normalized;
  const auto step  = 2. / (hbins-1);
  auto position = -1.;
  const double cos_norm =
    1. / std::accumulate(theCosineHist.begin(), theCosineHist.end(), 0.);
  normalized.reserve(hbins);
  for (int ih = 0; ih < hbins; ++ih) {
    normalized.emplace_back(position, theCosineHist[ih] * cos_norm);
    position+=step;
  }
  return normalized;
}

void BremTest::simulate() {
  SimSBTables theSBTables;
  theSBTables.LoadData(gInputDataDir, 1);
  Track primaryTrack;
  Track secondaryTrack;
  for (int is = 0; is < gNumPrimaries; ++is) {
    // set initial properties of the primary electron (that are updated in the
    // interaction) according to the input arguments
    // NOTE: material indes do not changed but keep them together for clarity
    primaryTrack.fEkin = gPrimaryEnergy;
    primaryTrack.fMatIndx = gMaterialIndex;
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
    const double theK = secondaryTrack.fEkin;
    const double theCost = secondaryTrack.fDirection[2];
    theCosineHist[(int) ((theCost * (hbins/2)) + (hbins/2))] += 1;

    // compute the reduced photon energy k/E_0
    double redPhEnergy = theK / gPrimaryEnergy;
    if (redPhEnergy > 0.0) {
      redPhEnergy = std::log10(redPhEnergy);
      theHist[(int) ((redPhEnergy - xmin) * ihdel)] += 1.0;
    }
  }
}

void BremTest::writeHists() {
  AbstractTest::writeHists("../output/res_brem_test_");
}

void MollerTest::simulate() {
  SimMollerTables theMollerTable;
  theMollerTable.LoadData(gInputDataDir, 1);
  Track primaryTrack;
  Track secondaryTrack;
  for (int is = 0; is < gNumPrimaries; ++is) {
    // set initial properties of the primary electron (that are updated in the
    // interaction) according to the input arguments
    // NOTE: material indes do not changed but keep them together for clarity
    primaryTrack.fEkin = gPrimaryEnergy;
    primaryTrack.fMatIndx = gMaterialIndex;
    primaryTrack.fDirection[0] = 0.;
    primaryTrack.fDirection[1] = 0.;
    primaryTrack.fDirection[2] = 1.;
    //
    // use the fuction from SimDPMLike to perform the brem intercation
    PerformMoller(primaryTrack, &theMollerTable);
    // NOTE: in each interactions, secondary tracks are pushed to the global
    //       TrackStack so get the secondary track
    TrackStack::Instance().PopIntoThisTrack(secondaryTrack);
    // get the secondary (emitted photon) energy
    const double theK = secondaryTrack.fEkin;
    const double theCost = secondaryTrack.fDirection[2];
    theCosineHist[(int) ((theCost * (hbins/2)) + (hbins/2))] += 1;

    // compute the reduced photon energy k/E_0
    double redPhEnergy = theK / gPrimaryEnergy;
    if (redPhEnergy > 0.0) {
      redPhEnergy = std::log10(redPhEnergy);
      theHist[(int) ((redPhEnergy - xmin) * ihdel)] += 1.0;
    }
  }
}
void MollerTest::writeHists() {
  AbstractTest::writeHists("../output/res_moller_test_");
}
void MSCAngularDeflectionTest::simulate() {
  SimGSTables theGSTables;
  theGSTables.LoadData(gInputDataDir, 1);
  Track primaryTrack;
  for (int is = 0; is < gNumPrimaries; ++is) {
    // set initial properties of the primary electron (that are updated in the
    // interaction) according to the input arguments
    // NOTE: material indes do not changed but keep them together for clarity
    primaryTrack.fEkin = gPrimaryEnergy;
    primaryTrack.fMatIndx = gMaterialIndex;
    primaryTrack.fDirection[0] = 0.;
    primaryTrack.fDirection[1] = 0.;
    primaryTrack.fDirection[2] = 1.;
    //
    // use the fuction from SimDPMLike to perform the brem intercation
    PerformMSCAngularDeflection(primaryTrack, gPrimaryEnergy, &theGSTables);
    // NOTE: in each interactions, secondary tracks are pushed to the global
    //       TrackStack so get the secondary track
    // get the secondary (emitted photon) energy
    const double theK = primaryTrack.fEkin;
    const double theCost = primaryTrack.fDirection[2];
    theCosineHist[(int) ((theCost * (hbins/2)) + (hbins/2))] += 1;

    // compute the reduced photon energy k/E_0
    double redPhEnergy = theK / gPrimaryEnergy;
    if (redPhEnergy > 0.0) {
      redPhEnergy = std::log10(redPhEnergy);
      theHist[(int) ((redPhEnergy - xmin) * ihdel)] += 1.0;
    }
  }
}
void MSCAngularDeflectionTest::writeHists() {
  AbstractTest::writeHists("../output/res_mscAngularDeflection_test_");
}
