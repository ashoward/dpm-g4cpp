//
// Created by mbarbone on 7/22/21.
//

#ifndef DPM_MODELVALIDATIONAPI_TESTAPI_H_
#define DPM_MODELVALIDATIONAPI_TESTAPI_H_

#include <string>
#include <tuple>
#include <utility>
#include <vector>


class AbstractTest {
 protected:
  const std::string gInputDataDir;
  const double gPrimaryEnergy;
  const int gNumPrimaries;
  const int gMaterialIndex;
  // histogram parameter
  int hbins;
  double xmin;
  double xmax;
  double hdel;
  double ihdel;
  std::vector<double> theHist;
  std::vector<double> theCosineHist;
 public:
  explicit AbstractTest(std::string gInputDataDir = "../data",
                        double gPrimaryEnergy = 12.345,
                        int gNumPrimaries = 1.0E+5,
                        int gMaterialIndex = 0);
  void writeHists(const std::string &filenamePrefix);
  __attribute__((unused)) std::vector<std::tuple<double, double>> getEnergyHist();
  __attribute__((unused)) std::vector<std::tuple<double, double>> getCosHist();
  virtual void simulate() = 0;
};

class BremTest : public AbstractTest {
 public:
  explicit BremTest(std::string
                    gInputDataDir = "../data",
                    double gPrimaryEnergy = 12.345,
                    int gNumPrimaries = 1.0E+5,
                    int gMaterialIndex = 0
  ) : AbstractTest(std::move(gInputDataDir),
                   gPrimaryEnergy,
                   gNumPrimaries,
                   gMaterialIndex
  ) {}
  void simulate() final;
  void writeHists();
};

class MollerTest : public AbstractTest {
 public:
  explicit MollerTest(std::string
                      gInputDataDir = "../data",
                      double gPrimaryEnergy = 12.345,
                      int gNumPrimaries = 1.0E+5,
                      int gMaterialIndex = 0
  ) : AbstractTest(std::move(gInputDataDir),
                   gPrimaryEnergy,
                   gNumPrimaries,
                   gMaterialIndex
  ) {}
  void simulate() final;
  void writeHists();
};

class MSCAngularDeflectionTest : public AbstractTest {
 public:
  explicit MSCAngularDeflectionTest(std::string
                                    gInputDataDir = "../data",
                                    double gPrimaryEnergy = 12.345,
                                    int gNumPrimaries = 1.0E+5,
                                    int gMaterialIndex = 0
  ) : AbstractTest(std::move(gInputDataDir),
                   gPrimaryEnergy,
                   gNumPrimaries,
                   gMaterialIndex
  ) {}
  void simulate() final;
  void writeHists();
};

#endif //DPM_MODELVALIDATIONAPI_TESTAPI_H_
