#ifndef SimLinAliasData_HH
#define SimLinAliasData_HH

//
// M. Novak: 2021
//
// Sampling table based on Walker's alias sampling and linear approximation of
// the discretized p.d.f.
// It is assumed that the data will be filled from file through the provided
// `FillData(int indx)` method. The `Sample()` method can then be used for
// producing samples form the represented (approximated) distribution without
// rejection.

#include <vector>

class SimLinAliasData {
public:
  // CTR and DTR
  SimLinAliasData(int size);
 ~SimLinAliasData();

  // method provided to fill the data: supposed to be called fNumData times
  void   FillData(int indx, double xval, double yval, double aliasw, int aliasi);

  // method provided to produce samples according to the represented distribution
  // `rndm1` and `rndm2` are uniformly random values on [0,1]
  double Sample(double rndm1, double rndm2);


private:
  // data that describes a single point in the table
  struct OnePoint {
    // the discretized stochastic variable value
    double fXdata;
    // the corresponding p.d.f. value (not necessarily normalised over fXdata)
    double fYdata;
    // the correspnding alias probability (not necessarily normalised over fXdata)
    double fAliasW;
    // the corresponding alias index
    int    fAliasIndx;
  };

  //
  // Data member declaration:
  //
  // #data points in the table
  int                    fNumData;
  // the fNumData points that build up the sampling table
  std::vector<OnePoint>  fTheTable;  // [fNumData]


};

#endif
