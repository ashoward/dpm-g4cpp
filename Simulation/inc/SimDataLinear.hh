#ifndef SimDataLinear_HH
#define SimDataLinear_HH

//
// M. Novak: 2021
//
// Utility class to hold discretised function in a form of (x_i,y_i) like data
// pairs (that are used during the simulation e.g. mfp), and to provide y values
// interpolated at any x using LINEAR interpolation.
//
// NOTE: it is assumed that the discrete {x_i}_{i=0}^{N} values are equally
//       spaced (on linear-scale) between x_0 and x_{N-1}. Therefore, a fast
//       calculation of the index `j`, such that x_j <= x < x_{j+1} is possible
//       given x, x_0 and the 1./[x_{i+1}-x_i] (that latter is const
//       for all `i`) as j = (int) (x - x_0)/[x_{i+1}-x_i]

#include <vector>
#include <cmath>
#include <iostream>
class SimDataLinear {

public:
  SimDataLinear() {}
  SimDataLinear(int size);
 ~SimDataLinear() {}

  // set the size, i.e. the number of data to be stored (must be called before
  // filling in the data by using the FillData method)
  void  SetSize(int size) { fData.resize(size); fNumData = size; }

  // method provided to fill the data: supposed to be called fNumData times
  void  FillData(int indx, double xval, double yval);

  //
  // interpolation methods:
  //
  // ilow is the `j` index such that x_j <= x < x_{j+1}
  double GetValueAt(double xval, int ilow) {
    return (fData[ilow+1].fY-fData[ilow].fY)*(xval-fData[ilow].fX)/(fData[ilow+1].fX-fData[ilow].fX) + fData[ilow].fY;
  }
  double GetValue(double xval) {
    int ilow = (int) ((xval-fMinX)*fInvDelta);
    ilow = std::max(0, std::min(fNumData-2, ilow));
    return GetValueAt(xval, ilow);
  }


private:

  // Number of data stored
  int     fNumData;
  // Minimum x-value
  double  fMinX;
  // Inverse delta x i.e. 1./[x_{i+1}-x_i]
  double  fInvDelta;

  // The x,y-values
  struct OnePoint {
    double fX;
    double fY;
  };
  std::vector<OnePoint>  fData;
};

#endif // SimDataLinear_HH
