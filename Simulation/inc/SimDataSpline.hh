#ifndef SimDataSpline_HH
#define SimDataSpline_HH

//
// M. Novak: 2021
//
// Utility class to hold discretised function in a form of (x_i,y_i) like data
// pairs (that are used during the simulation e.g. mfp), and to provide y values
// interpolated at any x using SPLINE interpolation.
//
// NOTE: it is assumed that the discrete {x_i}_{i=0}^{N} values are equally
//       spaced on log-scale between x_0 and x_{N-1}. Therefore, a fast
//       calculation of the index `j`, such that x_j <= x < x_{j+1} is possible
//       given ln[x], ln[x_0] and the 1./ln[x_{i+1}/x_i] (that latter is const
//       for all `i`) as j = (int) (ln[x] - ln[x_0])/ln[x_{i+1}/x_i]

#include <vector>
#include <cmath>

class SimDataSpline {

public:
  SimDataSpline() {}
  SimDataSpline(int size);
 ~SimDataSpline() {}

  // set the size, i.e. the number of data to be stored (must be called before
  // filling in the data by using the FillData method)
  void  SetSize(int size) { fData.resize(size); fNumData = size; }

  // method provided to fill the data: supposed to be called fNumData times
  void  FillData(int indx, double xval, double yval);

  //
  // interpolation methods:
  //
  // ilow is the `j` index such that x_j <= x < x_{j+1}
  double GetValueAt(double xval, int ilow) const {
    xval -= fData[ilow].fX;
    return fData[ilow].fY + xval*(fData[ilow].fB + xval*(fData[ilow].fC + xval*fData[ilow].fD));
  }
  double GetValue(double xval, double logxval) const {
    int ilow = (int) ((logxval-fLogMinX)*fInvLogDelta);
    ilow = std::max(0, std::min(fNumData-2, ilow));
    return GetValueAt(xval, ilow);
  }
  double GetValue(double xval) const {
    const double logxval = std::log(xval);
    return GetValue(xval, logxval);
  }


private:

  // will be invoked automatically when the last data is filled in to compute
  // the spline interpolation parameter values
  void   SetUp();


private:

  // Number of data stored
  int     fNumData;
  // Minimum x-value
  double  fMinX;
  // The logarithm of the minimum x-value
  double  fLogMinX;
  // Inverse delta log x i.e. 1./ln[x_{i+1}/x_i]
  double  fInvLogDelta;

  // The x,y-values and the 4 spline interpolation parameters
  struct OnePoint {
    double fX;
    double fY;
//    double fA;  // we could get rid of it since it's fY
    double fB;
    double fC;
    double fD;
  };
  std::vector<OnePoint>  fData;
};

#endif // SimDataSpline_HH
