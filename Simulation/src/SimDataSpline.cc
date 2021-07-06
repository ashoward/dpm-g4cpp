#include "SimDataSpline.hh"

#include <iostream>
#include <cmath>

SimDataSpline::SimDataSpline(int size) : fNumData(size) { fData.resize(size); }


void SimDataSpline::FillData(int indx, double xval, double yval) {
  // just a check
  if (indx>=fNumData) {
    std::cerr << " *** ERROR SimDataSpline::FillData: \n"
              << "     indx = " << indx << " >= " << " fNumData = " << fNumData
              << std::endl;
    exit(EXIT_FAILURE);
  }
  //
  fData[indx].fX = xval;
  fData[indx].fY = yval;
  if (indx==0) {
    fMinX    = xval;
    fLogMinX = std::log(xval);
  }
  if (indx==1) {
    fInvLogDelta = 1./(std::log(xval)-fLogMinX);
  }
  if (indx==fNumData-1) {
    SetUp();
  }
}


void SimDataSpline::SetUp() {
  int    m1  = 1;
  int    m2  = fNumData - 1;
  double s   = 0.0;
  double r   = 0.0;
  for (int m=0; m<m2; ++m) {
    fData[m].fD =  fData[m+1].fX - fData[m].fX;
    r           = (fData[m+1].fY - fData[m].fY) / fData[m].fD;
    fData[m].fC = r - s;
    s           = r;
  }
  r = 0.0;
  s = 0.0;
  fData[0].fC  = 0.0;
  fData[m2].fC = 0.0;
  for (int m=m1; m<m2; ++m) {
    fData[m].fC = fData[m].fC + r * fData[m-1].fC;
    fData[m].fB = 2.0 * (fData[m-1].fX - fData[m+1].fX) - r*s;
    s           = fData[m].fD;
    r           = s/fData[m].fB;
  }
  int mr = m2 - 1;
  for (int m=m1; m<m2; ++m) {
    fData[mr].fC = (fData[mr].fD*fData[mr+1].fC - fData[mr].fC)/fData[mr].fB;
    --mr;
  }
  for (int m=0; m<m2; ++m) {
    s           = fData[m].fD;
    r           = fData[m+1].fC - fData[m].fC;
    fData[m].fD = r/s;
    fData[m].fC = 3.0*fData[m].fC;
    fData[m].fB = (fData[m+1].fY - fData[m].fY)/s - (fData[m].fC + r)*s;
//    fData[m].fA = fData[m].fY;
  }
}
