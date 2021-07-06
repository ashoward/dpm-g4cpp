#include "SimDataLinear.hh"

#include <iostream>
#include <cmath>

SimDataLinear::SimDataLinear(int size) : fNumData(size) { fData.resize(size); }


void SimDataLinear::FillData(int indx, double xval, double yval) {
  // just a check
  if (indx>=fNumData) {
    std::cerr << " *** ERROR SimDataLinear::FillData: \n"
              << "     indx = " << indx << " >= " << " fNumData = " << fNumData
              << std::endl;
    exit(EXIT_FAILURE);
  }
  //
  fData[indx].fX = xval;
  fData[indx].fY = yval;
  if (indx==0) {
    fMinX    = xval;
  }
  if (indx==1) {
    fInvDelta = 1./(xval-fMinX);
  }
}
