#include "SimLinAliasData.hh"

#include <iostream>
#include <cmath>

SimLinAliasData::SimLinAliasData(int size) : fNumData(size) { fTheTable.resize(size); }
SimLinAliasData::~SimLinAliasData() { fTheTable.clear(); }


void SimLinAliasData::FillData(int indx, double xval, double yval, double aliasw, int aliasi) {
  // just a check
  if (indx>=fNumData) {
    std::cerr << " *** ERROR SimLinAliasData::FillData: \n"
              << "     indx = " << indx << " >= " << " fNumData = " << fNumData
              << std::endl;
    exit(EXIT_FAILURE);
  }
  fTheTable[indx].fXdata     = xval;
  fTheTable[indx].fYdata     = yval;
  fTheTable[indx].fAliasW    = aliasw;
  fTheTable[indx].fAliasIndx = aliasi;
}

// produce sample from the represented distribution u
double SimLinAliasData::Sample(double rndm1, double rndm2) {
  // get the lower index of the bin by using the alias part
  double rest  = rndm1*(fNumData-1);
  int    indxl = (int) (rest);
  if (fTheTable[indxl].fAliasW<rest-indxl)
    indxl = fTheTable[indxl].fAliasIndx;
  // sample value within the selected bin by using linear aprox. of the p.d.f.
  double xval   = fTheTable[indxl].fXdata;
  double xdelta = fTheTable[indxl+1].fXdata-xval;
  if (fTheTable[indxl].fYdata > 0.0) {
    double dum = (fTheTable[indxl+1].fYdata-fTheTable[indxl].fYdata)/fTheTable[indxl].fYdata;
    if (std::abs(dum)>0.1)
      return xval - xdelta/dum * (1.0 - std::sqrt(1.0+rndm2*dum*(dum+2.0)));
    else // use second order Taylor around dum = 0.0
      return xval + rndm2*xdelta*(1.0-0.5*dum*(rndm2-1.0)*(1.0+dum*rndm2));
  }
  return xval + xdelta*std::sqrt(rndm2);
}
