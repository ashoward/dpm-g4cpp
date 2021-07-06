
#include "Geom.hh"

#include <cmath>
#include <cstdio>
#include <iostream>


Geom::Geom(double lbox, SimMaterialData* matData, int geomIndex)
: fPreDefinedGeomIndex(geomIndex),
  fLBox(lbox),
  fInvLBox(1./lbox),
  fLHalfBox(0.5*lbox),
  fMaterialData(matData) {
  fEdepHist.resize(20,0.0);
  fStepHist.resize(20,0.0);
}


double Geom::DistanceToBoundary(double* r, double* v, int* i) {
  static double kHalfTolerance = 0.5*kTolerance;
  //
  // Let's say that kExtent is the max extent of our geometry except the -z
  // direction in which it's only half box
  if (std::abs(r[0])>kExtent || std::abs(r[1])>kExtent || r[2]>kExtent || r[2]<-fLBox) {
    // indicates out of the geometry
    return -1.0;
  }
  // compute the current x,y and z box/volume index based on the current r(rx,ry,rz)
  // round downward to integer
  i[0] = std::floor((r[0]+fLHalfBox)*fInvLBox);
  i[1] = std::floor((r[1]+fLHalfBox)*fInvLBox);
  i[2] = std::floor((r[2]+fLHalfBox)*fInvLBox);
  //
  // transform the r(rx,ry,rz) into the current local box
  const double trX = r[0]-i[0]*fLBox;
  const double trY = r[1]-i[1]*fLBox;
  const double trZ = r[2]-i[2]*fLBox;
  //
  // compute the distance to boundary in the local box (centered around 0,0,0)
  double pdist = 0.0;
  double stmp  = 0.0;
  double snext = 1.0E+20;
  //
  // calculate x
  if (v[0] > 0.0) {
    pdist = fLHalfBox - trX;
    // check if actually this location is on boudnary
    if (pdist<kHalfTolerance) {
      // on boundary: push to the next box/volume, i.e. to the otherside and recompute
      r[0] += kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      snext = pdist/v[0];
    }
  } else if (v[0] < 0.0) {
    pdist = fLHalfBox + trX;
    if (pdist<kHalfTolerance) {
      // push to the otherside
      r[0] -= kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      snext = -pdist/v[0];
    }
  }
  //
  // calcualte y
  if (v[1] > 0.0) {
    pdist = fLHalfBox - trY;
    if (pdist<kHalfTolerance) {
      r[1] += kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      stmp = pdist/v[1];
      if (stmp < snext) {
        snext = stmp;
      }
    }
  } else if (v[1] < 0.0){
    pdist = fLHalfBox + trY;
    if (pdist<kHalfTolerance) {
      r[1] -= kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      stmp = -pdist/v[1];
      if (stmp < snext) {
        snext = stmp;
      }
    }
  }
  //
  // calculate z
  if (v[2] > 0.0) {
    pdist = fLHalfBox - trZ;
    if (pdist<kHalfTolerance) {
      r[2] += kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      stmp = pdist/v[2];
      if (stmp < snext) {
        snext = stmp;
      }
    }
  } else if (v[2] < 0.0) {
    pdist = fLHalfBox + trZ;
    if (pdist<kHalfTolerance) {
      r[2] -= kTolerance;
      return DistanceToBoundary(r, v, i);
    } else {
      stmp = -pdist/v[2];
      if (stmp < snext) {
        snext = stmp;
      }
    }
  }
  return snext;
}


void Geom::Score(double edep, int iz) {
  if (iz>=0) {
    const std::size_t indx = (std::size_t)(iz);
    if (indx>=fEdepHist.size()) {
      fEdepHist.resize(indx+10);
      fStepHist.resize(indx+10);
    }
    fEdepHist[indx] += edep;
    fStepHist[indx] += 1.0;
  }
}


int Geom::GetMaterialIndex(int* i) {
  const int iz = i[2];
  // vacuum
  if (iz < 0) {
    return -1;
  }
  // iz >= 0
  switch (fPreDefinedGeomIndex) {
    // homogeneous material with index 0
    case 0:
            return 0;
    // 0 - 1 cm --> 0; 1 - 3 cm --> 1; 3 - ... --> 0;
    case 1:
            if (iz < 10.0*fInvLBox) { return 0; } else
            if (iz < 30.0*fInvLBox) { return 1; }
            return 0;
    // 0 - 2 cm --> 0; 2 - 4 cm --> 1; 4 - ... --> 0;
    case 2:
            if (iz < 20.0*fInvLBox) { return 0; } else
            if (iz < 40.0*fInvLBox) { return 1; }
            return 0;
    // 0 - 1 cm --> 0; 1 - 1.5 cm --> 1; 1.5 - 2.5 cm --> 2; 2.5 - ... ->0
    case 3:
            if (iz < 10.0*fInvLBox) { return 0; } else
            if (iz < 15.0*fInvLBox) { return 1; } else
            if (iz < 25.0*fInvLBox) { return 2; }
            return 0;
    // default: homogeneous material
    default:
            return 0;
  }
}


void Geom::Write(const std::string& fname, int nprimaries) {
  FILE *f = fopen(fname.c_str(),"w");
  if (!f) {
    std::cerr<< " *** ERROR in Geom::Write(): "
             << " file " << fname << " could not be created! "
             << std::endl;
  }
  const double toCm = 0.1;
  int      sizeHist = fEdepHist.size();
  double       norm = 1./(nprimaries*fLBox*toCm);
  double    sumEdep = 0.0;
  for (int i=0; i<sizeHist; ++i) { sumEdep += fEdepHist[i]; }
  fprintf(f, "# === Mean energy deposit in the target: %13.4e  [MeV/event]\n", sumEdep/nprimaries);
  fprintf(f, "# === Energy deposit as a function of the depth: \n");
  // write the histogram: depth dose [MeV/cm]/density[g/cm3]--> [MeV cm2/g] as a function of depth [cm]
  fprintf(f, "# === Index       Depth [cm]    Edep [MeV cm2/g]  #Simulation-Steps \n");
  int idumy[] = {0,0,0};
  for (int i=0; i<sizeHist; ++i) {
    // set the iz voxel index and get the corresponding material density
    idumy[2] = i;
    double matDensity = GetVoxelMaterialDensity(idumy);
    fprintf(f, " %10d    %13.4e    %13.4e    %13.4e\n", i, (i+0.5)*fLBox*toCm, fEdepHist[i]*norm/matDensity, fStepHist[i]/nprimaries);
  }
  fclose(f);
  std::cout << " === Energy deposit histogram is written to the file:  " << fname << "\n" << std::endl;
}
