#ifndef PhotonData_HH
#define PhotonData_HH

#include <iostream>

#include <cmath>
#include <cstdio>
#include <string>

#include "KNTableBuilder.hh"

//#include "G4SystemOfUnits.hh"

//#include "Spline.hh"

// Data for photons per material including:
// - macroscopic cross sections for Compton, Pair-production and Photoelectric
//   NOTE: it is assumed that the photon absorption energy or photon production
//         threshold is set to higher than the highest absorption edge energy so
//         the Photoelectric cross section is a smooth function of the photon
//         energy
//  NOTE: we will use linear interpolation here instead of the spline and even
//        the common discrete photon energy grid will be linearly spaced so a
//        relatively high (~1024-2048) data will be used
class PhotonData {
  struct PhotonDataPerMaterial {
    // CTR that allocates memory
    PhotonDataPerMaterial(int numdata) : fNumData(numdata) {
      fMaterialName    = "";

      fComptonMXsec    = new double[numdata];
      fPairMXsec       = new double[numdata];
      fPhoElecMXsec    = new double[numdata];
      fTotalMXsec      = new double[numdata];
    }
    //
    ~PhotonDataPerMaterial() {
      delete [] fComptonMXsec;
      delete [] fPairMXsec;
      delete [] fPhoElecMXsec;
      delete [] fTotalMXsec;
    }

    // number of discrete energy grid point i.e. the size of the arrays below
    int         fNumData;

    std::string fMaterialName;
    double      fMaterialDensity; // in [g/cm3]

    // Macroscopic cross sections for Compton, Pair-production and Photoelectric
    // and their sum over the common kinetic energy grid
    double* fComptonMXsec;
    double* fPairMXsec;
    double* fPhoElecMXsec;
    double* fTotalMXsec;
  };

public:
  // emin = gammcut
  PhotonData(int nmaterial, double emin=0.0001, double emax=100, int negrid=1024) {
    fNumMaterial       = nmaterial;
    //
    fNumPhEData        = negrid;
    fInvDeltaPhE       = 1.0/( (emax-emin)/ (fNumPhEData-1) );
    // init the common energy grid
    fPhEGrid           = new double[negrid];
    fPhEGrid[0]        = emin;
    fPhEGrid[negrid-1] = emax;
    for (int i=1; i<negrid-1; ++i) {
      fPhEGrid[i] = emin + i/fInvDeltaPhE;
    }
    fGlobalMaxTotalMXsec = new double[negrid];
    // allocate data per material
    fDataPerMaterial = new PhotonDataPerMaterial*[nmaterial];
    for (int imat=0; imat<fNumMaterial; ++imat) {
      fDataPerMaterial[imat] = new PhotonDataPerMaterial(negrid);
    }
    // construct the Klein-Nishina table builder
    fKNTables = new KNTableBuilder();
  }

 ~PhotonData() {
    // free the dynamic kinetic energy grid
    delete[] fPhEGrid;
    delete[] fGlobalMaxTotalMXsec;
    // free all dynamic arrays allocated to store data per material
    for (int imat=0; imat<fNumMaterial; ++imat) {
      delete fDataPerMaterial[imat];
    }
    delete[] fDataPerMaterial;
  }


  double GetGlobalMaxIMFP(double ekin) {
    double mxsec = 0.0;
    if (ekin <= fPhEGrid[0]) {
      mxsec = fGlobalMaxTotalMXsec[0];
    } else if (ekin >= fPhEGrid[fNumPhEData-1]) {
      mxsec = fGlobalMaxTotalMXsec[fNumPhEData-1];
    } else {
      int eindx = (int)((ekin-fPhEGrid[0])*fInvDeltaPhE);
      mxsec = (fGlobalMaxTotalMXsec[eindx+1]-fGlobalMaxTotalMXsec[eindx])*fInvDeltaPhE*(ekin-fPhEGrid[eindx]) + fGlobalMaxTotalMXsec[eindx];
    }
    return std::max(1.0E-20, mxsec);
  }


  double GetTotalIMFP(double ekin, int imat) {
    if (imat<0) return 1.0E-20;
    double mxsec = 0.0;
    if (ekin <= fPhEGrid[0]) {
      mxsec = fDataPerMaterial[imat]->fTotalMXsec[0];
    } else if (ekin >= fPhEGrid[fNumPhEData-1]) {
      mxsec = fDataPerMaterial[imat]->fTotalMXsec[fNumPhEData-1];
    } else {
      int eindx = (int)((ekin-fPhEGrid[0])*fInvDeltaPhE);
      mxsec = (fDataPerMaterial[imat]->fTotalMXsec[eindx+1]-fDataPerMaterial[imat]->fTotalMXsec[eindx])*fInvDeltaPhE*(ekin-fPhEGrid[eindx]) + fDataPerMaterial[imat]->fTotalMXsec[eindx];
    }
    return std::max(1.0E-20, mxsec);
  }


  double GetComptonIMFP(double ekin, int imat) {
    if (imat<0) return 1.0E-20;
    double mxsec = 0.0;
    if (ekin <= fPhEGrid[0]) {
      mxsec = fDataPerMaterial[imat]->fComptonMXsec[0];
    } else if (ekin >= fPhEGrid[fNumPhEData-1]) {
      mxsec = fDataPerMaterial[imat]->fComptonMXsec[fNumPhEData-1];
    } else {
      int eindx = (int)((ekin-fPhEGrid[0])*fInvDeltaPhE);
      mxsec = (fDataPerMaterial[imat]->fComptonMXsec[eindx+1]-fDataPerMaterial[imat]->fComptonMXsec[eindx])*fInvDeltaPhE*(ekin-fPhEGrid[eindx]) + fDataPerMaterial[imat]->fComptonMXsec[eindx];
    }
    return std::max(1.0E-20, mxsec);
  }

  double GetPairIMFP(double ekin, int imat) {
    if (imat<0) return 1.0E-20;
    const double kEMC2 = 0.510991;
    double mxsec = 0.0;
    if (ekin < 2.0*kEMC2) {
      mxsec = 0.0;
    } else if (ekin <= fPhEGrid[0]) {
      mxsec = fDataPerMaterial[imat]->fPairMXsec[0];
    } else if (ekin >= fPhEGrid[fNumPhEData-1]) {
      mxsec = fDataPerMaterial[imat]->fPairMXsec[fNumPhEData-1];
    } else {
      int eindx = (int)((ekin-fPhEGrid[0])*fInvDeltaPhE);
      mxsec = (fDataPerMaterial[imat]->fPairMXsec[eindx+1]-fDataPerMaterial[imat]->fPairMXsec[eindx])*fInvDeltaPhE*(ekin-fPhEGrid[eindx]) + fDataPerMaterial[imat]->fPairMXsec[eindx];
    }
    return std::max(1.0E-20, mxsec);
  }

  double GetPhotoElecIMFP(double ekin, int imat) {
    if (imat<0) return 1.0E-20;
    double mxsec = 0.0;
    if (ekin <= fPhEGrid[0]) {
      mxsec = fDataPerMaterial[imat]->fPhoElecMXsec[0];
    } else if (ekin >= fPhEGrid[fNumPhEData-1]) {
      mxsec = fDataPerMaterial[imat]->fPhoElecMXsec[fNumPhEData-1];
    } else {
      int eindx = (int)((ekin-fPhEGrid[0])*fInvDeltaPhE);
      mxsec = (fDataPerMaterial[imat]->fPhoElecMXsec[eindx+1]-fDataPerMaterial[imat]->fPhoElecMXsec[eindx])*fInvDeltaPhE*(ekin-fPhEGrid[eindx]) + fDataPerMaterial[imat]->fPhoElecMXsec[eindx];
    }
    return std::max(1.0E-20, mxsec);
  }

  // Writes the data into the 'phys.dat' file (only those) that are required for
  // the simulation
  void WriteData(const std::string& fname) {
    char name[512];
    //
    // the imfp for compton scattering
    sprintf(name, "%s/imfp_compton.dat", fname.c_str());
    FILE* f = fopen(name, "w");
    fprintf(f, "# The imfp for Compton scattering\n");
    fprintf(f, "# E [MeV]   Compton_IMFP [1/mm] \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumPhEData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumPhEData; ++ie) {
        double ekin = fPhEGrid[ie];
        double imfp = fDataPerMaterial[im]->fComptonMXsec[ie];
        fprintf(f, "%.8E  %.8E\n", ekin, imfp/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the imfp for pair-production
    sprintf(name, "%s/imfp_pairp.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The imfp for e- e+ pair production.\n");
    fprintf(f, "# E [MeV]   Pair_IMFP [1/mm] \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumPhEData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumPhEData; ++ie) {
        double ekin = fPhEGrid[ie];
        double imfp = fDataPerMaterial[im]->fPairMXsec[ie];
        fprintf(f, "%.8E  %.8E\n", ekin, imfp/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the imfp for photoelectric effect
    sprintf(name, "%s/imfp_photoe.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The imfp for photoelectric effect.\n");
    fprintf(f, "# E [MeV]   PhotoEl_IMFP [1/mm] \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumPhEData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumPhEData; ++ie) {
        double ekin = fPhEGrid[ie];
        double imfp = fDataPerMaterial[im]->fPhoElecMXsec[ie];
        fprintf(f, "%.8E  %.8E\n", ekin, imfp/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the total imfp
    sprintf(name, "%s/imfp_total.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The total imfp = compton + pair + photoe  .\n");
    fprintf(f, "# E [MeV]   Total_IMFP [1/mm] \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumPhEData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumPhEData; ++ie) {
        double ekin = fPhEGrid[ie];
        double imfp = fDataPerMaterial[im]->fTotalMXsec[ie];
        fprintf(f, "%.8E  %.8E\n", ekin, imfp/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the total imfp
    sprintf(name, "%s/imfp_globalMax.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The global max tota-imfp.\n");
    fprintf(f, "# E [MeV]   Max-Total_IMFP [1/mm] \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumPhEData);
    fprintf(f, "# -----------------------------\n");
    for (int ie=0; ie<fNumPhEData; ++ie) {
      double ekin = fPhEGrid[ie];
      double imfp = fGlobalMaxTotalMXsec[ie];
      fprintf(f, "%.8E  %.8E\n", ekin, imfp);
    }
    fclose(f);

    // write the Klein-Nishina smapling tables
    fKNTables->Write(fname);
 }

  int     fNumMaterial;
  //
  double  fElectronCut; // in MeV
  double  fGammaCut;    // in MeV
  //
  int     fNumPhEData;
  double  fInvDeltaPhE;
  double* fPhEGrid;
  double* fGlobalMaxTotalMXsec;
  //
  PhotonDataPerMaterial** fDataPerMaterial;

  KNTableBuilder*  fKNTables;
};

#endif
