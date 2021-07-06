

#ifndef ElectronData_HH
#define ElectronData_HH

#include <iostream>

#include <cmath>
#include <cstdio>
#include <string>

#include "G4SystemOfUnits.hh"

#include "Spline.hh"
#include "GSDtrData.hh"
#include "SBTableBuilder.hh"
#include "MollerTableBuilder.hh"

// Data for electrons per material including:
// - elastic scattering related data (i.e. elastic and first transport MFP,
//   screening parameter that reproduce their ratio)
// - restricted stopping power (i.e. both radiative and collisonal) as well as
//   the corresponding range
class ElectronData {
  // elastic, stopping power and range data for a given material
  struct ElectronDataPerMaterial {
    // CTR that allocates memory
    ElectronDataPerMaterial(int numdata) : fNumData(numdata) {
      fMaterialName          = "";

      fMaterialDensity       = -1.0;
      fAtomicNumber          = -1.0;
      fInelCorrZs            = -1.0;
      fInelCorrZ2            = -1.0;
      fInelCorrFactor        = -1.0;
      fAtomicWeight          = -1.0;
      fZAweight              = -1.0;
      fMeanIoni              = -1.0;

//      fNumElements           = -1;
//      fAtomicNumbers         = nullptr;
//      fElemDensity           = nullptr;

      fElasticMFP            = new double[numdata];
      fFirstTransportMFP     = new double[numdata];
      fScreeningParamater    = new double[numdata];
      fMaxAllowedScatPower   = new double[numdata];

      fTheGSData             = nullptr;

      fStoppingPower         = new double[numdata];
      fRange                 = new double[numdata];

      fBremIMFPEkin          = new double[numdata];
      fBremIMFP              = new double[numdata];
      fMollerIMFPEkin        = new double[numdata];
      fMollerIMFP            = new double[numdata];

      fSpElasticMFP          = nullptr;
      fSpFirstTransportMFP   = nullptr;
      fSpScreeningParamater  = nullptr;
      fSpStoppingPower       = nullptr;
      fSpRange               = nullptr;
      fSpInverseRange        = nullptr;
      fSpMaxAllowedScatPower = nullptr;
      fSpBremIMFP            = nullptr;
      fSpMollerIMFP          = nullptr;
    }

    //
    ~ElectronDataPerMaterial() {
      // delete[] fAtomicNumbers;
      // delete[] fElemDensity;

      delete[] fElasticMFP;
      delete[] fFirstTransportMFP;
      delete[] fScreeningParamater;
      delete[] fMaxAllowedScatPower;

      delete   fTheGSData;

      delete[] fStoppingPower;
      delete[] fRange;

      delete[] fBremIMFPEkin;
      delete[] fBremIMFP;
      delete[] fMollerIMFPEkin;
      delete[] fMollerIMFP;

      delete fSpElasticMFP;
      delete fSpFirstTransportMFP;
      delete fSpScreeningParamater;
      delete fSpStoppingPower;
      delete fSpRange;
      delete fSpInverseRange;
      delete fSpMaxAllowedScatPower;
      delete fSpBremIMFP;
      delete fSpMollerIMFP;

    }

    // number of discrete energy grid point i.e. the size of the arrays below
    int fNumData;

    //
    // material properties
    std::string fMaterialName;

    double  fMaterialDensity; // in [g/cm3] units !!!!
    // n_i: fractional number of the i-th atom (#atoms_i_per_volume/total_#atoms_pre_volume
    double  fAtomicNumber;    // sum_i {n_i Z_i}        (atno in DPM)
    double  fInelCorrZ2;      // sum_i {n_i Z_i^2}      (atno2  in DPM)
    double  fInelCorrZs;      // sum_i {n_i Z_i(Z_i+1)} (atno21 in DPM)
    double  fInelCorrFactor;  // the ratio of the above two fInelCorrZs/fInelCorrZ2
    double  fAtomicWeight;    // sum_i {n_i w_i} with w_i in [amu] (mass in DPM)
    double  fZAweight;        // sum_i {n_i w_iZ_i} in [amu]  (zasum in DPM)
    double  fMeanIoni;        // exp[ sum_i {n_i ln(I_i)} ] with I_i [MeV] is the mean ionisation(excitation) energy (excite in DPM)

    // material composition
//    int     fNumElements;     // number of elements this material is built up
//    int*    fAtomicNumbers;   // atomic number of each elements
//    double* fElemDensities;   // number of atoms in unit volumes


    // elastic and first transport MFP over the discrete kinetic energy grid
    // as well as the screening paraneter values that reproduce their ratio
    double* fElasticMFP;
    double* fFirstTransportMFP;
    double* fScreeningParamater;
    double* fMaxAllowedScatPower;

    GSDtrData* fTheGSData;

    // the total (radiative and collisonal) restricted stopping power and the
    // corresponding range values over the discrete kinetic energy grid
    double* fStoppingPower;
    double* fRange;

    double* fBremIMFPEkin;
    double* fBremIMFP;
    double* fMollerIMFPEkin;
    double* fMollerIMFP;

    // the spline interpolators for all these above quantiteis
    Spline* fSpElasticMFP;
    Spline* fSpFirstTransportMFP;
    Spline* fSpScreeningParamater;

    Spline* fSpStoppingPower;
    Spline* fSpRange;
    Spline* fSpInverseRange;

    Spline* fSpMaxAllowedScatPower;

    Spline* fSpBremIMFP;
    Spline* fSpMollerIMFP;
  };


public:
  ElectronData(const std::string& datadir, int nmaterial=1, double emin=0.0001, double emax=100, int negrid=128) {
    fNumMaterial = nmaterial;

    // secondary electron (elCut) and gamma (gamCut) production thresholds in energy (MeV)
    // They will be set at init
    fElectronCut = -1.0;
    fGammaCut    = -1.0;

    fNumEKinData = negrid;
    fLogEKinMin  = std::log(emin);
    fLogEKinMax  = std::log(emax);
    fInvLogDelta = 1.0/( (fLogEKinMax-fLogEKinMin)/ (fNumEKinData-1) );
    // init the common energy grid
    fEKinGrid           = new double[negrid];
    fEKinGrid[0]        = emin;
    fEKinGrid[negrid-1] = emax;
    for (int i=1; i<negrid-1; ++i) {
      fEKinGrid[i] = std::exp(fLogEKinMin + i/fInvLogDelta);
    }

    // allocate data per material
    fDataPerMaterial = new ElectronDataPerMaterial*[nmaterial];
    for (int imat=0; imat<fNumMaterial; ++imat) {
      fDataPerMaterial[imat] = new ElectronDataPerMaterial(negrid);
    }
    //
    fSBTables     = new SBTableBuilder(datadir);
    fMollerTables = new MollerTableBuilder();
  }

 ~ElectronData() {
    // free the dynamic kinetic energy grid
    delete[] fEKinGrid;
    // free all dynamic arrays allocated to store data per material
    for (int imat=0; imat<fNumMaterial; ++imat) {
      delete fDataPerMaterial[imat];
    }
    delete[] fDataPerMaterial;
    //
    delete fSBTables;
    delete fMollerTables;
  }


  const std::string& GetMaterialName(int imat) { return fDataPerMaterial[imat]->fMaterialName; }

  // in g/cm3
  double  GetMaterialDensity(int imat) { return fDataPerMaterial[imat]->fMaterialDensity; }

  double GetMinEnergy() const { return fEKinGrid[0]; }
  double GetMaxEnergy() const { return fEKinGrid[fNumEKinData-1]; }

  //
  // Some quantities below, such as the elastic and first transport mean free
  // paths, can be obtained by including an inelastic scattering power correction
  // that is approximated by the sum_i {n_i Z_i(Z_i+1)} / sum_i { n_i Z_i^2}
  // factor, i.e. as replacing the usual Z^2 dependence of the cross sections
  // with the Z(Z+1).
  // This correction factor for the given material can be obtained by this method.
  double GetTotalScatteringCorrection(int imat) {
    return fDataPerMaterial[imat]->fInelCorrFactor;
  }

  // get interpolated elastic mfp in [mm] for the i-th matrial
  double GetElasticMFP(double ekin, int imat) {
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fElasticMFP[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fElasticMFP[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpElasticMFP->GetValueAt(ekin);
  }
  // inevrse ealstic mean free path with the inelastic scattering correction
  double GetInvTotalElasticMFP(double ekin, int imat) {
    return GetTotalScatteringCorrection(imat)/GetElasticMFP(ekin, imat);
  }

  // get interpolated imfp for bremsstrahlung in [1/mm] for the i-th matrial
  double GetBremIMFP(double ekin, int imat) {
    double mxsec = 0.0;
    if (ekin <= fDataPerMaterial[imat]->fBremIMFPEkin[0]) {
      mxsec = fDataPerMaterial[imat]->fBremIMFP[0];
    } else if (ekin >= fDataPerMaterial[imat]->fBremIMFPEkin[fNumEKinData-1]) {
      mxsec = fDataPerMaterial[imat]->fBremIMFP[fNumEKinData-1];
    } else {
      mxsec = fDataPerMaterial[imat]->fSpBremIMFP->GetValueAt(ekin);
    }
    return std::max(1.0E-20, mxsec);
  }
  double GetBremIMFPPerDensity(double ekin, int imat) {
    if (imat>-1) {
      return GetBremIMFP(ekin, imat)/fDataPerMaterial[imat]->fMaterialDensity;
    } else {
      return 1.0E-20;
    }
  }

  // get interpolated imfp for Moller in [mm] for the i-th matrial
  double GetMollerIMFP(double ekin, int imat) {
    double mxsec = 0.0;
    if (ekin <= fDataPerMaterial[imat]->fMollerIMFPEkin[0]) {
      mxsec = fDataPerMaterial[imat]->fMollerIMFP[0];
    } else if (ekin >= fDataPerMaterial[imat]->fMollerIMFPEkin[fNumEKinData-1]) {
      mxsec = fDataPerMaterial[imat]->fMollerIMFP[fNumEKinData-1];
    } else {
      mxsec = fDataPerMaterial[imat]->fSpMollerIMFP->GetValueAt(ekin);
    }
    return std::max(1.0E-20, mxsec);
  }
  double GetMollerIMFPPerDensity(double ekin, int imat) {
    if (imat>-1) {
      return GetMollerIMFP(ekin, imat)/fDataPerMaterial[imat]->fMaterialDensity;
    } else {
      return 1.0E-20;
    }
  }

  // get interpolated first transprt mfp in [mm] for the i-th matrial
  double GetFirstTransportMFP(double ekin, int imat) {
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fFirstTransportMFP[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fFirstTransportMFP[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpFirstTransportMFP->GetValueAt(ekin);
  }
  // inevrse 1rst transport mean free path with the inelastic scattering correction
  double GetInvTotal1rstTransportMFP(double ekin, int imat) {
    return GetTotalScatteringCorrection(imat)/GetFirstTransportMFP(ekin, imat);
  }
  double GetInvTotal1rstTransportMFPPerDensity(double ekin, int imat) {
    if (imat>-1) {
      return GetInvTotal1rstTransportMFP(ekin, imat)/fDataPerMaterial[imat]->fMaterialDensity;
    } else {
      return 1.0E-20;
    }
  }

  // get interpolated (screened Rutherford) screening parameter for the i-th matrial
  double GetScreeningPar(double ekin, int imat) {
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fScreeningParamater[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fScreeningParamater[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpScreeningParamater->GetValueAt(ekin);
  }

  // get interpolated (total) stopping power in [MeV/mm] for the i-th matrial
  double GetDEDX(double ekin, int imat) {
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fStoppingPower[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fStoppingPower[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpStoppingPower->GetValueAt(ekin);
  }
  double GetDEDXPerDensity(double ekin, int imat) {
    if (imat>-1) {
      return GetDEDX(ekin, imat)/fDataPerMaterial[imat]->fMaterialDensity;
    } else {
      return 1.0E-20;
    }
  }

  // get interpolated (total) range in [mm] for the i-th matrial
  double GetRange(double ekin, int imat) {
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fRange[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fRange[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpRange->GetValueAt(ekin);
  }
  // get interpolated energy after a given path in [MeV] for the i-th matrial
  double GetInverseRange(double path, int imat) {
    if (path <= fDataPerMaterial[imat]->fRange[0]) {
      return fEKinGrid[0];
    }
    if (path >= fDataPerMaterial[imat]->fRange[fNumEKinData-1]) {
      return fEKinGrid[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpInverseRange->GetValueAt(path);
  }

  // get maximaly allowed scateering strength along the step (according to the
  // stepping parameters s-low, s-high and e-cross)
  double GetMaxAllowedScatLength(double ekin, int imat) {
    if (imat<0) return 1.0E+20;
    if (ekin <= fEKinGrid[0]) {
      return fDataPerMaterial[imat]->fMaxAllowedScatPower[0];
    }
    if (ekin >= fEKinGrid[fNumEKinData-1]) {
      return fDataPerMaterial[imat]->fMaxAllowedScatPower[fNumEKinData-1];
    }
    return fDataPerMaterial[imat]->fSpMaxAllowedScatPower->GetValueAt(ekin);
  }
  double GetMaxAllowedScatLengthPerDensity(double ekin, int imat) {
    if (imat>-1) {
      return GetMaxAllowedScatLength(ekin, imat)/fDataPerMaterial[imat]->fMaterialDensity;
    } else {
      return 1.0E+20;
    }
  }

  // same as above but computes the value instead the above interpolation
  // (see more on the details at `InitElectronData::InitScatteringData()`)
  double ComputeMaxAllowedScatPow(double ekin, int imat, double parShigh, double parSlow, double parEcross) {
    const double kInvPi =  1.0/3.1415926535897932;
    const double kSharp =  5.0;
    const double kinE = std::min(std::max(fEKinGrid[0], ekin), fEKinGrid[fNumEKinData-1]-1.0E-8);
    const double  fx  = std::atan(kSharp*std::log(kinE/parEcross))*kInvPi + 0.5;
    return GetInvTotal1rstTransportMFP(kinE, imat)*parSlow*std::pow(parShigh/parSlow,fx);
  }

  // Provide sample from the GS angular distribution specified by the material
  // index and the kinetic energy.
  // NOTE: in the current version, material and kinetic energy determines the
  //       all required quantities such as s/lambda_e, screening parameter and
  //       pdf transformation parameter. This is because the stepping parameters
  //       (sLow, sHigh, eCross) determines the maximally allowed step length
  //       at a given kinetic energy in a given material (through the maximally
  //       allowed scattering strength as computed in InitElectronData::InitScatteringData
  //       or above in ComputeMaxAllowedScatPow). This maximal step length is
  //       considered at kinetic energy (and material) when the GS angular
  //       distributions are pre-computed (in InitElectronData::InitGSData).
  //       So during the simulation, the MSC angualr defelection will always be
  //       given to this maximum step that corresponds to the given kinetic energy
  //       in the given material.
  //
  // ekin must be between [electron production threshold, emax ] otherwise
  // this limits are used.
  double SampleFromGSDtr(double ekin, int imat, double rndm1, double rndm2) {
    // find the discrete kinetic energy grid point to take
    const GSDtrData* gsData = fDataPerMaterial[imat]->fTheGSData;
    ekin = std::min(gsData->fEKinGrid[gsData->fNumEKin-1]-1.0E-7, std::max(gsData->fEKinGrid[0], ekin) );
    const double lekin = std::log(ekin);
    double pEip1    = (lekin-gsData->fLogEKinMin)*gsData->fInvLogDelta;
    // lower index of the kinetic energy bin
    const int ielow = (int)pEip1;
    int    iekin = ielow;
    // prob. of taking the higher kinetic energy grid
    pEip1       -= iekin;
    if (rndm1 < pEip1) {
      ++iekin;
    }
    const GSDtrData::OneGSDtr& theGSTable = gsData->fGSDtrData[iekin];
    // lower index of the (common) discrete cumulative bin and the residual fraction
    const double delta = 1./gsData->fInvDeltaCum;
    const int    indxl = (int)(rndm2/delta);
    const double resid = rndm2-indxl*delta;

    // compute value `u` by using ratin based numerical inversion
    const double  parA = theGSTable.fCumData[indxl].fParmA;
    const double  parB = theGSTable.fCumData[indxl].fParmB;
    const double    u0 = theGSTable.fCumData[indxl].fVarU;
    const double    u1 = theGSTable.fCumData[indxl+1].fVarU;
    const double  dum1 = (1.0 + parA + parB) * delta * resid;
    const double  dum2 = delta * delta + parA * delta * resid + parB * resid * resid;
    const double  theU = u0 + dum1 / dum2 * (u1 - u0);
    // transform back the sampled `u` to `mu(u)` using the transformation parameter `a`
    // mu(u) = 1 - 2au/[1-u+a] as given by Eq.(34)
    // interpolate (linearly) the transformation parameter to E
    const double a0 = gsData->fGSDtrData[ielow].fTransformParam;
    const double a1 = gsData->fGSDtrData[ielow+1].fTransformParam;
    const double e0 = gsData->fEKinGrid[ielow];
    const double e1 = gsData->fEKinGrid[ielow+1];
    const double parTransf = (a1-a0)/(e1-e0)*(ekin-e0)+a0;
    return 1.-2.*parTransf*theU/(1.-theU+parTransf);
  }

  // Writes the data into the 'phys.dat' file (only those) that are required for
  // the simulation
  void WriteData(const std::string& fname) {
    char name[512];
    // reference material is assumed to be the very first
    int imref = 0;
    // material related informations
    sprintf(name, "%s/mat.dat", fname.c_str());
//    std::cout << name << std::endl;
    FILE *f = fopen(name, "w");
    //
    if (!f) {
      std::cout << "\n *** ERROR in ElectronData::Write(): "
                << "\n     The file = " << name << " cannot be created!"
                << "\n     (The directory <" << fname << "> probability dosen't exist)!\n"
                << std::endl;
      exit(-1);
    }
    // write first the electron and gamma prodcution thresholds in [MeV]
    fprintf(f, "# Secondary electron and gamma production thresholds:\n %.8E %.8E\n",fElectronCut, fGammaCut);
    // write MSC maximum step parameters 'slow', 'shigh' and 'ecross'
    fprintf(f, "# MSC maximum step length parameters:\n %.8E %.8E %.8E\n",fMSCStepParSLow, fMSCStepParSHigh, fMSCStepParEcross);
    // number of materials
    fprintf(f, "# Number of materials:\n %d\n",fNumMaterial);
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "# \n");
      fprintf(f, "# Material index, name and density [g/cm3] data:\n");
      fprintf(f, " %d %16s %.8E\n", im,  fDataPerMaterial[im]->fMaterialName.c_str(), fDataPerMaterial[im]->fMaterialDensity);
      fprintf(f, "# Z, Z^2, Z(Z+1), w [amu], wZ [amu], I [MeV] (all are weighted for molecules)\n");
      fprintf(f, " %.8E  %.8E  %.8E  %.8E  %.8E  %.8E\n",
                   fDataPerMaterial[im]->fAtomicNumber,
                   fDataPerMaterial[im]->fInelCorrZ2,
                   fDataPerMaterial[im]->fInelCorrZs,
                   fDataPerMaterial[im]->fAtomicWeight,
                   fDataPerMaterial[im]->fZAweight,
                   fDataPerMaterial[im]->fMeanIoni
              );
      fprintf(f, "# [Z/w]/[Z/w]_ref, [Z^2/w]/[Z^2/w]_ref (all are weighted for molecules) ref: im=0\n");
      fprintf(f, " %.8E  %.8E\n",
                   fDataPerMaterial[im]->fAtomicNumber/fDataPerMaterial[im]->fAtomicWeight*fDataPerMaterial[imref]->fAtomicWeight/fDataPerMaterial[imref]->fAtomicNumber,
                   fDataPerMaterial[im]->fInelCorrZ2/fDataPerMaterial[im]->fAtomicWeight*fDataPerMaterial[imref]->fAtomicWeight/fDataPerMaterial[imref]->fInelCorrZ2
             );
    }
    fclose(f);
    //
    // the inverse elastic mean free path data (i.e. macroscopic cross section)
    sprintf(name, "%s/el_iemfp.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# Inverse elastic mean free path, i.e. macroscopic cross section, data\n");
    fprintf(f, "# The `_scpc` are the scattering power corrected, i.e. inelastic scattering correected, data) \n");
    fprintf(f, "# E [MeV]   1/Lambda_elastic [1/mm] 1/Lambda_elastic_scpc [1/mm]  \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumEKinData);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double emfp = fDataPerMaterial[im]->fElasticMFP[ie];
        fprintf(f, "%.8E  %.8E  %.8E\n", fEKinGrid[ie], 1.0/emfp, GetTotalScatteringCorrection(im)/emfp);
      }
    }
    fclose(f);
    //
    // the inverse first transport mean free path data (i.e. macroscopic cross section)
    sprintf(name, "%s/el_itr1mfp.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# Inverse 1-rst transport mean free path, i.e. macroscopic cross section, data\n");
    fprintf(f, "# The `scpc` are the scattering power corrected, i.e. inelastic scattering correected, data) \n");
    fprintf(f, "# E [MeV]   1/Lambda_1 [1/mm][cm3/g] 1/Lambda_1^{scpc} [1/mm][cm3/g]  \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumEKinData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double tr1mfp = fDataPerMaterial[im]->fFirstTransportMFP[ie];
        fprintf(f, "%.8E  %.8E  %.8E\n", fEKinGrid[ie], 1.0/tr1mfp/matDensityIn_gPerCm3, GetTotalScatteringCorrection(im)/tr1mfp/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the screening parameter of the screened Rutherford DCS
    // (i.e. \eta(E_i), such that the accurate G_1(\eta) = mfp_el/mfp_tr1)
    sprintf(name, "%s/el_screeningPar.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# Screening parameter of the screened Rutherford DCS\n");
    fprintf(f, "# (i.e. \\eta(E_i), such that the accurate G_1(\\eta) = mfp_el/mfp_tr1) \n");
    fprintf(f, "# E [MeV]   Screening_parameter [] \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumEKinData);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double eta = fDataPerMaterial[im]->fScreeningParamater[ie];
        fprintf(f, "%.8E  %.8E\n", fEKinGrid[ie], eta);
      }
    }
    fclose(f);
    //
    // the max scattering strength allowed in a given step
    sprintf(name, "%s/el_scatStrength.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The max scattering strength `K_1(E)` allowed in a given MSC step.\n");
    fprintf(f, "# Max K_1(E) is approximated as max-MSC-step-length(E)/tr1-mfp(E).\n");
    fprintf(f, "# E [MeV]   Max-scattering-strength [] \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumEKinData);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double maxScatStr = fDataPerMaterial[im]->fMaxAllowedScatPower[ie];
        fprintf(f, "%.8E  %.8E\n", fEKinGrid[ie], maxScatStr);
      }
    }
    fclose(f);
    //
    // Data related to the transformed GS angular distribution sampling table
    sprintf(name, "%s/el_GSDtrData.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# Data related to the transformed GS angular distributions:\n");
    fprintf(f, "# (optimal transformation parameters `a` and tables for fast,\n");
    fprintf(f, "#  rejection free sampling from the q(u;s/emfp,scrPar,a) pdf-s)\n");
    fprintf(f, "# #ekin, emin [MeV], emax [MeV] and #cumulative:\n");
    int gsNumEkin = fDataPerMaterial[0]->fTheGSData->fNumEKin;
    int gsNumCum  = fDataPerMaterial[0]->fTheGSData->fNumCumData;
    fprintf(f, "%d  %.14E  %.14E  %d\n",
               gsNumEkin,
               fDataPerMaterial[0]->fTheGSData->fEKinGrid[0],
               fDataPerMaterial[0]->fTheGSData->fEKinGrid[gsNumEkin-1],
               gsNumCum
            );
    fprintf(f, "# For each materials: #ekin tables and 3x#cumulative data in each table,\n");
    fprintf(f, "# i, E_i [MeV] and transformation parameter are given before each table.\n");
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      const std::vector<double>& vekin = fDataPerMaterial[im]->fTheGSData->fEKinGrid;
      for (int ie=0; ie<gsNumEkin; ++ie) {
        const GSDtrData::OneGSDtr& gsTable = fDataPerMaterial[im]->fTheGSData->fGSDtrData[ie];
        // ekin, transformation parameter
        fprintf(f, "%d  %.14E  %.14E\n", ie, vekin[ie], gsTable.fTransformParam);
        for (int iu=0; iu<gsNumCum; ++iu) {
          fprintf(f, "%22.14E  %22.14E  %22.14E ", gsTable.fCumData[iu].fVarU, gsTable.fCumData[iu].fParmA, gsTable.fCumData[iu].fParmB);
          if ( (iu+1)%2==0 || iu==gsNumCum-1) {
            fprintf(f, "\n");
          }
        }
      }
    }
    fclose(f);
    //
    // the restricted, total (collisonal and radiative) stopping power
    sprintf(name, "%s/eloss_rdedx.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The restricted, total (collisonal and radiative) stopping power");
    fprintf(f, "# E [MeV]   dE/dx [MeV/mm][cm3/g] \n");
    fprintf(f, "# number of data per materials and #material\n %d %d\n", fNumEKinData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double dedx = fDataPerMaterial[im]->fStoppingPower[ie];
        fprintf(f, "%.8E  %.8E\n", fEKinGrid[ie], dedx/matDensityIn_gPerCm3);
      }
    }
    fclose(f);
    //
    // the restricted, total (collisonal and radiative) range
    sprintf(name, "%s/eloss_rRange.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The restricted, total (collisonal and radiative) range");
    fprintf(f, "# E [MeV]   Range [mm] \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumEKinData);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double range = fDataPerMaterial[im]->fRange[ie];
        fprintf(f, "%.8E  %.8E\n", fEKinGrid[ie], range);
      }
    }
    fclose(f);
    //
    // the restricted mfp for bremsstrahlung
    sprintf(name, "%s/imfp_brem.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The restricted imfp for bremsstrahlung interaction\n");
    fprintf(f, "# E [MeV]   Brem_IMFP [1/mm][cm3/g] \n");
    fprintf(f, "# number of data per materials and #materials\n %d %d\n", fNumEKinData, fNumMaterial);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double ekin = fDataPerMaterial[im]->fBremIMFPEkin[ie];
        double imfp = fDataPerMaterial[im]->fBremIMFP[ie]/matDensityIn_gPerCm3;
        fprintf(f, "%.8E  %.8E\n", ekin, imfp);
      }
    }
    fclose(f);
    //
    // the restricted imfp for Moller (ionisation) interaction
    sprintf(name, "%s/imfp_moller.dat", fname.c_str());
    f = fopen(name, "w");
    fprintf(f, "# The restricted imfp for Moller interaction\n");
    fprintf(f, "# E [MeV]   Moller_IMFP [1/mm][cm3/g] \n");
    fprintf(f, "# number of data per materials\n %d\n", fNumEKinData);
    fprintf(f, "# -----------------------------\n");
    for (int im=0; im<fNumMaterial; ++im) {
      fprintf(f, "#\n# Material name: %s\n#\n", fDataPerMaterial[im]->fMaterialName.c_str());
      double matDensityIn_gPerCm3 = fDataPerMaterial[im]->fMaterialDensity;
      for (int ie=0; ie<fNumEKinData; ++ie) {
        double ekin = fDataPerMaterial[im]->fMollerIMFPEkin[ie];
        double imfp = fDataPerMaterial[im]->fMollerIMFP[ie]/matDensityIn_gPerCm3;
        fprintf(f, "%.8E  %.8E\n", ekin, imfp);
      }
    }
    fclose(f);

    //
    // write sampling tables data for rejection free Moller interaction
    fMollerTables->Write(fname);
    //
    // write sampling tables data for rejection free (Seltzer-Berger) brem. interaction
    fSBTables->Write(fname);

  }


  //
  // data members (all public for simplicity)
  //
  int      fNumMaterial;
  // secondary e- and gamma prodcution thresholds in [MeV]
  double   fElectronCut; // in MeV
  double   fGammaCut;    // in MeV
  // MSC maximum step length paraeters
  double   fMSCStepParSHigh;
  double   fMSCStepParSLow;
  double   fMSCStepParEcross;

  // this energy grid is used for the elastic, stopping power and range data
  int      fNumEKinData;
  double   fLogEKinMin;
  double   fLogEKinMax;
  double   fInvLogDelta;
  double*  fEKinGrid;

  // physics data for the individiual materials
  ElectronDataPerMaterial** fDataPerMaterial;
  SBTableBuilder*           fSBTables;
  MollerTableBuilder*       fMollerTables;
};


#endif // DATA_HH
