
#include "InitElectronData.hh"

#include "ElectronData.hh"
#include "Spline.hh"
#include "GLIntegral.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4ProductionCutsTable.hh"
#include "G4DataVector.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"

#include "G4eDPWAElasticDCS.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"

#include "SBTableBuilder.hh"

#include "GS.hh"
#include "AliasTable.hh"
#include "Spline.hh"

#include <algorithm>

void   InitElectronData(ElectronData& elData,
                        double electronCutInEnergy, double gammaCutInEnergy,
                        double parShigh, double parSlow, double parEcross) {
  // set secondary e-/gamma prodcution thresholds in MeV
  elData.fElectronCut = electronCutInEnergy;
  elData.fGammaCut    = gammaCutInEnergy;
  //
  elData.fMSCStepParSHigh  = parShigh;
  elData.fMSCStepParSLow   = parSlow;
  elData.fMSCStepParEcross = parEcross;
  //
  elData.fSBTables->BuildTables(elData.fNumMaterial, elData.fGammaCut, elData.GetMaxEnergy());
  elData.fMollerTables->BuildTables(elData.fElectronCut, elData.GetMaxEnergy());

//  std::cout << " === Electron related data init... " << std::endl;
  InitElasticData(elData);
  InitElossData(elData, electronCutInEnergy, gammaCutInEnergy);
  InitScatteringData(elData, parShigh, parSlow, parEcross);

  InitGSData(elData, electronCutInEnergy);
}


void   InitElasticData(ElectronData& elData) {
  // create G4eDPWAElasticDCS for electron
  G4eDPWAElasticDCS theDPWADCS;
  // get the list of materials from
  int numMat = elData.fNumMaterial;
//  std::cout << "    --- Elastic scattering related data init " << std::endl;
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
    elData.fDataPerMaterial[im]->fMaterialName    = mat->GetName();
    elData.fDataPerMaterial[im]->fMaterialDensity = mat->GetDensity()/(CLHEP::g/CLHEP::cm3);
//    std::cout << " material: " << mat->GetName() << std::endl;
    // get the element composition and:
    // - compute the macroscopic elastic quantities
    // - in each case of Z, make sure that the theDPWADCS is initialised for that
    const double* elemDensity       = mat->GetVecNbOfAtomsPerVolume();
    const G4ElementVector* elemVect = mat->GetElementVector();
    const int numElems = elemVect->size();
    // for each discrete kinetic energy grid point in the Electron data, compute
    // the elastic data and fill into the Electron data obejct for this material
    double* theEKinGrid = elData.fEKinGrid;
    int     numEKin     = elData.fNumEKinData;
    for (int ie=0; ie<numEKin; ++ie) {
      // kinetic energy
      double ekin     = theEKinGrid[ie];
      // some inelastic correction data cumulatives
      double theZ     = 0.0; // sum_i ni Z_i
      double theZs    = 0.0; // sum_i ni Z_i (Z_i+1)
      double theZ2    = 0.0; // sum_i ni Z_i^2
      double theW     = 0.0; // sum_i ni w_i
      double theZW    = 0.0; // sum_i ni w_iZ_i
      double theIexc  = 0.0; // exp (sum_i ni ln(I_i))
      // some ealstic data cumulatives
      double elMXsec  = 0.0;
      double tr1MXsec = 0.0;
      double tr2MXsec = 0.0;
      for (int iel=0; iel<numElems; ++iel) {
        const G4Element *elem = (*elemVect)[iel];
        double zet = elem->GetZ();
        int   izet = std::lrint(zet);
        if (ie==0) {
          theDPWADCS.InitialiseForZ(izet);
          theZ    += elemDensity[iel]*zet;
          theZs   += elemDensity[iel]*zet*(zet+1.0);
          theZ2   += elemDensity[iel]*zet*zet;
          theW    += elemDensity[iel]*elem->GetAtomicMassAmu();
          theZW   += elemDensity[iel]*elem->GetAtomicMassAmu()*zet;
          theIexc += elemDensity[iel]*std::log(elem->GetIonisation()->GetMeanExcitationEnergy());
        }
        double elXsec, tr1Xsec, tr2Xsec;
        theDPWADCS.ComputeCSPerAtom(izet, ekin, elXsec, tr1Xsec, tr2Xsec);
        elMXsec  += elemDensity[iel]*elXsec;
        tr1MXsec += elemDensity[iel]*tr1Xsec;
        tr2MXsec += elemDensity[iel]*tr2Xsec;
      }
      // write the inelastic correction data for this material
      if (ie==0) {
        const double totalNumberOfAtoms = mat->GetTotNbOfAtomsPerVolume();
        elData.fDataPerMaterial[im]->fAtomicNumber   = theZ/totalNumberOfAtoms;
        elData.fDataPerMaterial[im]->fInelCorrZs     = theZs/totalNumberOfAtoms;
        elData.fDataPerMaterial[im]->fInelCorrZ2     = theZ2/totalNumberOfAtoms;
        elData.fDataPerMaterial[im]->fInelCorrFactor = theZs/theZ2;
        elData.fDataPerMaterial[im]->fAtomicWeight   = theW/totalNumberOfAtoms;
        elData.fDataPerMaterial[im]->fZAweight       = theZW/totalNumberOfAtoms;
        elData.fDataPerMaterial[im]->fMeanIoni       = std::exp(theIexc/totalNumberOfAtoms);
      }
      // write the elastic and first transport MFP into the ElectronData
      elData.fDataPerMaterial[im]->fElasticMFP[ie]        = 1.0/elMXsec;
      elData.fDataPerMaterial[im]->fFirstTransportMFP[ie] = 1.0/tr1MXsec;
      // solve the mfp_el/mfp_tr1 = G1(x) equation for x, i.e. for the screening
      // parameter such that G1(x) = 2x[ln(1+1/x)(1+x)-1] reproduce this ratio
      double screeningPar = FindScreeingParameter(tr1MXsec/elMXsec);
      elData.fDataPerMaterial[im]->fScreeningParamater[ie] = screeningPar;
    }
    // set up the Spline interpolator for these 3 elastic quantities
    elData.fDataPerMaterial[im]->fSpElasticMFP         = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fElasticMFP, numEKin);
    elData.fDataPerMaterial[im]->fSpFirstTransportMFP  = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fFirstTransportMFP, numEKin);
    elData.fDataPerMaterial[im]->fSpScreeningParamater = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fScreeningParamater, numEKin);
  }
}

// val = G1_ref where G1(x) = 2x[ln(1+1/x)(1+x)-1]
double FindScreeingParameter(double val) {
  const double eps = 1.0E-10;
  // incorporate the factor 2 from the right
  double  half  = 0.5*val;
  double ihalf  = 1.0/half;
  // find appropriate initial value for x
  double x = half;
  while ( x*((1.0+x)*std::log(1.0+1.0/x)-1.0) > half ) {
    x *= 0.1;
  }
  // Newton's method to find the root
  while (true) {
    double lterm = std::log(1.0+1.0/x);
    double g1    = x*((1.0+x)*lterm - 1.0);
    if (std::abs(g1*ihalf - 1.0) < eps) {
      break;
    }
    g1 -= half;
    double del = -2.0 + (1.0+2.0*x)*lterm;
    x  -= g1/del;
  }
  return x;
}


void   InitElossData(ElectronData& elData, double electronCutInEnergy, double gammaCutInEnergy) {
  // create Penelope models for brem and inoni for electrons to get the total SP
  G4PenelopeIonisationModel theIoni;
  G4PenelopeBremsstrahlungModel theBrem;
  theIoni.SetHighEnergyLimit(CLHEP::GeV);
  theBrem.SetHighEnergyLimit(CLHEP::GeV);

  G4ParticleDefinition* part = G4Electron::Electron();

  // get cuts for secondary e-
  const G4DataVector* theElCuts = static_cast<const G4DataVector*>(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(1));
  theIoni.Initialise(part, *theElCuts);
  // get cuts for secondary gamma
  const G4DataVector* theGamCuts = static_cast<const G4DataVector*>(G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(0));
  theBrem.Initialise(part, *theGamCuts);

  // create GL integral on 0,1 and a spline util
  int ngl = 16;
  GLIntegral* gl = new GLIntegral(ngl, 0.0, 1.0);
  Spline      sp;

  // get the list of materials from
  int numMat = elData.fNumMaterial;
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
//    std::cout << " material: " << mat->GetName() << std::endl;
    // generate a kinetic energy grid starting from a low (100 eV) for the dEdx
    // and range computations
    int    nekins = 1024;
    double  ekin0 = 0.0001; // 100 eV
    double  ekin1 = elData.fEKinGrid[elData.fNumEKinData-1];
    double lekin0 = std::log(ekin0);
    double ldekin = std::log(ekin1/ekin0)/(nekins-1);
    std::vector<double> vekins(nekins);
    std::vector<double> vdedx(nekins);
    std::vector<double> vrange(nekins);
    vekins[0]        = ekin0;
    vekins[nekins-1] = ekin1;
    for (int je = 1; je<nekins-1; ++je) {
      vekins[je] = std::exp(lekin0 + je*ldekin);
    }
    // compute the restricted dEdx (collisonal and radiative) over this ekin grid
    for (int ie=0; ie<nekins; ++ie) {
      double  ekin = vekins[ie];
      double cDEDX = theIoni.ComputeDEDXPerVolume(mat, part, ekin, electronCutInEnergy);
      double rDEDX = theBrem.ComputeDEDXPerVolume(mat, part, ekin, gammaCutInEnergy);
      vdedx[ie] = cDEDX + rDEDX;
    }
    // set a Spline util for this stopping power
    sp.SetUpSpline(vekins.data(), vdedx.data(), nekins);
    // use this spline to fill in the restr. dEdx data (that has a different ekin grid)
    double* theEKinGrid = elData.fEKinGrid;
    int     numEKin     = elData.fNumEKinData;
    for (int ie=0; ie<numEKin; ++ie) {
      double ekin = theEKinGrid[ie];
      elData.fDataPerMaterial[im]->fStoppingPower[ie] = sp.GetValueAt(ekin);
    }
    // set the Spline for this stopping power data
    elData.fDataPerMaterial[im]->fSpStoppingPower = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fStoppingPower, numEKin);
    //
    //now integrate the computed 1/dedx to get the range: use range = 0 at the lowest energy!
    const std::vector<double>& glW = gl->GetWeights();
    const std::vector<double>& glX = gl->GetAbscissas();
    double range = 0.0;
    vrange[0]    = range;
    for (int je=1; je<nekins; ++je) {
      const double ek0 = vekins[je-1];
      const double del = vekins[je] - ek0;
      double res = 0.0;
      for (int igl=0; igl<ngl; ++igl) {
        const double xi  = ek0 + del*glX[igl];
        const double val = sp.GetValueAt(xi);
        if (val>0.0) {
          res += glW[igl]/val;
        }
      }
      vrange[je] = vrange[je-1] + res*del;
    }
    // set a Spline util on this range data
    sp.SetUpSpline(vekins.data(), vrange.data(), nekins);
    // use this spline to fill in the restr. range data (that has a different ekin grid)
    for (int ie=0; ie<numEKin; ++ie) {
      double ekin = theEKinGrid[ie];
      elData.fDataPerMaterial[im]->fRange[ie] = sp.GetValueAt(ekin);
    }
    // set the Spline on the range
    elData.fDataPerMaterial[im]->fSpRange        = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fRange, numEKin);
    // set the Spline for the inverse range
    elData.fDataPerMaterial[im]->fSpInverseRange = new Spline(elData.fDataPerMaterial[im]->fRange, theEKinGrid,  numEKin);
    //
    // Compute the the MFP for bremsstrahlung
    // build the kinetic energy grid :
    // NOTE: emin = max (cut_e-, cut_gam.) since cut_e- is also a tracking cut for e-
    nekins = elData.fNumEKinData;
    ekin0  = std::max(electronCutInEnergy, gammaCutInEnergy);
    ekin1  = elData.fEKinGrid[elData.fNumEKinData-1];
    lekin0 = std::log(ekin0);
    ldekin = std::log(ekin1/ekin0)/(nekins-1);
    for (int je = 0; je<nekins; ++je) {
      double ekin  = (je==0) ? ekin0 : (je==nekins-1 ? ekin1 : std::exp(lekin0 + je*ldekin));
      // compute the restricted MFP for bremsstrahlung interaction
      double mxsec = std::max(1.0E-20, theBrem.CrossSectionPerVolume(mat, part, ekin, gammaCutInEnergy));
      elData.fDataPerMaterial[im]->fBremIMFPEkin[je] = ekin;
      elData.fDataPerMaterial[im]->fBremIMFP[je]     = mxsec;
    }
    elData.fDataPerMaterial[im]->fSpBremIMFP = new Spline(elData.fDataPerMaterial[im]->fBremIMFPEkin, elData.fDataPerMaterial[im]->fBremIMFP, nekins);

    //
    // compute the MFP for Moller(ionisation)
    // build the kinetic energy grid :
    // NOTE: emin = 2xcut_e- since both the 2 final e- needs to be above cut_e-
    nekins = elData.fNumEKinData;
    ekin0  = 2.0*electronCutInEnergy;
    ekin1  = elData.fEKinGrid[elData.fNumEKinData-1];
    lekin0 = std::log(ekin0);
    ldekin = std::log(ekin1/ekin0)/(nekins-1);
    for (int je = 0; je<nekins; ++je) {
      double ekin  = (je==0) ? ekin0 : (je==nekins-1 ? ekin1 : std::exp(lekin0 + je*ldekin));
//      std::cout << "ekin = " << ekin << std::endl;
      // compute the restricted MFP for Moller interaction
      double mxsec = std::max(1.0E-20, theIoni.CrossSectionPerVolume(mat, part, ekin, electronCutInEnergy));
      elData.fDataPerMaterial[im]->fMollerIMFPEkin[je] = ekin;
      elData.fDataPerMaterial[im]->fMollerIMFP[je]     = mxsec;
    }
    elData.fDataPerMaterial[im]->fSpMollerIMFP = new Spline(elData.fDataPerMaterial[im]->fMollerIMFPEkin, elData.fDataPerMaterial[im]->fMollerIMFP, nekins);

  }


  delete gl;
}


// computes the (total i.e. elastic with the inelastic correction) scattering
// strength that is allowed (according to the given stepping paraneters) in
// the given material at the given energy
// K_1(E) is computed according to the K_1(E_0) = int_0^s{ds' 1/lambda_1(s)} =
// = s\lambda_1(E_0) since energy dependence along the step is neglected
// Morever, different `s` (i.e. maximally allowed step lenght) values are used
// at low and high energies and connected smoothly:
// s(E) = s_low [s_high/s_low]^f(x) with x = log(E_0/ e_cross) and
// f(x) = arctang[sharp x]/pi + 1/2 such that f(x) is a smooth transition
// between 0 and 1 reaching 1/2 at x = 0. The sharpness of the transition is
// controlled by `sharp` (>0) parameter: the lower this value the faster the
// transition form 0  to 1 around x = 0. Since x = log(E/E_cross), the `E_cross`
// parameter (energy in MeV) will determine the enrgy around which and above
// s will be closer to `s_high` parameter than `s_low` since:
// log(E/E_cross) is  < 0, ~0 and > 0 when E < E_corss, E ~ E_cross and
// E > E_cross and the f[log(E/E_corss)] alues will be ~0, ~1/2 and ~1.
// Therefore, these 3 cases corresponds to s ~ s_low when E << E_cross,
// s ~ s_high when E >> E_corss and s ~ s_low [s_high/s_low]^1/2 when E ~ E_cross.
void   InitScatteringData(ElectronData& elData, double parShigh, double parSlow, double parEcross) {
  const double kInvPi =  1.0/3.1415926535897932;
  const double kSharp =  5.0;
  // get the list of materials from
  int numMat = elData.fNumMaterial;
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
//    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
//    std::cout << " material: " << mat->GetName() << std::endl;
    // get the scattering power correction
    double totScPowCor  = elData.GetTotalScatteringCorrection(im);
    double* theEKinGrid = elData.fEKinGrid;
    int     numEKin     = elData.fNumEKinData;
    for (int ie=0; ie<numEKin; ++ie) {
      double ekin  = theEKinGrid[ie];
      // compute the maximally allowed scatering power
      const double fx  = std::atan(kSharp*std::log(ekin/parEcross))*kInvPi +0.5;
//      double allowedK1 = elData.GetInvTotal1rstTransportMFP(ekin, im)*parSlow*std::pow(parShigh/parSlow,fx);
      double allowedK1 = totScPowCor*parSlow*std::pow(parShigh/parSlow,fx)/elData.fDataPerMaterial[im]->fFirstTransportMFP[ie];
      elData.fDataPerMaterial[im]->fMaxAllowedScatPower[ie] = allowedK1;
    }
    // set a Spline on this
    elData.fDataPerMaterial[im]->fSpMaxAllowedScatPower = new Spline(theEKinGrid, elData.fDataPerMaterial[im]->fMaxAllowedScatPower, numEKin);
  }
}

/*
void   InitGSData(ElectronData& elData, double electronCutInEnergy) {
  double eMin = electronCutInEnergy;
  double eMax = elData.fEKinGrid[elData.fNumEKinData-1];
  // number of discrete energy grid points and (common) cumulative values
  const int kNumEkin     = 256;
  const int kNumCumVals  =  64;
  // compute the delta of the (common) cumulative grid
  const double kCumDelta = 1./(kNumCumVals-1);

  int numMat = elData.fNumMaterial;
  // a (dense) initial `u` transformed variable grid on [0,1] and the q2+ transformed
  // GS pdf over it: this will be used to generate the final representation
  const int kNumUvals = 128;
  // utility containers
  std::vector<double>  theUgrid(kNumUvals, 0.0);
  std::vector<double>  theGSqPDF(kNumUvals, 0.0);
  std::vector<double>  theCumDistr(kNumUvals, 0.0);
  std::vector<double>  theParA(kNumUvals, 0.0);
  std::vector<double>  theParB(kNumUvals, 0.0);
  std::vector<double>  aDDum(kNumUvals, 0.0);
  std::vector<int>     aIDum(kNumUvals, 0.0);
  // fill the initial, dense u-grid
  double du = 1./(kNumUvals-1.);
  int iu = 0;
  for (; iu<kNumUvals-1; ++iu) {
    theUgrid[iu] = iu*du;
  }
  theUgrid[iu] = 1.0;
  //
  // these will store the final representation of a given GS q distribution (i.e. as
  // `a` and `b` parameters of a numerical inversion of the corresponding cumulative
  // at discrete `u` values that corresponds to a common, fixed cumulatove grid)
  std::vector<double>  theFinalUgrid(kNumCumVals, 0.0);
  std::vector<double>  theFinalParA(kNumCumVals, 0.0);
  std::vector<double>  theFinalParB(kNumCumVals, 0.0);
  // the common cumulatove grid
  std::vector<double>  theCumValueGrid(kNumCumVals, 0.0);
  int icum = 0;
  for ( ; icum < kNumCumVals-1; ++icum) {
    theCumValueGrid[icum] = icum*kCumDelta;
  }
  theCumValueGrid[icum] = 1.0;
  //
  // an AliasTable for peparing rational approximation based inverson
  AliasTable theAlias;
  Spline     sp;
  //
  // compute the transformed q2+ GS pdf-s at each E_i kinetic energy point for
  // each material
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
    std::cout << " GS data for material: " << mat->GetName() << std::endl;
    GSDtrData* theGSData = new GSDtrData(eMin, eMax, kNumEkin, kNumCumVals);
    // go over the kinetic energy grid
    for (int ie=0; ie<kNumEkin; ++ie) {
      // for each discrete E_i value, compute
      double ekin   = theGSData->fEKinGrid[ie];
//      std::cout << " imat = " << im << " ie = " << ie << " E_i = " << ekin << std::endl;
      // 1. get the maximally allowed scattering strength
      // NOTE: InitScatteringData must have been invoked before
      double theK1  = elData.GetMaxAllowedScatPow(ekin, im);
      // 2. compute G1 = lambda_el/lambda_1 ( see Eq.(7))
      double theG1  = elData.GetInvTotal1rstTransportMFP(ekin, im)/elData.GetInvTotalElasticMFP(ekin, im);
      // 3. since K1 is approximated as s/lambda_1 (such that s is a function of the enrgy
      //    and the 3 stepping parameter), one can get the corresponding s step length as
      //    s = K1*lambda_1. Moreover, K1/G1 = (s/lambda_1)/(lambda_e/lambda_1)= s/lambda_e
      double theNel = theK1/theG1;
      // 4. to find the screening parameter in the screened Rutherford DCS that reproduce this G1(x)
      double theParScreening = FindScreeingParameter(theG1);
      // 5. find the optimal transformation parameter that corresponds to this s/lambda_e and
      //    screening parameter
//      std::cout << "    theNel = " << theNel << " theParScreening = " << theParScreening  << std::endl;
      double theParTransform = GS::Instance().ComputeOptimalTransformationParameter(theNel, theParScreening);
      theGSData->fGSDtrData[ie].fTransformParam = theParTransform;
      // 6. compute the transformed q2+ GS function at this s/lambda_e, parScreening and optimal
      //    transformation parameter `a` over the `u` transformed variable grid
      iu = 0;
      for (; iu<kNumUvals; ++iu) {
        theGSqPDF[iu]= GS::Instance().ComputeTransformedGSDistribution(theNel, theParScreening, theParTransform, theUgrid[iu]);
      }
      // 7. prepare a rational approximation based inversion of the (dense) cumulative
      theAlias.PreparRatinTable(theUgrid.data(), theGSqPDF.data(), theCumDistr.data(), theParA.data(), theParB.data(), aDDum.data(),
                                aIDum.data(), kNumUvals);
      // 8. use this inversion to find the `u` values that correspond to the (common) discrete
      //    cumulative values
      theFinalUgrid[0]             = 0.0;
      theFinalUgrid[kNumCumVals-1] = 1.0;
      icum=1;
      for (; icum<kNumCumVals-1; ++icum) {
        double cumValue = theCumValueGrid[icum];
        // i such that cum_i <= cumValue < cum_{i+1}
        int indxl = std::lower_bound(theCumDistr.begin(), theCumDistr.end(), cumValue) - theCumDistr.begin() - 1;
        // compute the corresponding `u` value by using ratin based numerical inversion
        const double delta = theCumDistr[indxl + 1] - theCumDistr[indxl];
        const double aval  = cumValue-theCumDistr[indxl];
        const double dum1  = (1.0 + theParA[indxl] + theParB[indxl]) * delta * aval;
        const double dum2  = delta * delta + theParA[indxl] * delta * aval + theParB[indxl] * aval * aval;
        theFinalUgrid[icum]= theUgrid[indxl] + dum1 / dum2 * (theUgrid[indxl + 1] - theUgrid[indxl]);
      }
      // 9. recompute now the parameters of the rational inversion for this `u` grid
      // set-up a spline on the pdf in order to get values at new `u` points
      sp.SetUpSpline(theUgrid.data(), theGSqPDF.data(), kNumUvals);
      icum = 0;
      for ( ; icum<kNumCumVals-1; ++icum) {
        // the theCumValueGrid[icum+1] - theCumValueGrid[icum] is a const kCumDelta
        double dum1 = kCumDelta / (theFinalUgrid[icum + 1] - theFinalUgrid[icum]);
        double pj   = sp.GetValueAt(theFinalUgrid[icum]);
        double pjp1 = sp.GetValueAt(theFinalUgrid[icum+1]);
        double parb = 1.0 - dum1 * dum1 / (pj*pjp1);
        double para = dum1 / pj - 1.0 - parb;
        theFinalParA[icum] = para;
        theFinalParB[icum] = parb;
      }
      // 10. save final data for this E_i into the GSDtrData
      theGSData->fGSDtrData[ie].fCumData.resize(kNumCumVals);
      icum = 0;
      for ( ; icum<kNumCumVals; ++icum) {
        theGSData->fGSDtrData[ie].fCumData[icum].fVarU  = theFinalUgrid[icum];
        theGSData->fGSDtrData[ie].fCumData[icum].fParmA = theFinalParA[icum];
        theGSData->fGSDtrData[ie].fCumData[icum].fParmB = theFinalParB[icum];
      }
    }
    // store the pointer to this GSDtrData in the ElectronData for this material
    elData.fDataPerMaterial[im]->fTheGSData = theGSData;
  }
}
*/

void   InitGSData(ElectronData& elData, double electronCutInEnergy) {
  double eMin = electronCutInEnergy;
  double eMax = elData.fEKinGrid[elData.fNumEKinData-1];
  // number of discrete energy grid points and (common) cumulative values
  const int kNumEkin     = 256;
  const int kNumCumVals  =  64;
  // compute the delta of the (common) cumulative grid
  const double kCumDelta = 1./(kNumCumVals-1);

  int numMat = elData.fNumMaterial;
  // a (dense) initial `u` transformed variable grid on [0,1] and the q2+ transformed
  // GS pdf over it: this will be used to generate the final representation
  const int kNumUvals = 1024;
  // utility containers
  std::vector<double>  theUgrid(kNumUvals, 0.0);
  std::vector<double>  theGSqPDF(kNumUvals, 0.0);
  std::vector<double>  theCumDistr(kNumUvals, 0.0);
  std::vector<double>  theParA(kNumUvals, 0.0);
  std::vector<double>  theParB(kNumUvals, 0.0);
  std::vector<double>  aDDum(kNumUvals, 0.0);
  std::vector<int>     aIDum(kNumUvals, 0.0);
  // fill the initial, dense u-grid
  double du = 1./(kNumUvals-1.);
  int iu = 0;
  for (; iu<kNumUvals-1; ++iu) {
    theUgrid[iu] = iu*du;
  }
  theUgrid[iu] = 1.0;
  //
  // these will store the final representation of a given GS q distribution (i.e. as
  // `a` and `b` parameters of a numerical inversion of the corresponding cumulative
  // at discrete `u` values that corresponds to a common, fixed cumulatove grid)
  std::vector<double>  theFinalUgrid(kNumCumVals, 0.0);
  std::vector<double>  theFinalParA(kNumCumVals, 0.0);
  std::vector<double>  theFinalParB(kNumCumVals, 0.0);
  // the common cumulatove grid
  std::vector<double>  theCumValueGrid(kNumCumVals, 0.0);
  int icum = 0;
  for ( ; icum < kNumCumVals-1; ++icum) {
    theCumValueGrid[icum] = icum*kCumDelta;
  }
  theCumValueGrid[icum] = 1.0;
  //
  // an AliasTable for peparing rational approximation based inverson
  AliasTable theAlias;
  Spline     sp;
  //
  // compute the transformed q2+ GS pdf-s at each E_i kinetic energy point for
  // each material
  for (G4int im=0; im<numMat; ++im) {
    // the first material is the 'Universe' so skipp that one
//    G4Material*  mat = G4NistManager::Instance()->GetMaterial(im+1);
//    std::cout << " GS data for material: " << mat->GetName() << std::endl;
    GSDtrData* theGSData = new GSDtrData(eMin, eMax, kNumEkin, kNumCumVals);
    // go over the kinetic energy grid
    for (int ie=0; ie<kNumEkin; ++ie) {
      // for each discrete E_i value, compute
      double ekin   = theGSData->fEKinGrid[ie];
//      std::cout << " imat = " << im << " ie = " << ie << " E_i = " << ekin << std::endl;
      // 1. get the maximally allowed scattering strength
      // NOTE: InitScatteringData must have been invoked before
      double theK1  = elData.GetMaxAllowedScatLength(ekin, im);
      // 2. compute G1 = lambda_el/lambda_1 ( see Eq.(7))
      double theG1  = elData.GetInvTotal1rstTransportMFP(ekin, im)/elData.GetInvTotalElasticMFP(ekin, im);
      // 3. since K1 is approximated as s/lambda_1 (such that s is a function of the enrgy
      //    and the 3 stepping parameter), one can get the corresponding s step length as
      //    s = K1*lambda_1. Moreover, K1/G1 = (s/lambda_1)/(lambda_e/lambda_1)= s/lambda_e
      double theNel = theK1/theG1;
      // 4. to find the screening parameter in the screened Rutherford DCS that reproduce this G1(x)
      double theParScreening = FindScreeingParameter(theG1);
      // 5. find the optimal transformation parameter that corresponds to this s/lambda_e and
      //    screening parameter
//      std::cout << "    theNel = " << theNel << " theParScreening = " << theParScreening  << std::endl;
      double theParTransform = GS::Instance().ComputeOptimalTransformationParameter(theNel, theParScreening);
      theGSData->fGSDtrData[ie].fTransformParam = theParTransform;
      // 6. compute the transformed q2+ GS function at this s/lambda_e, parScreening and optimal
      //    transformation parameter `a` over the `u` transformed variable grid
      iu = 0;
      for (; iu<kNumUvals; ++iu) {
        theGSqPDF[iu]= GS::Instance().ComputeTransformedGSDistribution(theNel, theParScreening, theParTransform, theUgrid[iu]);
      }
      // 7. prepare a rational approximation based inversion of the cumulative such that the
      //    cumulative values are (the common) equally spaced ones on the [0,1]
      theAlias.PreparRatinTable(theUgrid.data(), theGSqPDF.data(), kNumUvals,  theFinalUgrid.data(), theFinalParA.data(), theFinalParB.data(), kNumCumVals);
      // 10. save final data for this E_i into the GSDtrData
      theGSData->fGSDtrData[ie].fCumData.resize(kNumCumVals);
      icum = 0;
      for ( ; icum<kNumCumVals; ++icum) {
        theGSData->fGSDtrData[ie].fCumData[icum].fVarU  = theFinalUgrid[icum];
        theGSData->fGSDtrData[ie].fCumData[icum].fParmA = theFinalParA[icum];
        theGSData->fGSDtrData[ie].fCumData[icum].fParmB = theFinalParB[icum];
      }
    }
    // store the pointer to this GSDtrData in the ElectronData for this material
    elData.fDataPerMaterial[im]->fTheGSData = theGSData;
  }
}
