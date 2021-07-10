#include "SimDPMLike.hh"


#include "Geom.hh"


#include "SimMaterialData.hh"


#include "SimElectronData.hh"

#include "SimITr1MFPElastic.hh"
#include "SimMaxScatStrength.hh"

#include "SimIMFPMoller.hh"
#include "SimIMFPBrem.hh"
#include "SimStoppingPower.hh"

#include "SimMollerTables.hh"
#include "SimSBTables.hh"
#include "SimGSTables.hh"



#include "SimPhotonData.hh"

#include "SimIMFPMaxPhoton.hh"
#include "SimIMFPPhoton.hh"

#include "SimKNTables.hh"


#include "Random.hh"
#include "Track.hh"
#include "TrackStack.hh"

#include <cstdio>


void   Simulate(int nprimary, double eprimary, bool iselectron, double lbox, SimMaterialData& matData, SimElectronData& elData, SimPhotonData& phData, int geomIndex) {
  const double kPI            = 3.1415926535897932;
  const double kEMC2          = 0.510991;
  //
  // create the simple geometry
  Geom geom(lbox, &matData, geomIndex);
  //
  const double theElectronCut = matData.fElectronCut;
  const double theGammaCut    = matData.fGammaCut;

  //
  // `theRZ0` such that the particle hits the first voxel/box with `iz` index of 0
  // from the left (i.e. --> |0|1|...)
  // NOTE: it will be located inside the box with z-index of 0 when `DistanceToBoundary`
  //       will be called since theRZ0 is within (half) tolerance outside.
  const double theRZ0 = -0.5*lbox;
  //
  // simulate `nprimary` histories: prite some infomation reagrding the progress
  const int kFraction = nprimary/10;
  std::cout << "\n === Start simulation of N = " << nprimary << " events === \n" << std::endl;
  for (int iPrimary=0; iPrimary<nprimary; ++iPrimary) {
    if (iPrimary%kFraction==0) {
      std::cout << " == Starting Primary = " << iPrimary << std::endl;
    }
    //
    // insert a track into the global track-stack with the properties of the primary
    // (note, that the Track is Reset() inside before we get back a refernce to that)
    Track& aPrimaryTrack = TrackStack::Instance().Insert();
    // set the primary kinetic energy (position, direcion, type etc is fine)
    aPrimaryTrack.fEkin         = eprimary;             // Primary kinetic energy in [MeV]
    aPrimaryTrack.fType         = iselectron ? -1 : 0 ; // e-(-1) or photon(0)
    //
    aPrimaryTrack.fDirection[0] = 0.0;                  // initial direction is [0,0,1]
    aPrimaryTrack.fDirection[1] = 0.0;
    aPrimaryTrack.fDirection[2] = 1.0;
    //
    aPrimaryTrack.fPosition[0] = 0.0;                   // initial position is [0,0, theRZ0]
    aPrimaryTrack.fPosition[1] = 0.0;
    aPrimaryTrack.fPosition[2] = theRZ0;
    //
    // While the stack becomes empty: pop the next track (secondary except the very firts case)
    // and keep tracking till its kinetic energy becomes zero.
    // NOTE: the stack might be populated with secondary tracks while keepd tracking a given track
    Track track; // this is the current Track that is under tracking
    while (-1 < TrackStack::Instance().PopIntoThisTrack(track)) {
          // compute the distance to the boundary
          // (This also sets the box indices so the material index can be obtained)
          // init the step length to this distance to boundary
//          double step =
          geom.DistanceToBoundary(track.fPosition, track.fDirection, track.fBoxIndx);
          // get the current material index: i.e. the base material of the voxel
          // NOTE: vacuum is indicated by the special voxel material index of -1
          //       Stop tracking if the track is in a voxel filled with vacuum.
          int theVoxelMatIndx = geom.GetMaterialIndex(track.fBoxIndx);
          if (theVoxelMatIndx<0) {
            continue;
          }
          // set the material index otherwise
          track.fMatIndx      = theVoxelMatIndx;
          // Use the dedicated tracking for photons if we have one in hand now:
          if (track.fType == 0) {
            KeepTrackingPhoton(phData, matData, geom, track);
            continue;
          }
          //
          // Track e-/e+ ortherwise:
          //
          // WE ASSUME HERE NOW THAT EACH VOXEL IS A CLEAR MATERIAL SO WE WILL
          // USE theVoxelMaterialDensity = theVoxelBaseMaterialDensity. HOWEVER,
          // THIS MIGH BE CHANGED ANYTIME WHEN THE GEOMETRY CAN PROVIDE A SUITABLE
          // VOXEL MATERIAL DENSITY.
          //
          // NOTE: density must be in g/cm3 units. (cannot be vacuum at this point)
          double theVoxelMatDensity = geom.GetVoxelMaterialDensity(theVoxelMatIndx);
          //
          // this will be used to alter between an hinge and remaining sub MSC-step
          bool   isMSCHinge    = true;
          // this will be used to keep track of the pre-step point energy that we need
          // when we reach the msc hinge point (init. is not important cause will be set below)
          double theEkin0      = track.fEkin;
          //
          //
          // Compute the initial number of mfp left till the different interactions
          // which is -ln(R) = s/mfp or U(R) = s/tr1-mfp for the MSC hinge interaction:
          //
          // 1. elastic interaction with msc: s/tr1mfp and sample the hinge position as well
          //    NOTE: data are used for the reference material and not scalled by the density
          //          so numTr1MFP ~ s_max(E)/tr1-mfp(E) [no units]
          double numTr1MFP    = elData.GetMaxScatStrength()->GetMaxScatStrength(track.fEkin);
          // travell this #tr1-mfp after the MSC-hinge took place
          double numTr1MFP0   = Random::UniformRand()*numTr1MFP;
          // travel this #tr1-mfp till the MSC-hinge
          numTr1MFP          -= numTr1MFP0;
          //
          // 2. Moller:
          // NOTE: as in DPM, the mfp for Moller will be assumed to have no energy dependence!!!
          //       So the update of the number-of-mfp-left will contain only the material
          //       dependence related scaling (that's again approximate but DPM does this).
          //       This material dependent update(scaling) relies on the \lam ~ [A/(Z\rho)]
          //       dependence: \lam' = \lam_ref(E') [Z\rho/A]_ref [A/(Z\rho)]_actual (such
          //       that in case of Moller (in DPM) even E'=E_0 in each step as dicussed above).
          double numMollerMFP = -std::log(Random::UniformRand());
          // again, the reference material Moller IMFP value is used
          double invMollerMFP = elData.GetIMFPMoller()->GetIMFPPerDensity(track.fEkin);
          //
          // 3. bremsstrahlung:
          // NOTE: IMFP for bremsstrahlung is important so it's interpolated by using values for
          //       the actual material and kinetic energy in the `KeepTracking` code. Here we
          //       sample a `number-of-interaction-left` value, i.e. -ln(R) only (=s/imfp(E,mat))
          double numBremMFP = -std::log(Random::UniformRand());
          //
          // Start and keep tracking this e- or e+ track while its stopped, i.e. its energy becomes
          // zero or goes away (in vacuum)
          while (track.fEkin > 0.0) {
            //
            // Now we can keep tracking this e- (even trhough boxes with different materials),
            // by decreasing these 3 above initial number-of-mfp/tr1-mfp-left at each step by the
            // correspnding n' = n - dS/mfp' or - dS/tr1-mfp as long as:
            // a.  any of these above 3 number of mfp/tr1-mfp goes to zero ==> the correspnding
            //     discrete Brem.(1), Moller(2) interaction happens or (3) either the MSC
            //     hinge or MSC step end point is reached
            // b.  (4) the e- energy drops below the e- tracking cut which is equal to the
            //     secondary e- production threshold (end of this history)
            // c.  the e- leaves the geometry, i.e. goes into vacuum (its energy set to 0 so return with 4)
            int whatHappend = KeepTrackingElectron(elData, matData, geom, track, numTr1MFP, numMollerMFP, invMollerMFP, numBremMFP);
            //
            theVoxelMatIndx = track.fMatIndx;
            // terminate if the tarck is in vacuum (its energy has been set to zero)
            if (theVoxelMatIndx<0) {
              continue;
            }
            theVoxelMatDensity = geom.GetVoxelMaterialDensity(theVoxelMatIndx);
            switch (whatHappend) {
              // (1) discrete bremsstrahlung interaction should be sampled:
              //     - sample energy transfer to the photon (if any)
              case 1 : {
                         // perform bremsstrahlung
                         if (track.fEkin > theGammaCut) {
                           PerformBrem(track, elData.GetTheSBTables());
                         }
                         // check if the post-interaction electron energy dropped
                         // below the tracking cut and stop tracking if yes
                         if (track.fEkin < theElectronCut) {
                           // deposit the whole energy and stop tracking
                           geom.Score(track.fEkin, track.fBoxIndx[2]);
                           track.fEkin = 0.0;
                           // perform annihilation in case of e+
                           if (track.fType == +1) {
                             PerformAnnihilation(track);
                           }
                           break;
                         }
                         // if the primary electron (or e+) survived, i.e. if we are here
                         // then, re-sample the #mfp to travel till the next brem event
                         // and break
                         numBremMFP = -std::log(Random::UniformRand());
                         break;
                       }
              // (2) discrete Moller interaction is sampled:
              // NOTE: no energy dependence is considered in case of Moller in DPM so
              //       the numMollerMFP is evaluated only at this point and assumed to be
              //       constant along the entire `KeepTracking` part (only material scaling is applied).
              //       Furthermore, the kinetic energies of both the post interaction priary
              //       and seconday electrons are guarantied to be above the secondary electron
              //       production threshold so no need to check if their energy droppe below after
              //       the interaction
              //       Furthermore, note that Moller interaction is independent from Z
              case 2 : {
                          // energy transfer is sampled in E0 units and E0 must be at least 2cut
                          const double kcut = theElectronCut/track.fEkin;
                          double secEkin = 0.0;
                          double secCost = 1.0;
                          if (kcut<0.5) {
                            const double k2EMC2 = 2.0*kEMC2;
                            secEkin = elData.GetTheMollerTables()->SampleEnergyTransfer(track.fEkin, Random::UniformRand(), Random::UniformRand(), Random::UniformRand());
                            secCost = std::min(1.0, std::sqrt(secEkin*(track.fEkin+k2EMC2)/(track.fEkin*(secEkin+k2EMC2))));
                          }
                          // if interaction was possible. i.e. both seconday and post-intecation primary ekin > ecut
                          if (secEkin>0.0) {
                            // insert the secondary e- track into the stack
                            Track& aTrack        = TrackStack::Instance().Insert();
                            aTrack.fType         = -1;
                            aTrack.fEkin         = secEkin;
                            aTrack.fMatIndx      = track.fMatIndx;
                            aTrack.fPosition[0]  = track.fPosition[0];
                            aTrack.fPosition[1]  = track.fPosition[1];
                            aTrack.fPosition[2]  = track.fPosition[2];
                            aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
                            aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
                            aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
                            const double sint    = std::sqrt((1.0+secCost)*(1.0-secCost));
                            const double phi     = 2.0*kPI*Random::UniformRand();
                            aTrack.fDirection[0] = sint*std::cos(phi);
                            aTrack.fDirection[1] = sint*std::sin(phi);
                            aTrack.fDirection[2] = secCost;
                            RotateToLabFrame(aTrack.fDirection, track.fDirection);
                            // decrease primary energy: DMP do not deflect the primary
                            track.fEkin -= secEkin;
                         }
                         // Resample #mfp left and interpolate the IMFP since the enrgy has been changed.
                         // Again, the reference material Moller IMFP value is used
                         numMollerMFP = -std::log(Random::UniformRand());
                         invMollerMFP = elData.GetIMFPMoller()->GetIMFPPerDensity(track.fEkin);
                         break;
                       }
              // (3) msc interaction happend: either hinge or just or end of an MSC step
              case 3 : {
                         if (isMSCHinge) {
                           // Sample angular deflection from GS distr. and apply it
                           // -----------------------------------------------------
                           const double dum0 = elData.GetTheGSTables()->SampleAngularDeflection(theEkin0, Random::UniformRand(), Random::UniformRand());
                           const double cost = std::max(-1.0, std::min(1.0, dum0));
                           const double sint = std::sqrt((1.0-cost)*(1.0+cost));
                           // smaple \phi: uniform in [0,2Pi] <== spherical symmetry of the scattering potential
                           const double phi  = 2.0*kPI*Random::UniformRand();
                           // compute new direction (relative to 0,0,1 i.e. in the scattering frame)
                           double u1 = sint*std::cos(phi);
                           double u2 = sint*std::sin(phi);
                           double u3 = cost;
                           // rotate new direction from the scattering to the lab frame
                           RotateToLabFrame(u1, u2, u3, track.fDirection[0], track.fDirection[1], track.fDirection[2]);
                           // update track direction
                           track.fDirection[0] = u1;
                           track.fDirection[1] = u2;
                           track.fDirection[2] = u3;
                           // -----------------------------------------------------
                           // set the #tr1-mfp left to the remaining, i.e. after hinge part
                           numTr1MFP   = numTr1MFP0;
                           // the end point is the next msc stop and not the hinge
                           isMSCHinge = false;
                         }  else {
                           // end point so resample #tr1-mfp left and the hinge point
                           theEkin0   = track.fEkin;
                           // again, the reference material K_1(E) is used
                           numTr1MFP  = elData.GetMaxScatStrength()->GetMaxScatStrength(theEkin0);
                           // travell this #tr1-mfp after the MSC-hinge took place
                           numTr1MFP0 = Random::UniformRand()*numTr1MFP;
                           // travell this #tr1-mfp till the MSC-hinge
                           numTr1MFP -= numTr1MFP0;
                           // hinge will be the next msc stop
                           isMSCHinge = true;
                         }
                         break;
                      }
              // (4) the kinetic energy dropped below the tracking cut so the particle is stopped
              case 4 :  // nothng to do now: track.fEkin should be zero et this point so the
                        // "traking" while loop should be terminated after this `break`
                        break;
            }
          } // end of tracking while loop
    } // end of this history (i.e. end of tracking even the last primary track from this history)
  } // end of tracking all the primaries (i.e. end of tracking the last primary with all its secondaries)
  std::cout << "\n === End simulation of N = " << nprimary << " events === \n" << std::endl;
  // write histograms
  geom.Write("hist.sim", nprimary);
}


int KeepTrackingElectron(SimElectronData& elData, SimMaterialData& matData, Geom& geom, Track& track, double& numTr1MFP, double& numMollerMFP, double invMollerMFP, double& numBremMFP) {
  int whatHappend = 0;
  // compute the distance to boundary: this will be the current value of the maximal step length
  double stepGeom = geom.DistanceToBoundary(track.fPosition, track.fDirection, track.fBoxIndx);
  // When the particle goes further than the pre-defined (+-10 [cm] along the xy plane, +10 cm
  // in the +z and the voxel-size in the -z direction) then `DistanceToBoundary` returns -1.0
  // indicating that the particle left the geometry: we stop tracking this particle
  if (stepGeom < 0.0) {
    // left the geometry
    track.fEkin = 0.0;
    return 4;
  }
  //
  const double theElectronCut = matData.fElectronCut;
  //
  // The Moller mfp energy dependence is not considered in DPM so we do the same:
  // we will keep using the Moller mfp evaluated at the initial energy for the reference
  // material. Only the material scaling will be applied below:
  // \lam' = \lam_ref(E') [Z\rho/A]_ref [A/(Z\rho)]_actual (for 1/\lam' to compute delta #mfp)
  while (whatHappend==0) {
    // init the current step lenght to that of the distance to boundary: we might
    // or might not go that far (depending on what the 3 physics tells) but for
    // sure that is the maximum step length becasue 'boundary crossing interaction'
    // happens after travelling that far. So set what happen to that (i.e. = 0.)
    double stepLength = stepGeom;
    whatHappend = 0;
    // get the current material index (when computing distance to boundary track.fBoxIndx
    // is updated if needed, i.e. if we crossed a boundary)
    int theVoxelMatIndx = geom.GetMaterialIndex(track.fBoxIndx);
    // stop tracking if the track entered to vacuum i.e. voxel with material index of -1
    if (theVoxelMatIndx < 0) {
      track.fEkin = 0.0;
      return 4;
    }
    // set the material index otherwise
    track.fMatIndx      = theVoxelMatIndx;
    // the material scaling factor for the Moller inverse-mf: [A/(Z\rho/)]_ref [(Z\rho)/A]_actual
    // or more exactly its [A/Z)]_ref [(Z)/A]_actual part
    double scalMolMFP = matData.fMollerIMFPScaling[theVoxelMatIndx];
    // WE ASSUME HERE NOW THAT EACH VOXEL IS A CLEAR MATERIAL SO WE WILL
    // USE theVoxelMaterialDensity = theVoxelBaseMaterialDensity. HOWEVER,
    // THIS MIGH BE CHANGED ANYTIME WHEN THE GEOMETRY CAN PROVIDE A SUITABLE
    // VOXEL MATERIAL DENSITY.
    //
    // NOTE: density must be in g/cm3 units !!!!
    double theVoxelMatDensity = geom.GetVoxelMaterialDensity(theVoxelMatIndx);
    //
    // Here we compute the decrese of the #mfp/#tr1-mfp for the 3 interaction with the current,
    // maximum step length (i.e. the distance to boundary): #mfp' = #mfp - ds/mfp' or - ds/tr1_mfp
    // for MSC (where ' indicates values at the end point)
    //
    // compute the mid-point energy along this step by assuming:
    // - constant dEdx along the step, i.e. dEdx=dEdx(E_0) and dE = s dEdx --> E_mid = E_0 - 0.5 s dEdx
    // - the step equal to the current one, i.e. `stepLength` (dist. to boundary)
    // the restricted stopping power for this material: for the referecne material and scalled with the current density
    double theDEDX    = elData.GetDEDX()->GetDEDXPerDensity(track.fEkin, theVoxelMatIndx)*theVoxelMatDensity;
    // make sure that do not go below the minim e- energy
    double midStepE   = std::max(track.fEkin-0.5*stepLength*theDEDX, theElectronCut );
    // elastic: #tr1-mfp' = #tr1-mfp - ds/tr1-mfp' so the change in #tr1-mfp is ds/tr1-mfp' and
    //          1/mfp' is computed here
    double delNumTr1MFP    = elData.GetITr1MFPElastic()->GetITr1MFPPerDensity(midStepE, theVoxelMatIndx)*theVoxelMatDensity;
    // moller: see above the details
    double delNumMollerMFP = invMollerMFP*scalMolMFP*theVoxelMatDensity;
    // brem: #mfp' = #mfp - ds/mfp' with mfp = brem_mfp so the change in #mfp is ds/mfp' and
    //       1/mfp' is computed here
    double delNumBremMFP   = elData.GetIMFPBrem()->GetIMFPPerDensity(midStepE, theVoxelMatIndx)*theVoxelMatDensity;
    //
    //
    // Now we will see how far actually we go by trying to decrese each of the 3 #mfp/#tr1-mfp
    // by the changes in the number of mfp/tr1-mfp computed above as `delNum` :
    // - if we could manage to decrese all the 3 #mfp/tr1-mfp such that they are > 0
    //   then actually we reached the boundary: we cross and perform an other step.
    // - if any of the 3 #mfp/#tr1-mfp goes down to zero, then the lowest (i.e. shortest
    //   path to) will be considered to happen: the given inetraction need to beinvoked
    // In all cases, the energy loss along the given step needs to be computed!
    //
    // In DMP, the #mfp are tried to decresed with the current step length (that
    // is always the current shortest) as #mfp: n' = n - ds/mfp' then corrected
    // back if goes below zero, etc...
    // We compute the step length ds = n x mfp' (or n x tr1-mfp') i.e. that
    // would bring the given number #mfp/#tr1-mfp down to 0 (since n'=n-ds/mfp' or ds/tr1-mfp).
    // If this step is the shortest then we take this as current step length and we determine
    // the overal shortest and the corresponding interaction will happen (including boundary crossing).
    //
    double stepBrem    = numBremMFP/delNumBremMFP;
    if (stepBrem < stepLength) {
      // discrete bremsstrahlung (might) happen before reaching the boundary:
      stepLength  = stepBrem;
      whatHappend = 1;
    }
    double stepMoller  = numMollerMFP/delNumMollerMFP;
    if (stepMoller < stepLength) {
      // discrete Moller (might) happen even before bremsstrahlung:
      stepLength  = stepMoller;
      whatHappend = 2;
    }
    double stepElastic = numTr1MFP/delNumTr1MFP;
    if (stepElastic < stepLength) {
      // elastic interaction happens:
      // - either the hinge: sample and apply deflection and update numElMFP to the
      //                     remaining part
      // - end of 2nd step : nothing to do just resample the #mfp left since all
      //                     has been eaten up by the step lenght travelled so far
      // Before anything, refine the comutation of the mfp (i.e. the first transprt
      // mean free path in acse of elastic) regarding its energy dependence.
      // NOTE: that the 1/mfp values were computed at the mid-point energy (brem
      //       and elatsic since Moller is assumed to be constant), assuming that
      //       the geometry step will be taken and the dEdx is constant along this
      //       step (i.e. no energy dependence).
      //       Here we know that actually not the geometry step, but the stepElastic
      //       is taken since that is the shortest. So we recompute the mid-step-point
      //       energy according to the step lenght of stepElastic and re-evaluate
      //       the 1./mfp i.e. 1/tr1mfp at this energy value
      stepElastic = numTr1MFP/(elData.GetITr1MFPElastic()->GetITr1MFPPerDensity(track.fEkin, theVoxelMatIndx)*theVoxelMatDensity);
      midStepE    = std::max( track.fEkin-0.5*stepElastic*theDEDX, theElectronCut );
      delNumTr1MFP = elData.GetITr1MFPElastic()->GetITr1MFPPerDensity(midStepE, theVoxelMatIndx)*theVoxelMatDensity;
      // don't let longer than the original in order to make sure that it is still the
      // minimum of all step lenghts
      stepElastic = std::min(stepLength, numTr1MFP/delNumTr1MFP);
      stepLength  = stepElastic;
      whatHappend = 3;
    }
    //
    // At this point, we know the step lenght so we can decrease all #mfp by
    // substracting the delta #mfp = ds/#mfp that correspond to this final stepLength
    numBremMFP   = (whatHappend == 1) ? 0 : numBremMFP   - stepLength*delNumBremMFP;
    numMollerMFP = (whatHappend == 2) ? 0 : numMollerMFP - stepLength*delNumMollerMFP;
    numTr1MFP    = (whatHappend == 3) ? 0 : numTr1MFP    - stepLength*delNumTr1MFP;


    //
    // Compte the (sub-treshold, i.e. along step) energy loss:
    // - first the mid-step energy using the final value of the step lenght and the
    //   pre-step point dEdx (assumed to be constant along the step).
    midStepE     = std::max( track.fEkin-0.5*stepLength*theDEDX, theElectronCut );
    // - then the dEdx at this energy
    theDEDX      = elData.GetDEDX()->GetDEDXPerDensity(midStepE, theVoxelMatIndx)*theVoxelMatDensity;
    // - then the energy loss along the step using the mid-step dEdx (as constant)
    //   and the final energy
    double deltE = stepLength*theDEDX;
    double eFinal= track.fEkin-deltE;
    // check if energy dropped below tracking cut, i.e. below seconday e- production threshold
    // NOTE: HERE THERE IS A SUB-THRESHOLD TRACKING CONDITION IN DPM BUT WE WILL NEED TO SEE THAT !!!
    // ALSO: IF THE SELECTED STEP LENGHT BRINGS THE EKIN BELOW THRESHOLD WE DON'T TRY TO FIND THE
    //       STEP LENGTH (a smaller than selected) THAT ACTUALLY BRINGS EKIN EXACTLY TO THE THRESHOLD.
    //       SO THE TRACK LENGTH IS NOT PRECISE, SINCE WE KNOW THAT THE e- WAS STOPPED IN THIS VOLUME/BOX
    //       BUT WE DON'T COMPUTE THE EXCT POSITION
    if (eFinal < theElectronCut) {
      // deposit the whole energy and stop tracking
      track.fEdep  = track.fEkin;
      track.fEkin  = 0.0;
      // perform e+e- annihilation in case of e+
      if (track.fType == +1) {
        // annihilate the e+ at the correct position !!!
        track.fPosition[0] += track.fDirection[0]*stepLength;
        track.fPosition[1] += track.fDirection[1]*stepLength;
        track.fPosition[2] += track.fDirection[2]*stepLength;
        // make sure that the mat index is up to date
        geom.DistanceToBoundary(track.fPosition, track.fDirection, track.fBoxIndx);
        theVoxelMatIndx = geom.GetMaterialIndex(track.fBoxIndx);
        // check if we are in vacuum: nothing happens in that case
        if (theVoxelMatIndx < 0) {
          // kinetic energy is already zero
          return 4;
        }
        track.fMatIndx  = theVoxelMatIndx;
        PerformAnnihilation(track);
      }
      whatHappend  = 4;
    } else {
      track.fEkin  = eFinal;
      track.fEdep  = deltE;
    }
    //
    // Update particle position, track length etc.
    track.fPosition[0] += track.fDirection[0]*stepLength;
    track.fPosition[1] += track.fDirection[1]*stepLength;
    track.fPosition[2] += track.fDirection[2]*stepLength;
    // update cummulative track length
    track.fStepLenght   = stepLength;
    track.fTrackLength += stepLength;
    //
    // Score the continuous energy loss before going back to perform the discrete
    // interaction (whatHappend={1,2,3}) OR to terminate the tracking (whatHappend=4)
    // NOTE: we need to score before calling DistanceToBoundary again because that
    //       might positon the particle to the next volume.
    geom.Score(track.fEdep, track.fBoxIndx[2]);
    //
    // Compute distance to boundary if geometry limited this step:
    // - if geometry limited the step, the current position above is on a
    //   volume/box boundary
    // - when calling DistanceToBoundary such that the position is within half
    //   tolerance to a boudnary, the track position is updated to be on the
    //   other side, the distance to boundary from the new positon is computed
    //   and the x,y and z box coordinate indices are updated (so the material
    //   index will be the new one)
    if (whatHappend==0) {
      stepGeom = geom.DistanceToBoundary(track.fPosition, track.fDirection, track.fBoxIndx);
      if (stepGeom<0.0) {
        // left the geometry
        track.fEkin = 0.0;
        return 4;
      }
    }
  }
  return whatHappend;
}


void KeepTrackingPhoton(SimPhotonData& phData, SimMaterialData& matData, Geom& geom, Track& track) {
  const double kPI      = 3.1415926535897932;
  const double kEMC2    = 0.510991;
  const double kInvEMC2 = 1.0/kEMC2;
  //
  const double theElectronCut = matData.fElectronCut;
  const double theGammaCut    = matData.fGammaCut;
  //
  while (track.fEkin > 0.0) {
    // get the global max-macroscopic cross section and use it for samppling the
    // the length till the next photon intercation (that includes delta interaction
    // as well)
    double globalMaxMFP   = 1.0/phData.GetIMFPTotalMax()->GetIMFP(track.fEkin);
    double     stepLength = -globalMaxMFP*std::log(Random::UniformRand());
    // Update particle position, track length etc.
    track.fPosition[0] += track.fDirection[0]*stepLength;
    track.fPosition[1] += track.fDirection[1]*stepLength;
    track.fPosition[2] += track.fDirection[2]*stepLength;
    // update cummulative track length
    track.fStepLenght   = stepLength;
    track.fTrackLength += stepLength;
    // determine currecnt voxel index
    if (geom.DistanceToBoundary(track.fPosition, track.fDirection, track.fBoxIndx) < 0) {
      // left the geometry
      track.fEkin = 0.0;
      return;
    }
    //
    // check if any interaction happened
    int theVoxelMatIndx = geom.GetMaterialIndex(track.fBoxIndx);
    if (theVoxelMatIndx < 0) {
      // terminate because its in the vacuum
      track.fEkin = 0.0;
      return;
    }
    track.fMatIndx = theVoxelMatIndx;
    double theVoxelMatDensity = geom.GetVoxelMaterialDensity(theVoxelMatIndx);
    //
    double totalIMFP = phData.GetIMFPTotal()->GetIMFPPerDensity(track.fEkin, theVoxelMatIndx)*theVoxelMatDensity;
    //
    // P(no-inetaction) = 1.0-mxsecTotal/mxsecGlobalMax
    const double r1 = Random::UniformRand();
    double theProb = 1.0-totalIMFP*globalMaxMFP;
    if (r1 < theProb) {
      continue; // with the same globalMaxMFP since the enrgy did not changed !!!
    }
    //
    // otherwise: check which interaction happend P(i) = mxsec-i/mxsecTotal
    // compute cumulated probability of adding Compton prob
    theProb += phData.GetIMFPCompton()->GetIMFPPerDensity(track.fEkin, theVoxelMatIndx)*theVoxelMatDensity*globalMaxMFP;
    if (r1 < theProb) {
      // Compton interaction: Klein-Nishina like
      // the photon scattering angle and post-interafctin energy fraction
      const double theEps = phData.GetTheKNTables()->SampleEnergyTransfer(track.fEkin, Random::UniformRand(), Random::UniformRand(), Random::UniformRand());
      const double kappa  = track.fEkin*kInvEMC2;
      const double phCost = 1.0-(1.0-theEps)/(theEps*kappa); // 1- (1-cost)
      const double phEner = theEps*track.fEkin;
      const double phPhi  = 2.0*kPI*Random::UniformRand();
      const double elEner = track.fEkin-phEner;
      // insert the secondary e- only if ist energy is above the tracking cut
      // and deposit the corresponding enrgy otherwise
      if (elEner < theElectronCut) {
        geom.Score(elEner, track.fBoxIndx[2]);
        //geom.Score(elEner, track.fPosition[2]);
      } else {
        // insert secondary e- but first compute its cost
        const double e0 = track.fEkin*kInvEMC2;
        double elCost   = (1.0+e0)*std::sqrt((1.0-theEps)/(e0*(2.0+e0*(1.0-theEps))));
        //
        Track& aTrack        = TrackStack::Instance().Insert();
        aTrack.fType         = -1;
        aTrack.fEkin         = elEner;
        aTrack.fMatIndx      = track.fMatIndx;
        aTrack.fPosition[0]  = track.fPosition[0];
        aTrack.fPosition[1]  = track.fPosition[1];
        aTrack.fPosition[2]  = track.fPosition[2];
        aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
        aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
        aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
        const double sint    = std::sqrt((1.0+elCost)*(1.0-elCost));
        const double phi     = phPhi+kPI;
        aTrack.fDirection[0] = sint*std::cos(phi);
        aTrack.fDirection[1] = sint*std::sin(phi);
        aTrack.fDirection[2] = elCost;
        RotateToLabFrame(aTrack.fDirection, track.fDirection);
      }
      // update the photon properties:
      // stop the photon if its energy dropepd below the photon absorption threshold
      track.fEkin = phEner;
      if (track.fEkin < theGammaCut) {
        geom.Score(track.fEkin, track.fBoxIndx[2]);
        //geom.Score(track.fEkin, track.fPosition[2]);
        track.fEkin = 0.0;
        return;
      } else {
        double phSint = std::sqrt((1.0-phCost)*(1.0+phCost));
        double u1 = phSint*std::cos(phPhi);
        double u2 = phSint*std::sin(phPhi);
        double u3 = phCost;
        // rotate new direction from the scattering to the lab frame
        RotateToLabFrame(u1, u2, u3, track.fDirection[0], track.fDirection[1], track.fDirection[2]);
        // update track direction
        track.fDirection[0] = u1;
        track.fDirection[1] = u2;
        track.fDirection[2] = u3;
      }
      continue;
    }

    // compute cumulated probability of adding Pair-production prob
    theProb += phData.GetIMFPPairProd()->GetIMFPPerDensity(track.fEkin, theVoxelMatIndx)*theVoxelMatDensity*globalMaxMFP;
    if (r1 < theProb) {
      // Pair-production interaction:
      const double sumEkin = track.fEkin-2.0*kEMC2;
      // simple uniform share of the enrgy between the e- and e+ going to the
      // same direction as the original photon.
      // no difference between the e- and e+ transport till the end:
      // - when the e+ stops, 2 photons are emitted
      // we will assume that e1 is the e+
      double e1 = Random::UniformRand()*sumEkin;
      double e2 = sumEkin-e1;
      // insert the e- and e+ only if their energy is above the tracking cut
      // the e-
      if (e2 < theElectronCut) {
        geom.Score(e2, track.fBoxIndx[2]);
        //geom.Score(e2, track.fPosition[2]);
      } else {
        Track& aTrack        = TrackStack::Instance().Insert();
        aTrack.fType         = -1;
        aTrack.fEkin         = e2;
        aTrack.fMatIndx      = track.fMatIndx;
        aTrack.fPosition[0]  = track.fPosition[0];
        aTrack.fPosition[1]  = track.fPosition[1];
        aTrack.fPosition[2]  = track.fPosition[2];
        aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
        aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
        aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
        aTrack.fDirection[0] = track.fDirection[0];
        aTrack.fDirection[1] = track.fDirection[1];
        aTrack.fDirection[2] = track.fDirection[2];
      }
      // the e+
      if (e1 < theElectronCut) {
        geom.Score(e1, track.fBoxIndx[2]);
        //geom.Score(e1, track.fPosition[2]);
        PerformAnnihilation(track);
      } else {
        Track& aTrack        = TrackStack::Instance().Insert();
        aTrack.fType         = +1;
        aTrack.fEkin         = e1;
        aTrack.fMatIndx      = track.fMatIndx;
        aTrack.fPosition[0]  = track.fPosition[0];
        aTrack.fPosition[1]  = track.fPosition[1];
        aTrack.fPosition[2]  = track.fPosition[2];
        aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
        aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
        aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
        aTrack.fDirection[0] = track.fDirection[0];
        aTrack.fDirection[1] = track.fDirection[1];
        aTrack.fDirection[2] = track.fDirection[2];
      }
      // kill the primary photon
      track.fEkin = 0.0;
      return;
    }

    // if we are here then Photoelectric effect happens that absorbs the photon:
    // - score the current phton energy and stopp the photon
    geom.Score(track.fEkin, track.fBoxIndx[2]);
    //geom.Score(track.fEkin, track.fPosition[2]);
    track.fEkin = 0.0;
  };
}


void RotateToLabFrame(double &u, double &v, double &w, double u1, double u2, double u3) {
  double up = u1*u1 + u2*u2;
  if (up>0.) {
    up = std::sqrt(up);
    double px = u;
    double py = v;
    double pz = w;
    u = (u1*u3*px - u2*py)/up + u1*pz;
    v = (u2*u3*px + u1*py)/up + u2*pz;
    w =    -up*px +             u3*pz;
  } else if (u3<0.) {       // phi=0  teta=pi
    u = -u;
    w = -w;
  }
}

void RotateToLabFrame(double* dir, double* refdir) {
  RotateToLabFrame(dir[0], dir[1], dir[2], refdir[0], refdir[1], refdir[2]);
}


// it is assumed that track.fEkin > gamma-cut
void PerformBrem(Track& track, SimSBTables* theSBTable) {
  const double kPI            = 3.1415926535897932;
  const double kEMC2          = 0.510991;
  const double kHalfSqrt2EMC2 = kEMC2 * 0.7071067812;
  // sample energy transferred to the emitted gamma photon
  const double eGamma = theSBTable->SampleEnergyTransfer(track.fEkin,
                                                         track.fMatIndx,
                                                         Random::UniformRand(),
                                                         Random::UniformRand(),
                                                        Random::UniformRand());
 // insert the secondary gamma track into the stack
 Track& aTrack        = TrackStack::Instance().Insert();
 aTrack.fType         = 0;
 aTrack.fEkin         = eGamma;
 aTrack.fMatIndx      = track.fMatIndx;
 aTrack.fPosition[0]  = track.fPosition[0];
 aTrack.fPosition[1]  = track.fPosition[1];
 aTrack.fPosition[2]  = track.fPosition[2];
 aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
 aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
 aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
 //
 // compute emission direction (rough approximation in DPM by the mean)
 // and no deflection of the primary e-
 const double dum0    = kHalfSqrt2EMC2/(track.fEkin+kEMC2);
 const double cost    = std::max(-1.0, std::min(1.0, 1.0-dum0*dum0));
 const double sint    = std::sqrt((1.0+cost)*(1.0-cost));
 const double phi     = 2.0*kPI*Random::UniformRand();
 aTrack.fDirection[0] = sint*std::cos(phi);
 aTrack.fDirection[1] = sint*std::sin(phi);
 aTrack.fDirection[2] = cost;
 RotateToLabFrame(aTrack.fDirection, track.fDirection);
 // decrease the primary energy:
 track.fEkin = track.fEkin-eGamma;
}


void PerformAnnihilation(Track& track) {
  const double kPI      = 3.1415926535897932;
  const double kEMC2    = 0.510991;
  // isotropic direction
  const double cost = 1.0-2.0*Random::UniformRand();
  const double sint = std::sqrt((1.0-cost)*(1.0+cost));
  const double phi  = 2.0*kPI*Random::UniformRand();
  const double rx   = sint*cos(phi);
  const double ry   = sint*sin(phi);
  const double rz   = cost;

  Track& aTrack        = TrackStack::Instance().Insert();
  aTrack.fType         = 0;
  aTrack.fEkin         = kEMC2;
  aTrack.fMatIndx      = track.fMatIndx;
  aTrack.fPosition[0]  = track.fPosition[0];
  aTrack.fPosition[1]  = track.fPosition[1];
  aTrack.fPosition[2]  = track.fPosition[2];
  aTrack.fBoxIndx[0]   = track.fBoxIndx[0];
  aTrack.fBoxIndx[1]   = track.fBoxIndx[1];
  aTrack.fBoxIndx[2]   = track.fBoxIndx[2];
  aTrack.fDirection[0] = rx;
  aTrack.fDirection[1] = ry;
  aTrack.fDirection[2] = rz;

  Track& aTrack1        = TrackStack::Instance().Insert();
  aTrack1.fType         = 0;
  aTrack1.fEkin         = kEMC2;
  aTrack1.fMatIndx      = track.fMatIndx;
  aTrack1.fPosition[0]  = track.fPosition[0];
  aTrack1.fPosition[1]  = track.fPosition[1];
  aTrack1.fPosition[2]  = track.fPosition[2];
  aTrack1.fBoxIndx[0]   = track.fBoxIndx[0];
  aTrack1.fBoxIndx[1]   = track.fBoxIndx[1];
  aTrack1.fBoxIndx[2]   = track.fBoxIndx[2];
  aTrack1.fDirection[0] = -rx;
  aTrack1.fDirection[1] = -ry;
  aTrack1.fDirection[2] = -rz;
}
