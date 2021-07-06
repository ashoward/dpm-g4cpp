#ifndef SimDPMLike_HH
#define SimDPMLike_HH

//
// M. Novak: 2021
//
// A DMP like simulation of e-/e+ and gamma particles across a voxelised geometry.
// The main differences between DMP and this version:
// - some integrated data were extracted from Geant4 at the pre-init phase
// - other data have been computed and written into files during the pre-inint
//   phase, that make possible the rejection free sampling (i.e. without non-
//   deterministic loop counters) all of the e-/photon interactions including
//   sampling of the angular deflections in MSC according to Goudsmit-Saunderson
//   angular distributions.
//
// Note: that the geomerty implementation is rather simple since this is a
//       simulation prototype only but the simulation itself do not depend on
//       any speciality of the current geometry implementation.

class Geom;
class SimMaterialData;
class SimElectronData;
class SimPhotonData;

class Track;

//
// The top level function that performs the complete simulation of `nprimary`
// histories of electrons (`iselectron=true`) or photons (`iselectron`=false)
// hitting one of the predefined (voxelised) geometies perpendicularly with a
// primary energy of `eprimary` in [MeV]. The size of the geomerty voxles (boxes)
// is determined by the `lbox` input parameter in [mm] units.
// The additional `Material-`, `Electron-` and `PhotonData` input paraneters are supposed to
// be already loaded into the memory. The `geomIndex` is one of the pre-defied geometrical
// configuration (see more at `Geom.hh`).
//
// The resulted depth dose histogram of the simulation is written into the `hist.sim`
// file at the termination of the simulation.
//
void  Simulate(int nprimary, double eprimary, bool iselectron, double lbox, SimMaterialData& matData, SimElectronData& elData, SimPhotonData& phData, int geomIndex);



// Keeps tracking an e- till one of the following condition is reached:
// 1. a discrete brem interaction take place
// 2. a discrete ioni interaction take place
// 3. MSC hinge or end of MSC step take place
// 4. the electron energy drops below zero so stops
int   KeepTrackingElectron(SimElectronData& elData, SimMaterialData& matData, Geom& geom, Track& track, double& numElTr1MFP, double& numMollerMFP, double invMollerMFP, double& numBremMFP);

//
// Keeps tracking a photon till one of the following condition is reached:
// 1. pair-production take place
// 2. Compton sattering take place
// 3. Photoelectric absoprtion take place
// Unlike in case of electrons, this function also performs the interactions
// themself since the photon interactions are very simple in a DPM like simulation
void   KeepTrackingPhoton(SimPhotonData& phData, SimMaterialData& matData, Geom& geom, Track& track);

//
// A set of utility methods:
//

//
// Roate the direction [u,v,w] given in the scattering frame to the lab frame.
// Details: scattering is described relative to the [0,0,1] direction (i.e. scattering
// frame). Therefore, after the new direction is computed relative to this [0,0,1]
// original direction, the real original direction [u1,u2,u3] in the lab frame
// needs to be accounted and the final new direction, i.e. in the lab frame is
// computed.
void   RotateToLabFrame(double &u, double &v, double &w, double u1, double u2, double u3);
void   RotateToLabFrame(double* dir, double* refdir);

// Auxiliary funtion for simple e+e- annihilation
void   PerformAnnihilation(Track& track);

#endif // SimDPMLike_HH
