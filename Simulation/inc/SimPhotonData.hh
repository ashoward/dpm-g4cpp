#ifndef SimPhotonData_HH
#define SimPhotonData_HH

//
// M. Novak: 2021
//
// Top level photon related data collection used during the simulation.
//
// It is assumed, that the corresponding data has been generated and written
// into files before the simulation. Therefore, this class will load and provide
// at run-time all the data that are required to perform the photon related part
// of the simulation.


#include <string>

class SimIMFPMaxPhoton;
class SimIMFPPhoton;

class SimKNTables;

class SimPhotonData {

public:

  SimPhotonData();
 ~SimPhotonData();

  // loads all photon related data
  void Load(const std::string& dataDir, int verbose=0);

  // for description see the member docs.
  SimIMFPMaxPhoton*    GetIMFPTotalMax () const  { return fIMFPTotalMax; }
  SimIMFPPhoton*       GetIMFPTotal()     const  { return fIMFPTotal;    }
  SimIMFPPhoton*       GetIMFPCompton()   const  { return fIMFPCompton;  }
  SimIMFPPhoton*       GetIMFPPairProd()  const  { return fIMFPPairProd; }

  SimKNTables*         GetTheKNTables()   const  { return fTheKNTables;  }


private:

  // global maximum of the total IMFP over the entire geometry (i.e. across the
  // different materials) at each individual discrete photon energy in [1/mm]
  // units.
  SimIMFPMaxPhoton*    fIMFPTotalMax;

  // total (sum of Compton, Pair-production and Photoelectric) inverse MFP for
  // each of the individual material and scalled (divided) by the material density
  // in [g/cm3]. So it will have the proper [1/mm] units only after scaling back
  // (multiplying) by a given (voxel) material density in [g/cm3] units.
  SimIMFPPhoton*       fIMFPTotal;

  // IMFP for Compton scattering for each of the individual material and scalled
  // (divided) by the material density in [g/cm3]. So it will have the proper
  // [1/mm] units only after scaling back (multiplying) by a given (voxel) material
  // density in [g/cm3] units.
  SimIMFPPhoton*       fIMFPCompton;

  // IMFP for Pair-production for each of the individual material and scalled
  // (divided) by the material density in [g/cm3]. So it will have the proper
  // [1/mm] units only after scaling back (multiplying) by a given (voxel) material
  // density in [g/cm3] units.
  SimIMFPPhoton*       fIMFPPairProd;


  // sampling tables to provide rejection free sampling of the (reduced) energy
  // transferred to the secondary electron in Compton scattering according to
  // the Klein-Nishina DCS
  SimKNTables*         fTheKNTables;
};

#endif // SimPhotonData_HH
