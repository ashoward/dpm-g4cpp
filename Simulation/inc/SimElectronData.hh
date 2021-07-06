#ifndef SimElectronData_HH
#define SimElectronData_HH

//
// M. Novak: 2021
//
// Top level electron/positron related data collection used during the simulation.
//
// It is assumed, that the corresponding data has been generated and written
// into files before the simulation. Therefore, this class will load and provide
// at run-time all the data that are required to perform the electron/positron
// related part of the simulation.


#include <string>

class SimITr1MFPElastic;
class SimMaxScatStrength;

class SimIMFPMoller;
class SimIMFPBrem;

class SimStoppingPower;

class SimMollerTables;
class SimSBTables;
class SimGSTables;

class SimMaterialData;

class SimElectronData {

public:

  SimElectronData();
 ~SimElectronData();

  // loads all electron related data
  void Load(const std::string& dataDir, int verbose=0);

  // getters without further descriptions (see the correspnding member doc.)
  SimITr1MFPElastic*   GetITr1MFPElastic()   const { return fITr1MFPElastic;  }
  SimMaxScatStrength*  GetMaxScatStrength()  const { return fMaxScatStrength; }
  SimIMFPMoller*       GetIMFPMoller()       const { return fIMFPMoller;      }
  SimIMFPBrem*         GetIMFPBrem()         const { return fIMFPBrem;        }
  SimStoppingPower*    GetDEDX()             const { return fDEDX;            }

  SimMollerTables*     GetTheMollerTables()  const { return fTheMollerTables; }
  SimSBTables*         GetTheSBTables()      const { return fTheSBTables;     }
  SimGSTables*         GetTheGSTables()      const { return fTheGSTables;     }


private:

  // inverse first transprt mean free path for elastic interaction only for
  // the reference material and scalled (divided) by the material density in [g/cm3]
  // So it will have the proper [1/mm] units only after scaling back (multiplying)
  // by a given material density in [g/cm3] units.
  SimITr1MFPElastic*   fITr1MFPElastic;

  // maximum scattering strength (K_1(E)) in MSC step allowed in the reference
  // material that is approximated as K_1(E) ~ S_max(E)/tr1-mfp(E) with S_max(E)
  // being the max-MSC-step length determined by the `slow`, `shigh` and `ecross`
  // parameters.
  // Note, data are loaded only for the reference material and not scalled by
  // the density (i.e. it has no units).
  SimMaxScatStrength*  fMaxScatStrength;

  // restriced IMFP (macroscopic cross section) for Moller interaction only for
  // the reference material and scalled (divided) by the material density in [g/cm3]
  // So it will have the proper [1/mm] units only after scaling back (multiplying)
  // by a given material density in [g/cm3] units.
  SimIMFPMoller*       fIMFPMoller;

  // restriced IMFP (macroscopic cross section) for Brem. interaction for each of
  // the individual material and scalled (divided) by the material density in [g/cm3]
  // So it will have the proper [1/mm] units only after scaling back (multiplying)
  // by a given material density in [g/cm3] units.
  SimIMFPBrem*         fIMFPBrem;

  // restriced total (radiative plus collisonal) stopping power for each of the
  // individual material and scalled (divided) by the material density in [g/cm3]
  // So it will have the proper [MeV/mm] units only after scaling back (multiplying)
  // by a given material density in [g/cm3] units.
  SimStoppingPower*    fDEDX;


  // sampling tables to provide rejection free sampling of the energy transferred
  // to the secondary electron in Moller interaction
  SimMollerTables*     fTheMollerTables;

  // sampling tables to provide rejection free sampling of the energy that is
  // transferred to the secondary gamma in electron bremsstrahlung according to
  // the Seltzer-Berger DCS
  SimSBTables*         fTheSBTables;

  // sampling tables to provide rejection free sampling of the angular deflection
  // due to multiple Coulomb scattering in the reference material along the
  // maximally allowed MSC step length that is a function of the energy and
  // determined by a (sigmoid-like) function of the `slow`, `shigh` and `ecross`
  // parameters.
  SimGSTables*         fTheGSTables;
};

#endif // SimElectronData_HH
