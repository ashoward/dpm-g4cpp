
#ifndef G4Setup_HH
#define G4Setup_HH

#include <vector>
#include <string>

class G4Material;

// This builds a fake Geant4 geometry having the given list of NIST materials
// in the geometry with the given production threshold. The production threshold
// given in energy is chnage to the corresponding length before construction.

// A vector fo G4 NIST material names and the production cut value
// in energy {MeV}
void   FakeG4Setup(const std::vector<std::string>& g4Materials, double electronCutInEnergy, double gammaCutInEnergy, int verbose=0);

double ConvertCutToRange(double electronCutInEnergy, const G4Material* mat);
double ConvertCutToAbsLength(double gammaCutInEnergy, const G4Material* theMaterial);
double ComputeApproximatedDEDX(double ekin, double zet);
double ComputeApproximatedAbsXsecPerAtom(double ekin, double zet);


#endif // G4Setup_HH
