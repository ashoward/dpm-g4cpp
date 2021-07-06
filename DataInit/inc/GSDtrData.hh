
#ifndef GSDtrData_HH
#define GSDtrData_HH


#include <vector>

// The GS angular distributions, more precisely, a representation of the
// transformd q1+ GS pdf-s that make possible the efficient (rejection free
// sampling) of the corresponding `mu = cos(theta)` angular deflections over
// a given s-path lenght.
//

struct GSDtrData {

  GSDtrData (double eMin, double eMax, int numEkin=128, int numCumVals=51) {
    //
    // set the kinetic energy grid related variables
    fLogEKinMin   = std::log(eMin);
    fNumEKin      = numEkin;
    double dlEkin = std::log(eMax/eMin)/(numEkin-1);
    fInvLogDelta  = 1./dlEkin;
    // generate the kinetic energy grid
    fEKinGrid.resize(numEkin, 0.);
    fEKinGrid[0]  = eMin;
    int ie = 0;
    for (; ie<numEkin-1; ++ie) {
      fEKinGrid[ie] = std::exp(fLogEKinMin+ie*dlEkin);
    }
    fEKinGrid[ie] = eMax;
    //
    // set the (common) cumulative grid related variables
    fNumCumData   = numCumVals;
    fInvDeltaCum  = (numCumVals-1);  // i.e. 1./[(1-0)/(N-1)]
    //
    // set the data vector size
    fGSDtrData.resize(numEkin);
  }

  // at each kinetic energy grid point (that determines s/lambda_e and parScreening)
  // a GS q1+ pdf distribution is computed over a given `u` transformed variable
  // grid. Then this is used to generate an approximate inversion of the corresponding
  // cumulative, such that the cumulative values are equaly spaced over the [0,1]
  // inteval. This makes possible a fast, rejection free sampling of the transformed
  // variable `u`, and the transformation variable `a` (that belongs to a given s/lambda_e
  // and parScreening) can be used to transform it back to get \mu = \cos(\theta) at
  // this E_i.
  // At run time at a given energy `E`, first the discrete kinetic energy index `i`
  // is computed such that E_{i} <= E < E_{i+1}, then statistical interpolation
  // (corresponding to a linera interpolation in log-energy scale) is used to select
  // either the distribution at the 'E_i' or E_{i+1}, and a mu=cos(theta) value
  // is drawn according to these GS F2+ pdf.

  struct OneCumData {
    double fVarU;     // the transfomred variable `u` that corresponds this cumulative value
    double fParmA;    // interpolation parameters for the approximate inversion of the cumulative
    double fParmB;
  };
  struct OneGSDtr {
    double  fTransformParam;          // transformation paraneter `a` used to generate the corresponding smooth pdf
    std::vector<OneCumData> fCumData; // the representation of the inverse cumulative of this pdf
  };




  // the kinetic energy grid: at each e-kin point,  transformed GS dtr (q1+ pdf)
  //                          is computed and a representation of the inverse
  //                          cumulative is computed. This later is used at run-time
  //                          for sampling the `mu(u) = cos(theta)` deflection at
  //                          the given E_i kinetic energy point.
  int      fNumEKin;
  double   fLogEKinMin;
  double   fInvLogDelta;
  std::vector<double>   fEKinGrid;
  // the (common) discrete cumulative grid: a discrete cumulative, equaly spaced
  //                          on [0,1]. The representation of an individual inverse
  //                          cumulative contains the `u` grid, tat corresponds to
  //                          these cumulative values.
  int      fNumCumData;
  double   fInvDeltaCum;


  // the GS q1+ distributions over the kinetic energy grid: each of them corresponds
  //                          the given E_i kinetic energy and contains a representation
  //                          of the inverse cumulative of the corresponding tranformed
  //                          q1+ GS pdf.
  std::vector<OneGSDtr> fGSDtrData;
};


#endif // GSDtrData_HH
