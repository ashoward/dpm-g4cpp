

#ifndef GS_HH
#define GS_HH

#include <vector>

class GS {
public:
  static GS& Instance();

  // Computes the GS series coefs (exp(-s/lambda_i)) related xi_i quantities.
  // These are the xi_i coefs given by Eq.(40) in my notes and used later in
  // Eqs.(38) and (39) to get the optimal screening parameters (smoothing
  // the pdf.) and to compute the pdf given by Eq.(15).
  //
  // Note, that each time a number of coefficients are pre-computed and as long
  // as they are fine they are just returned. Whenever further coefs. are needed
  // we recompute them
  //
  // Also note:in DPM, the real GS coefs, i.e. the exp(-s/lambda_i) = exp(-s/lambda_e G_i)
  // (or exp(-s/lambda_e G_i) - exp(s/lambda_e) if the zero scatteriing part is
  //  required to be removed) coefficients are computed in 'libelastic::kcoef'.
  // This is becasue in my case the q2+ GS pdf-s are computed i.e. the zero and
  // single scattering terms are not considered. Moreover, the normalisation
  // (division by 1 - exp(-s/lambda_e) - s/lambda_e exp(-s/lambda_e)) is already
  // included in these coefficients, that is down only in the 'libelastic::gs'
  // (called from 'libelastic::qsurf') with the appropriate division by 1 or 1-exp(-s/lambda_e)
  // depending if zero scattering has been (yes in the second case as usual) requested
  // to be removed.
  // This might be changed later depending on how zero and single scattering
  // cases are dealt during the simulation (it seems DPM includes single scatering
  // in the GS pdf-s and exludes only no-scattering. So I expect now that a random
  // sampling happens and angular deflection, i.e. any elastic scattering, is
  // considered only if 'urand > exp(-s/lambda_e) i.e. p(no-scateirng) = exp(-s/lambda_e)).
  //
  // According to Eq.(13), if we remove only the no-scattering, then the at-least
  // one scattering distribution is:
  //
  // sum_l (2l+1)/2 P_l(\mu) [e^{-s/\lambda_e G_l} - e^{-s/\lambda_e}]
  //
  // and by (since this above must give the integral of 1-exp(-s/lambda_e) since
  // the total Eq.(13) is a pdf i.e. integral 1) dividing with the
  // probability of having at least 1 elastic scatteirng along s, i.e. by
  // 1 - exp(-s/lambda_e) gives the
  //
  // sum_l (2l+1)/2 P_l(\mu) [e^{-s/\lambda_e G_l} - e^{-s/\lambda_e}]/[1 - exp(-s/lambda_e)]
  //
  // which is now a proper pdf (with integral of 1) of having at least 1 elastic
  // scatteirng along the path s, instead of Eq.(15). Note, that this normalisation
  // of the xi_i coefficientis not important when computing the optimal transformation
  // paraneters (smoothing this pdf) since only ratios of xi_i (squares) appears so
  //  any such normalisation would cancel each other.
  // (This is why this normalisation is done only in the `qsurf` computations by DPM
  //  but it can already be included here in the xi_i computation without a problem.)
  //
  //
  // nel = s/lambda_el and parScreening is the screening parameter
  double GetGSCoef(double nel, double parScreening, int i, bool isRemoveOnlyNoScattering=true);

  // relies on the above GetGSCoef() to obtain the GS series coefficients
  // x_i given by Eq.(40) and using them to compute the optimal GS angular pdf
  // transformation paraneter (i.e. the one that results in a transformed pdf
  // that is the closest to the uniform dfstribution) given by Eqs.(37-39).
  //
  // This optimal transformation parameter can then be used in Eq.(41) to compute
  // the transformed, flat pdf.
  //
  // nel = s/lambda_el and parScreening is the screening parameter
  double ComputeOptimalTransformationParameter(double nel, double parScreening, double acc=1.0E-10);

  // Computes the q2+ GS pdf value (as given by Eq.(41) at the given input parameters:
  // nel = s/lambda_e
  // parScreening = screening parameter of the screened Rutherford DCS
  // parTransf    = optimal transformation parameter (as given by GetOptimalTransformationParameter())
  // uValue       = value of the transformed variable in [0,1]
  //
  // NOTE: since the GS coefficients, xi_i (as well the optimal transformation parameter
  //       value `a` through this xi_ values) depends on a given s/lambda_e and screening paraneter
  //       value, and these xi_i values are computed and reused as long as the s/lambda and screening
  //       paraneter value do not change: IT IS STONGLY ADVISED to use this function with a fixed `nel`
  //       and `parScreening` while changing only the `uVale`! (otherwise, several re-computations of
  //       the xi_i GS series will happen that makes the coputation slower)
  double ComputeTransformedGSDistribution(double nel, double parScreening, double parTransf, double uVale);



private:

  GS ();

  // utils to obtain the n-th harmonic number - Euler-const as \sum_{i=1}^N 1/i - gamma
  double GetNthHarmonicNumebr(int n);



private:

  int fCurMaxEl;
  int fAdvanceBy;
  int fMaxEl;

  double fCurNel;
  double fCurParScreening;

  // container for x_i (x, lambda, eta) values computed/obtained by GetGSCoef
  std::vector<double> fGSCoefs;


};


#endif // GS_HH
