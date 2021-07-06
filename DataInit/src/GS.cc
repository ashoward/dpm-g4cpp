
#include "GS.hh"


#include "Cyl_Bessel_K1.hh"

#include <cmath>
#include <iostream>

GS& GS::Instance() {
  static GS gs;
  return gs;
}


double GS::GetGSCoef(double nel, double parScreening, int i, bool isRemoveOnlyNoScattering) {
  const double expNel = std::exp(-nel);
  // recompute from the begining
  if (std::abs(nel-fCurNel)>1.0E-14 || std::abs(parScreening-fCurParScreening)>1.0E-14 ) {
    //ResetGSCoefs();
    fCurNel          = nel;
    fCurParScreening = parScreening;
    fCurMaxEl        = -1;
  }
  while ( i>fCurMaxEl ) {
    for (int j=0; j<fAdvanceBy; ++j) {
      int l = j + fCurMaxEl +1;
      // the first transprt coefficient is zero and the corresponding xi_0 = 1.0;
      double thisXi = 1.0;
      if (l>0) {
        // compute the l-th transprt coefficient G_l using Kawrakow's approximation Eq.(42)
        const double el1  = l*(l+1.);
        const double yel = 2.*std::sqrt(el1*parScreening);
        double    theGel = 1.0;
        // K1(yel=30) ~ 1E-14 i.e. already near zero so G_l ~ 1.0 as init. above
        if (yel<30.0) {
          const double besselK1 = GSL::Ir_mod_cyl_Bessel_K1(yel);
          const double dum      = GetNthHarmonicNumebr(l) - 0.5*std::log(el1);
          theGel = 1. - yel*besselK1*( 1. + 0.5*yel*yel*dum );
        }
        // compute xi_l
        double xi1 = std::exp(-nel*theGel) - expNel;
        double xi2 = 1. - expNel;
        if (!isRemoveOnlyNoScattering) {
          xi1 -= expNel*nel*(1.-theGel);
          xi2 -= nel*expNel;
        }
        thisXi = xi1/xi2;
      }
      if (l+1>fMaxEl) {
        std::cout << " ***Error in GS::GetGSCoef: increase fMaxEl = " << fMaxEl << " and repeate current computations."<< std::endl;
      }
      fGSCoefs[l] = thisXi;
    }
    fCurMaxEl += fAdvanceBy;
  }
  return fGSCoefs[i];
}



double GS::ComputeOptimalTransformationParameter(double nel, double parScreening, double acc) {
  const int kMinEl   = 1000;
  double sumAlpha    = 0.0;
  double sumBetha    = 0.0;
  bool   isDoneAlpha = false;
  bool   isDoneBetha = false;
  // computing the optimal transformation parameter requires 3 onsecutive xi_i
  // GS coefficientsat each `el` (i.e. `i`) -value
  double xi0 = GS::Instance().GetGSCoef(nel, parScreening, 0);
  double xi1 = GS::Instance().GetGSCoef(nel, parScreening, 1);
  double i0  = 0.0;
  while (! (isDoneBetha && isDoneAlpha) ) {
    double i1  = i0 + 1.0;
    double i2  = i1 + 1.0;
    double xi2 = GS::Instance().GetGSCoef(nel, parScreening, i2);
    // \betha = (l+1) xi_i x_{i+1}
    double iBetha = i1*xi0*xi1;
    if (i0==0.0) {
      sumBetha      = iBetha;
    } else if (!isDoneBetha) {
      isDoneBetha   = i0 > kMinEl && iBetha<=acc*sumBetha;
      sumBetha     += iBetha;
    }
    // \alpha = ... Eq.(38)
    const double dum0 = (1.5*i0+0.0625/(i0+1.5)+0.0625/(i0-0.5)+0.75)*xi0;
    const double dum1 = i1*i2/(2.0*i0+3.0)*xi2 - 2.0*i1*xi1;
    double iAlpha = xi0*( dum0 + dum1);
    if (i0==0.0) {
      sumAlpha     = iAlpha;
    } else if (!isDoneAlpha) {
      isDoneAlpha  = i0 > kMinEl && iAlpha<acc*sumAlpha;
      sumAlpha    += iAlpha;
    }
//    std::cout << i0 << " "<<sumAlpha << " " << sumBetha << " iBetha = " << iBetha<< " " << xi0 << " " << xi1<<std::endl;
    xi0 = xi1;
    xi1 = xi2;
    i0 += 1.0;
  }
  //sumBetha = std::max(sumBetha, 1.0E-20);
  const double kZero = 1.0E-14;
  const double dum = (sumBetha>kZero) ? 0.25*sumAlpha/sumBetha : 1./kZero;
  //const double dum = 0.25*sumAlpha/sumBetha;
  return dum + std::sqrt(dum*(dum+1.0));
}


double GS::ComputeTransformedGSDistribution(double nel, double parScreening, double parTransf, double uVale) {
  const double kAcc = 1.0E-10;
  const double kDum = 1./(1.-uVale+parTransf);
  // if u=0.0 then mu(u) = +1.0 and since P_el(+1) = x^ell
  // if u=1.0 then mu(u) = -1.0 and since P_el(-1) = x^ell
  // the `mu(u)` value (i.e. the inverse transform) as given by Eq. (34)
  const double x = 1. - 2.*parTransf*uVale*kDum;
  // compute the sum of the Legendre polynomials as given in the summation of Eq.(41)
  // NOTE: using Bonnet’s recursion formula P_{n+1}(x) =  [x(2n+1)P_{n}(x) - nP_{n-1}(x)]/(n+1)
  double p0    = 1.;  // P_0(x) = 1.0
  double p1    = x;   // P_1(x) = x with x = mu(u)
  // sum of the n = 0  and n = 1
  const double term0 = 0.5*p0*GS::GetGSCoef(nel, parScreening, 0);
  const double term1 = 1.5*p1*GS::GetGSCoef(nel, parScreening, 1);
  double sum   = term0+term1;
  // compute the n+1 -th term starting with n+1=2
  double n1    = 2.;  // n+1 = 2 (n=1 and n-1=0)
  bool   isDone = false;
  do {
    const double n  = n1 - 1.;
    // the n+1 -th Legendre polynomial
    const double p2 = (x*(2.*n+1.)*p1 - n*p0)/n1;
    // the n+1 -th term of the sum
    const double iTerm = (n1+0.5)*p2*GS::GetGSCoef(nel, parScreening, n1);
    // is convergred
    isDone =  std::abs(iTerm) < kAcc*std::abs(sum);
    // add this n+1 -th term to the sum
    sum   += iTerm;
    // update n1, p0 and p1 fo rthe next (n+2 -th) step
    n1 += 1.;
    p0  = p1;
    p1  = p2;
  } while (!isDone);
  // the final multiplicative factor in Eq.(41) and return
  return 2.*parTransf*(1.+parTransf)*kDum*kDum*sum;
}



GS::GS () {
  fCurMaxEl  = -1;
  fAdvanceBy = 20;
  fMaxEl     = 100000;

  fCurNel          = -1.0;
  fCurParScreening = -1.0;

  fGSCoefs.resize(fMaxEl, 0.0);
}


double GS::GetNthHarmonicNumebr(int n) {
    const int    kLimit = 8;
    const double kGamma = 0.577215664901; // Euler number
    const double kC2    = 1.0/12.0;
    const double kC3    = 1.0/120.0;
    double res = 1.0;
    if (n<kLimit) {
      for (int i=2; i<n+1; ++i) {
         res += 1./i;
       }
    } else {
      double inn = 1./(n*n);
      res = std::log(n) + kGamma + 0.5/n -kC2*inn - kC3*inn*inn;
    }
    return res - kGamma;
  }
