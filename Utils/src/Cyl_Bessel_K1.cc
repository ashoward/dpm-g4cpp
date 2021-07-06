
#include "Cyl_Bessel_K1.hh"

#include <iostream>
#include <cmath>

namespace GSL {

static inline int cheb_eval_e(const cheb_series* cs, const double x, sf_result* result) {
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += std::fabs(y2*temp) + std::fabs(dd) + std::fabs(cs->c[j]);
    dd = temp;
  }

  {
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += std::fabs(y*temp) + std::fabs(dd) + 0.5 * std::fabs(cs->c[0]);
  }

  result->val = d;
  result->err = DBL_EPSILON * e + std::fabs(cs->c[cs->order]);

  return SUCCESS;
}



int sf_bessel_I1_scaled_e(const double x, sf_result* result) {
  const double xmin    = 2.0 * DBL_MIN;
  const double x_small = ROOT_EIGHT * SQRT_DBL_EPSILON;
  const double y       = std::fabs(x);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return SUCCESS;
  }
  else if(y < xmin) {
    std::cerr << " UNDERFLOW_ERROR Error in sf_bessel_I1_scaled_e !" << std::endl;
    //UNDERFLOW_ERROR(result);
  }
  else if(y < x_small) {
    result->val = 0.5*x;
    result->err = 0.0;
    return SUCCESS;
  }
  else if(y <= 3.0) {
    const double ey = std::exp(-y);
    sf_result c;
    cheb_eval_e(&bi1_cs, y*y/4.5-1.0, &c);
    result->val  = x * ey * (0.875 + c.val);
    result->err  = ey * c.err + y * DBL_EPSILON * std::fabs(result->val);
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else if(y <= 8.0) {
    const double sy = sqrt(y);
    sf_result c;
    double b;
    double s;
    cheb_eval_e(&ai1_cs, (48.0/y-11.0)/5.0, &c);
    b = (0.375 + c.val) / sy;
    s = (x > 0.0 ? 1.0 : -1.0);
    result->val  = s * b;
    result->err  = c.err / sy;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else {
    const double sy = std::sqrt(y);
    sf_result c;
    double b;
    double s;
    cheb_eval_e(&ai12_cs, 16.0/y-1.0, &c);
    b = (0.375 + c.val) / sy;
    s = (x > 0.0 ? 1.0 : -1.0);
    result->val  = s * b;
    result->err  = c.err / sy;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  return SUCCESS;
}


int sf_bessel_I1_e(const double x, sf_result* result) {
  const double xmin    = 2.0 * DBL_MIN;
  const double x_small = ROOT_EIGHT * SQRT_DBL_EPSILON;
  const double y = std::fabs(x);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return SUCCESS;
  }
  else if(y < xmin) {
    std::cerr << " UNDERFLOW_ERROR Error in sf_bessel_I1_e !" << std::endl;
    //UNDERFLOW_ERROR(result);
  }
  else if(y < x_small) {
    result->val = 0.5*x;
    result->err = 0.0;
    return SUCCESS;
  }
  else if(y <= 3.0) {
    sf_result c;
    cheb_eval_e(&bi1_cs, y*y/4.5-1.0, &c);
    result->val  = x * (0.875 + c.val);
    result->err  = y * c.err;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else if(y < LOG_DBL_MAX) {
    const double ey = std::exp(y);
    sf_result I1_scaled;
    sf_bessel_I1_scaled_e(x, &I1_scaled);
    result->val  = ey * I1_scaled.val;
    result->err  = ey * I1_scaled.err + y * DBL_EPSILON * std::fabs(result->val);
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else {
    std::cerr << " OVERFLOW_ERROR Error in sf_bessel_I1_e !" << std::endl;
    //OVERFLOW_ERROR(result);
  }
  return SUCCESS;
}



int sf_bessel_K1_scaled_e(const double x, sf_result* result) {
  if(x <= 0.0) {
    std::cerr << " DOMAIN_ERROR Error in sf_bessel_K1_scaled_e: x = " << x << " must be > 0 !" << std::endl;
    //DOMAIN_ERROR(result);
  }
  else if(x < 2.0*DBL_MIN) {
    std::cerr << " OVERFLOW_ERROR Error in sf_bessel_K1_scaled_e !" << std::endl;
//    OVERFLOW_ERROR(result);
  }
  else if(x <= 2.0) {
    const double lx = std::log(x);
    const double ex = std::exp(x);
    int stat_I1;
    sf_result I1;
    sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = sf_bessel_I1_e(x, &I1);
    result->val  = ex * ((lx-M_LN2)*I1.val + (0.75 + c.val)/x);
    result->err  = ex * (c.err/x + std::fabs(lx)*I1.err);
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return stat_I1;
  }
  else if(x <= 8.0) {
    const double sx = std::sqrt(x);
    sf_result c;
    cheb_eval_e(&ak1_cs, (16.0/x-5.0)/3.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else {
    const double sx = std::sqrt(x);
    sf_result c;
    cheb_eval_e(&ak12_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  return SUCCESS;
}


int sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, sf_result* result) {
  const double ay  = std::fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = std::fabs(dy * std::exp(x));
    return SUCCESS;
  }
  else if(   ( x < 0.5*LOG_DBL_MAX   &&   x > 0.5*LOG_DBL_MIN)
          && (ay < 0.8*SQRT_DBL_MAX  &&  ay > 1.2*SQRT_DBL_MIN)
    ) {
    double ex = std::exp(x);
    result->val  = y * ex;
    result->err  = ex * (std::fabs(dy) + std::fabs(y*dx));
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return SUCCESS;
  }
  else {
    const double ly  = std::log(ay);
    const double lnr = x + ly;

    if(lnr > LOG_DBL_MAX - 0.01) {
      std::cerr << " OVERFLOW_ERROR Error in sf_exp_mult_err_e !" << std::endl;
      //OVERFLOW_ERROR(result);
    }
    else if(lnr < LOG_DBL_MIN + 0.01) {
      std::cerr << " UNDERFLOW_ERROR Error in sf_exp_mult_err_e !" << std::endl;
      //UNDERFLOW_ERROR(result);
    }
    else {
      const double sy  = SIGN(y);
      const double M   = std::floor(x);
      const double N   = std::floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      const double eMN = std::exp(M+N);
      const double eab = std::exp(a+b);
      result->val  = sy * eMN * eab;
      result->err  = eMN * eab * 2.0*DBL_EPSILON;
      result->err += eMN * eab * std::fabs(dy/y);
      result->err += eMN * eab * std::fabs(dx);
      return SUCCESS;
    }
  }
  return SUCCESS;
}


int sf_bessel_K1_e(const double x, sf_result* result) {
  if (x <= 0.0) {
    std::cerr << " DOMAIN_ERROR Error in sf_bessel_K1_e(x): x = " << x << " must be > 0 !" << std::endl;
  }
  else if(x < 2.0*DBL_MIN) {
    std::cerr << " OVERFLOW_ERROR Error in sf_bessel_K1_e(x): x = " << x << " must be > 2.0*DBL_MIN !" << std::endl;
  }
  else if(x <= 2.0) {
    const double lx = std::log(x);
    int stat_I1;
    sf_result I1;
    sf_result c;
    cheb_eval_e(&bk1_cs, 0.5*x*x-1.0, &c);
    stat_I1 = sf_bessel_I1_e(x, &I1);
    result->val  = (lx-M_LN2)*I1.val + (0.75 + c.val)/x;
    result->err  = c.err/x + std::fabs(lx)*I1.err;
    result->err += 2.0 * DBL_EPSILON * std::fabs(result->val);
    return stat_I1;
  }
  else {
    sf_result K1_scaled;
    int stat_K1 = sf_bessel_K1_scaled_e(x, &K1_scaled);
    int stat_e  = sf_exp_mult_err_e(-x, 0.0,
                                           K1_scaled.val, K1_scaled.err,
                                           result);
    result->err = fabs(result->val) * (DBL_EPSILON*fabs(x) + K1_scaled.err/K1_scaled.val);
    return ERROR_SELECT_2(stat_e, stat_K1);
  }
  return SUCCESS;
}


double Ir_mod_cyl_Bessel_K1(const double x)
{
  sf_result result;
  sf_bessel_K1_e(x, &result);
  return result.val;
}


} // namespace GSL
