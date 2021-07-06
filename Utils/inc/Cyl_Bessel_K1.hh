

#ifndef CYL_BESSEL_K1
#define CYL_BESSEL_K1



// The irregular modified cylindrical Bessel function of first order: K1(x) for
// x > 0 taken from GSL

/* specfunc/bessel_K1.c
 * Author:  G. Jungman
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


namespace GSL {

#define 	DBL_MIN          2.2250738585072014e-308
#define   DBL_MAX          1.7976931348623157e+308
#define   LOG_DBL_MAX      7.0978271289338397e+02
#define   LOG_DBL_MIN     (-7.0839641853226408e+02)
#define   SQRT_DBL_MAX     1.3407807929942596e+154
#define   SQRT_DBL_MIN     1.4916681462400413e-154
#define 	DBL_EPSILON      2.2204460492503131e-16
#define   SQRT_DBL_EPSILON 1.4901161193847656e-08
#define 	SUCCESS          0

#define 	ERROR_SELECT_2(a, b)   ((a) != SUCCESS ? (a) : ((b) != SUCCESS ? (b) : SUCCESS))
#define   SIGN(x)    ((x) >= 0.0 ? 1 : -1)


struct sf_result {
  double val;
  double err;
};


// -------------------------------------------------------------------------- //
// "cheb_eval.c"
struct cheb_series {
  double*  c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};

static inline int cheb_eval_e(const cheb_series* cs, const double x, sf_result* result);



// -------------------------------------------------------------------------- //
// "bessel_I1.c"

#define ROOT_EIGHT (2.0*1.4142135623730951)

static double bi1_data[11] = {
  -0.001971713261099859,
   0.407348876675464810,
   0.034838994299959456,
   0.001545394556300123,
   0.000041888521098377,
   0.000000764902676483,
   0.000000010042493924,
   0.000000000099322077,
   0.000000000000766380,
   0.000000000000004741,
   0.000000000000000024
};
static cheb_series bi1_cs = {
  bi1_data,
  10,
  -1, 1,
  10
};

static double ai1_data[21] = {
  -0.02846744181881479,
  -0.01922953231443221,
  -0.00061151858579437,
  -0.00002069971253350,
   0.00000858561914581,
   0.00000104949824671,
  -0.00000029183389184,
  -0.00000001559378146,
   0.00000001318012367,
  -0.00000000144842341,
  -0.00000000029085122,
   0.00000000012663889,
  -0.00000000001664947,
  -0.00000000000166665,
   0.00000000000124260,
  -0.00000000000027315,
   0.00000000000002023,
   0.00000000000000730,
  -0.00000000000000333,
   0.00000000000000071,
  -0.00000000000000006
};
static cheb_series ai1_cs = {
  ai1_data,
  20,
  -1, 1,
  11
};

static double ai12_data[22] = {
   0.02857623501828014,
  -0.00976109749136147,
  -0.00011058893876263,
  -0.00000388256480887,
  -0.00000025122362377,
  -0.00000002631468847,
  -0.00000000383538039,
  -0.00000000055897433,
  -0.00000000001897495,
   0.00000000003252602,
   0.00000000001412580,
   0.00000000000203564,
  -0.00000000000071985,
  -0.00000000000040836,
  -0.00000000000002101,
   0.00000000000004273,
   0.00000000000001041,
  -0.00000000000000382,
  -0.00000000000000186,
   0.00000000000000033,
   0.00000000000000028,
  -0.00000000000000003
};

static cheb_series ai12_cs = {
  ai12_data,
  21,
  -1, 1,
  9
};


int sf_bessel_I1_scaled_e(const double x, sf_result* result);
int sf_bessel_I1_e(const double x, sf_result* result);





static double bk1_data[11] = {
   0.0253002273389477705,
  -0.3531559607765448760,
  -0.1226111808226571480,
  -0.0069757238596398643,
  -0.0001730288957513052,
  -0.0000024334061415659,
  -0.0000000221338763073,
  -0.0000000001411488392,
  -0.0000000000006666901,
  -0.0000000000000024274,
  -0.0000000000000000070
};

static cheb_series bk1_cs = {
  bk1_data,
  10,
  -1, 1,
  8
};

static double ak1_data[17] = {
   0.27443134069738830,
   0.07571989953199368,
  -0.00144105155647540,
   0.00006650116955125,
  -0.00000436998470952,
   0.00000035402774997,
  -0.00000003311163779,
   0.00000000344597758,
  -0.00000000038989323,
   0.00000000004720819,
  -0.00000000000604783,
   0.00000000000081284,
  -0.00000000000011386,
   0.00000000000001654,
  -0.00000000000000248,
   0.00000000000000038,
  -0.00000000000000006
};
static cheb_series ak1_cs = {
  ak1_data,
  16,
  -1, 1,
  9
};

static double ak12_data[14] = {
   0.06379308343739001,
   0.02832887813049721,
  -0.00024753706739052,
   0.00000577197245160,
  -0.00000020689392195,
   0.00000000973998344,
  -0.00000000055853361,
   0.00000000003732996,
  -0.00000000000282505,
   0.00000000000023720,
  -0.00000000000002176,
   0.00000000000000215,
  -0.00000000000000022,
   0.00000000000000002
};

static cheb_series ak12_cs = {
  ak12_data,
  13,
  -1, 1,
  7
};


int sf_bessel_K1_scaled_e(const double x, sf_result* result);
int sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, sf_result* result);
int sf_bessel_K1_e(const double x, sf_result* result);



// irregular modified cylindrical Bessel funtion of the first kind
double Ir_mod_cyl_Bessel_K1(const double x);

} // namespace GSL

#endif // CYL_BESSEL_K1
