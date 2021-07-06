
#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

/**
 * @brief   Spline interpolation utility.
 * @class   Spline
 * @author  M Novak, A Ribon
 * @date    january 2016
 *
 * Utility class to performe natural cubic spline interpolation and constrained natural cubic spline interpolation.
 */

class Spline {
public:
  /**
  * @name Constructor, destructor:
  */
  //@{
  /**
    * @brief Constructor.
    *
    * @param[in] xdata          Pointer to the first x-value. Spline does not own the data.
    * @param[in] ydata          Pointer to the first y-value. Spline does not own the data.
    * @param[in] numdata        Number of data points from the first x,y values to be used to set up the interpolator.
    * @param[in] isconstrained  Flag to indicate if constrained cubic spline interpolation is requested.
    * @param[in] isacsecderis   Flag to indicate that a more accurate second derivative computations is need to use.
    *                           Considered only if isconstrained=false and numdata>4.
    */
  Spline(double *xdata, double *ydata, int numdata, bool isacsecderis = true, bool isconstrained = false);

  /** @brief Default constructor. Spline::SetUpSpline method needs to be invoked before using the object created by this
   *         default constructor.
   */
  Spline() {}

  /** @brief Destructor. Spline does not own the x and y data used to set up the interpolator. */
  ~Spline() {}
  //@}

  /**
  * @name Public setters/getters.
  */
  //@{
  /**
    * @brief Public method to reset an already existing interpolator object.
    *
    * @param[in] xdata          Pointer to the new first x-value. Spline does not own the data.
    * @param[in] ydata          Pointer to the new first y-value. Spline does not own the data.
    * @param[in] numdata        Number of data points from the first x,y values to be used to set up the interpolator.
    * @param[in] isconstrained  Flag to indicate if constrained cubic spline interpolation is requested.
    * @param[in] isacsecderis   Flag to indicate that a more accurate second derivative computations is need to use.
    *                           Considered only if isconstrained=false and numdata>4.
    */
  void SetUpSpline(double *xdata, double *ydata, int numdata, bool isacsecderis = true, bool isconstrained = false);

  /** @brief Public method to get an interpolated y-value at a given x value.
    *
    * A binary search will be performed on the x-data grid to find \f$i\f$ such that \f$x_i \leq val < x_{i+1}\f$ or
    * \f$x_{i=0}\f$ is used if \f$val \leq x_{0}\f$ and \f$x_{i=numdata-1}\f$ is used if \f$val \leq x_{numdata-1}\f$.
    *
    * @param[in] val The x-value where the interpolated y-value is required.
    * @return    The interpolated y-value at the given x-value.
    */
  double GetValueAt(double val);

  /** @brief Public method to get an interpolated y-value at a given x value without binary search.
    *
    * The binary search can be skipped if we know x data grid index \f$i\f$ such that \f$x_i \leq val < x_{i+1}\f$.
    *
    * @param[in] x    The x-value where the interpolated y-value is required.
    * @param[in] indx The x-value grid index such that \f$x_i \leq val < x_{i+1}\f$.
    * @return    The interpolated y-value at the given x-value.
    */
  double GetValueAt(double x, int indx);
  //@}

  //
  //  double SplineInterpolation(double x);
  //  double SplineInterpolation(double x, int indx);

  // should be public or building with root fails
  struct Parameters {
    double fA;
    double fB;
    double fC;
    double fD;
  };

private:
  /** @brief Internal method to compute natural cubic spline interpolation parameters at set up. */
  void SetParameters();
  /** @brief Internal method to compute constrained natural cubic spline interpolation parameters at set up. */
  void SetParametersConstrained();

  ///
  void FillSecondDerivatives();

  // data members
private:
  /** @brief Flag to indicate if constrained natural cubic spline interpolator is required. */
  bool fIsConstrained; // flag to set the type of the cubic spline
  /** @brief Flag to indicate if more accurate second derivative computation is required. Considered only if
             fIsConstrained = false and fNPoints>4.*/
  bool fIsAcSecDerivs; // flag to use more accurate second derivatives; only in case of non-constrained and N>4
  /** @brief Some index value for the binary search. It's value depends on if x-grid is increasing or decreasing. */
  int fDirection; // index depending on if x grid is increasing or decreasing
  /** @brief Index value of the highest x-value in the x-grid. Used in the binary search. */
  int fUpperm; // highest x value index
  /** @brief Index value of the lowest x-value in the x-grid. Used in the binary search. */
  int fLowerm; // lowest x value index
  /** @brief Number of x,y data points used to set up the interpolator. */
  int fNPoints;   // number of XY points
                  /** @brief Pointer to the first x-value. Spline::fNPoints x-value will be used for the interpolation.
                    *        Spline does not own these x data.
                    */
  double *fXdata; // starting point of X(fNPoints will be considered from this as X data points)
                  /** @brief Pointer to the first y-value. Spline::fNPoints y-value will be used for the interpolation.
                    *        Spline does not own these y data.
                    */
  double *fYdata; // starting point of Y(fNPoints will be considered from this as Y data points)

  /** @brief Parameters of the cubic spline interpolation. Computed at set up. Size of the vector is Spline::fNPoints.*/
  std::vector<Parameters> fData;

  ///
  //  std::vector<double> fSecondDerivatives;
};


#endif // SPLINE_H
