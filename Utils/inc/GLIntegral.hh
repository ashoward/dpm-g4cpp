
#ifndef GLINTEGRAL_H
#define GLINTEGRAL_H

#include <vector>

/**
 * @brief   Numerical integral utility.
 * @class   GLIntegral
 * @author  M Novak, A Ribon
 * @date    january 2016
 *
 * Utility class to generates abscissas \f$\{x_i\}_N\f$ and weights \f$\{w_i\}_N\f$ of an N-point Gauss-Legendre
 * quadrature for perform numerical integration of a function \f$f(x)\f$ between \f$x_{min}\f$ and \f$x_{max}\f$
 * \f[
 *   \int_{x_{min}}^{x_{max}} f(x)\mathrm{d}x \approx \sum_{i=0}^{N-1} w_{i} f(x_i)
 * \f]
 */

class GLIntegral {
public:
  /**
  * @name Constructor, destructor:
  */
  //@{
  /**
    * @brief Constructor.
    *
    * Abscissas and weights are generated at construction so the constructed object is ready to use.
    *
    * @param[in] npoints  Defines number of points used in the Gauss-Legendre quadrature.
    * @param[in] xmin     Lower limit of the integral.
    * @param[in] xmax     Upper limit of the integral.
    */
  GLIntegral(int npoints, double xmin, double xmax);
  /** @brief Destructor. */
  ~GLIntegral() {}
  //@}

  /**
  * @name Public getters.
  */
  //@{
  /** @brief Public method to get the weight vector. Size of the vector is GLIntegral::fNPoints. */
  const std::vector<double> &GetWeights() const { return fWeights; }
  /** @brief Public method to get the abscissa vector. Size of the vector is GLIntegral::fNPoints. */
  const std::vector<double> &GetAbscissas() const { return fAbscissas; }
  //@}

private:
  /** @brief Intenal method to compute abscissas and weights. Used at construction. */
  void SetParameters();

private:
  /** @brief Order of appoximation i.e. number of abscissas and weights. */
  int fNPoints;
  /** @brief Lower limit of the integral. */
  double fXmin;
  /** @brief Upper limit of the integral. */
  double fXmax;
  /** @brief Container to store the generated weights. Size of the vector is GLIntegral::fNPoints. */
  std::vector<double> fWeights;
  /** @brief Container to store the generated abscissas. Size of the vector is GLIntegral::fNPoints. */
  std::vector<double> fAbscissas;
};


#endif // GLINTEGRAL_H
