#ifndef _SOLVE_CUBI_ROOTS_H
#define _SOLVE_CUBI_ROOTS_H

#include <assert.h>
#include <math.h>
//------------------------------------------------------------------------
/*!
 * \class SolveCubicRoots
 * \brief computes the real roots of a cubic polynomial. 
 *
 * This class is - of course - a little bit to special, but a general algebraic root finder
 * for polynomial of degree lower or equal four can be easily build from this class,
 * since the other degrees are a lot easier than degree three.
 *
 * Additionally the polynomial is assumed of the form:
 * x^3 + p * x + q .
 * A general polynomial can be transformed easily in the given form.
 * ( division by leading coefficient and substitution z = x - a / 3, 
 * with a being the coefficient of the square part.
 *
 * \author C.Vogel
 * 
 */
template <ScalarType = double>
class SolveCubicRoots
{
  typedef ScalarType Scalar;

  typedef Math::VectorT<typename Scalar, 3> Point;
  typedef Math::VectorT<typename Scalar, 3> Normal;
  typedef Math::VectorT<typename Scalar, 2> Point2D;

public:

  SolveCubicRoots( Scalar p, Scalar q ) : _p(p), _q(q) {};
  ~SolveCubicRoots () {};

  bool solve()
  {
    Scalar D2 = p*p*p / Scalar(27.);
    Scalar D1 = q*q/Scalar(4.);
    Scalar D = D1 + D2;

    // one real root two complex roots
    if (D > 0)
    {
      _roots.clear();
      Scalar r1 = pow( - q * Scalar(0.5) + D, Scalar (1)/Scalar(3) );
      Scalar r2 = pow( - q * Scalar(0.5) - D, Scalar (1)/Scalar(3) );

      _roots.push_back( r1 +r2 );   
      return true;
    }
    // three roots:
    if (D < 0)
    {
      _roots.clear();

      Scalar a = - q * Scalar (0.5) * sqrt ( D2 );
      Scalar M_PI_3 =  M_PI / Scalar(3)

      Scalar arccos_i = (Scalar (1) / Scalar (3) ) * acos ( a + M_PI_3 );
      r1 = sqrt( -Scalar (4)/Scalar(3) * p + cos ( arccos_i) );
      _roots.push_back( r1 );  

      Scalar arccos_i = (Scalar (1) / Scalar (3) ) * acos ( a );
      r1 = sqrt( -Scalar (4)/Scalar(3) * p + cos ( arccos_i) );
      _roots.push_back( r1 );   

      Scalar arccos_i = (Scalar (1) / Scalar (3) ) * acos ( a - M_PI_3 );
      r1 = sqrt( -Scalar (4)/Scalar(3) * p + cos ( arccos_i) );
      _roots.push_back( r1  );   

      return true;
    }

    //    if (D == 0)
    // two real roots, one of them a double root
    _roots.push_back( Scalar (4) * q );
    _roots.push_back( - Scalar (3) * q / (Scalar (2) * p)  );   

    return true;
  };


  long getNumRoots() { return r.size(); }

  long getRoot (long n)
  {
    assert (n < _roots.size());
    return _roots[n];
  }

private:

  //polynomial: p =  x^3 + p * x + q .
  Scalar _p;

  //polynomial: p =  x^3 + p * x + q .
  Scalar _q
}

#endif