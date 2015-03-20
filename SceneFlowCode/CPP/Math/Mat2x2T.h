//=============================================================================
//
//  CLASS Mat2x2T
//
//=============================================================================

#ifndef MAT2X2_HH
#define MAT2X2_HH


//== INCLUDES =================================================================

#include "VectorT.h"
#include <vector>
#include <iostream>

//== CLASS DEFINITION =========================================================

namespace Math {
//-----------------------------------------------------------------------------
/*! \brief defines a 2x2 matrix. Several often used algorithms are implemented as well. 
 *
 * \author C.Vogel
*/
template<class Scalar>
class Mat2x2T
{
public:

  // COLUMN_WISE stored

  typedef VectorT<Scalar,2>    Vec2;
  
  /// Empty constructor - object remains uninitialized.
  Mat2x2T( void ) {};

  Mat2x2T( const Scalar a00, const Scalar a01,
        	 const Scalar a10, const Scalar a11 )
  {
    elem_[0][0] = a00; elem_[1][0] = a01;
    elem_[0][1] = a10; elem_[1][1] = a11;
  };


  Mat2x2T( const Vec2& col0 , const Vec2& col1 )
  {
    elem_[0] = col0; elem_[1] = col1;
  };


  /// Copy-Constructor
  Mat2x2T( const Mat2x2T<Scalar> &p )
  {
    elem_[0] = p.elem_[0];
    elem_[1] = p.elem_[1];
  };
  
  /// Destructor
  ~Mat2x2T( void ) {};

  
  //
  // Operators
  //
  
  Mat2x2T<Scalar>& operator=( const Mat2x2T<Scalar>& _rhs) 
  {
    for (unsigned int i=0; i<2; ++i)
      operator()(i) = _rhs(i);
    return *this;
  }

  /// +=-Operator
  Mat2x2T<Scalar>& operator+=( const Mat2x2T<Scalar>& p )
  {
    elem_[0] += p.elem_[0];
    elem_[1] += p.elem_[1];
    return *this;
  };

  /// -=-Operator
  Mat2x2T<Scalar>& operator-=( const Mat2x2T<Scalar>& p )
  {
    elem_[0] -= p.elem_[0];
    elem_[1] -= p.elem_[1];
    return *this;
  };

  /// /=-Operator
  Mat2x2T<Scalar>& operator/=(const Scalar &s ) 
  {
    elem_[0] /= s;
    elem_[1] /= s;
    return *this;
  };

  /// *=-Operator : Matrix * Scalar
  Mat2x2T<Scalar>& operator*=( Scalar s ) 
  {
    elem_[0] *= s;  elem_[1] *= s;
    return *this;
  };

  /// *=-Operator : Matrix * Matrix
  Mat2x2T<Scalar>& operator*=( const Mat2x2T<Scalar>& p )
  {
    return ( *this = *this * p );
  };

  Mat2x2T<Scalar>& left_mult( const Mat2x2T<Scalar>& p )
  {
    return ( *this = p * *this );
  };

  /// +-Operator 
  Mat2x2T<Scalar>  operator+( const Mat2x2T<Scalar>& p ) const
  {
    return Mat2x2T<Scalar>( elem_[0]+p.elem_[0],
			                      elem_[1]+p.elem_[1] );
  };

  /// --Operator
  Mat2x2T<Scalar>  operator-( const Mat2x2T<Scalar>& p ) const 
  {
    return Mat2x2T<Scalar>( elem_[0]-p.elem_[0],
			    elem_[1]-p.elem_[1] );
  };
 
  ///-Operator 
  Mat2x2T<Scalar>  operator/(const Scalar &s ) const 
  {
    Mat2x2T<Scalar> m( elem_[0]/s, elem_[1]/s);
    return m;
  };
  
  /// *-Operator : Matrix * Scalar
  Mat2x2T<Scalar>  operator*(const Scalar &s ) const
  {
    return Mat2x2T<Scalar>( elem_[0]*s, elem_[1]*s);
  };

  /// friend *-Operator : Scalar * Matrix
  friend Mat2x2T operator*(const Scalar s, Mat2x2T& p )
  {
    return Mat2x2T<Scalar>( p.elem_[0]*s, p.elem_[1]*s );
  };
  
  /// *-Operator : Matrix * Vector
  Vec2 operator*( const Vec2& vec ) const 
  {
    return Vec2(elem_[0]*vec[0] + elem_[1]*vec[1] );
  };

  /// *-Operator : Matrix * Matrix
  Mat2x2T<Scalar>  operator*( const Mat2x2T<Scalar>& p ) const
  {
    Mat2x2T<Scalar> result;
    result.elem_[0] = ( elem_[0]*p.elem_[0][0] +
			elem_[1]*p.elem_[0][1] );
    result.elem_[1] = ( elem_[0]*p.elem_[1][0] +
			elem_[1]*p.elem_[1][1] );
    return result;
  };


  /// read access for matrix elements
  Scalar operator() (int _row, int _col) const { return elem_[_col][_row]; };

  /// write access for matrix elements
  Scalar& operator() (int _row, int _col) { return elem_[_col][_row]; };

  /// read access for matrix elements
  Vec2 operator() (int _col) const { return elem_[_col]; };

  /// write access for matrix elements
  Vec2& operator() (int _col) { return elem_[_col]; };

  bool operator==( const Mat2x2T<Scalar>& p )
  {
    const Scalar epsilon = 0.0001;
    if ( ( fabs (elem_[0][0] - p.elem_[0][0]) > epsilon )  )return false;
    if ( ( fabs (elem_[0][1] - p.elem_[0][1]) > epsilon )  )return false;
    if ( ( fabs (elem_[1][0] - p.elem_[1][0]) > epsilon )  )return false;
    if ( ( fabs (elem_[1][1] - p.elem_[1][1]) > epsilon )  )return false;
    return true;
  };


  //
  // Methods
  //

   /// determinant of 2x2T Matrix
  Scalar  det()
  {
    return ( (elem_[0][0]*elem_[1][1] - elem_[0][1]*elem_[1][0]) );
  };


  /// Transposed Matrix
  Mat2x2T<Scalar> transpose()
  {
    return( Mat2x2T( elem_[0][0], elem_[0][1],
                     elem_[1][0], elem_[1][1] ) );
  }


  ///// Eigenvalues and Eigenvectors
  //int symm_eigenv(Vec2& eigenvals, 
		//  Vec2& eigenvec1, 
		//  Vec2& eigenvec2 );

  
  bool invert()
  {
    Scalar dt = det();

    if (dt == Scalar(0)) return false;
    
    dt = 1.0/ dt;
    Scalar e00 = elem_[0][0];
    elem_[0] = Vec2 ( elem_[1][1] * dt, -elem_[0][1] * dt );
    elem_[1] = Vec2 (-elem_[1][0] * dt,  e00 * dt );

    return true;
  }
  
  // frobenius norm
  Scalar sqr_norm()
  {
    return elem_[0].sqrnorm() + elem_[1].sqrnorm();
  }

  Scalar norm()
  {
    return sqrt (sqr_norm());
  }

  Scalar trace()
  {
    return (elem_[0][0] + elem_[1][1]);
  }

  void identity()
  {
    elem_[0] = Vec2(1,0);
    elem_[1] = Vec2(0,1);
  }

  friend std::ostream& operator<<(std::ostream& os, Mat2x2T<Scalar> const &mat)
  {
    for(int j=0; j<2; ++j)
    {
      for(int i=0; i<2; ++i) 
        os << mat.elem_[i][j] << " ";
      os << " |  ";
    }
    return os;
  }

private:

  Vec2 elem_[2];
};


/// typedef
typedef Mat2x2T< float >  Mat2x2f;
/// typedef
typedef Mat2x2T< double > Mat2x2d;

//=============================================================================
} //namespace Math
//=============================================================================

//=============================================================================
#endif // MAT2X2T_HH defined
