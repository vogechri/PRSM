//=============================================================================
//
//  CLASS Mat3x3T
//
//=============================================================================

#ifndef MAT3X3_HH
#define MAT3X3_HH


//== INCLUDES =================================================================

#include "VectorT.h"
#include <vector>
#include <iostream>
#include <limits>
//== CLASS DEFINITION =========================================================

//TODO : LUDSolve3x3Eqn	has optimal matrix equation solvers -> use them for inversion, etc.

namespace Math {
//-----------------------------------------------------------------------------
/*! \brief defines a 3x3 matrix. Several often used algorithms are implemented as well. 
 *
 * \author C.Vogel
*/
template<class Scalar>
class Mat3x3T
{
public:

  // COLUMN_WISE stored

  typedef VectorT<Scalar,3>    Vec3;
  
  /// Empty constructor - object remains uninitialized.
  Mat3x3T( void ) {};

  /// init with axi meaning column i
  Mat3x3T( const Scalar a00, const Scalar a01,const Scalar a02, 
	   const Scalar a10, const Scalar a11,const Scalar a12, 
	   const Scalar a20, const Scalar a21,const Scalar a22)
  {
    elem_[0][0] = a00; elem_[1][0] = a01; elem_[2][0] = a02;
    elem_[0][1] = a10; elem_[1][1] = a11; elem_[2][1] = a12;
    elem_[0][2] = a20; elem_[1][2] = a21; elem_[2][2] = a22;
  };

  Mat3x3T(const Vec3 a0, const Vec3 a1,const Vec3 a2)
  {
    elem_[0] = a0; elem_[1] = a1; elem_[2] = a2;
  };

  Mat3x3T(const Scalar* _values) 
  {
    elem_[0] = Vec3(&_values[0]);
    elem_[1] = Vec3(&_values[3]);
    elem_[2] = Vec3(&_values[6]);
  }


  /// Copy-Constructor
  Mat3x3T( const Mat3x3T<Scalar> &p )
  {
    elem_[0] = p.elem_[0];
    elem_[1] = p.elem_[1];
    elem_[2] = p.elem_[2];
  };
  

  /// Destructor
  ~Mat3x3T( void ) {};

  
  //
  // Operators
  //
  
  Mat3x3T<Scalar>& operator=( const Mat3x3T<Scalar>& _rhs) 
  {
    for (unsigned int i=0; i<3; ++i)
      operator()(i) = _rhs(i);
    return *this;
  }

  /// +=-Operator
  Mat3x3T<Scalar>& operator+=( const Mat3x3T<Scalar>& p )
  {
    elem_[0] += p.elem_[0];
    elem_[1] += p.elem_[1];
    elem_[2] += p.elem_[2];
    return *this;
  };

  /// -=-Operator
  Mat3x3T<Scalar>& operator-=( const Mat3x3T<Scalar>& p )
  {
    elem_[0] -= p.elem_[0];
    elem_[1] -= p.elem_[1];
    elem_[2] -= p.elem_[2];
    return *this;
  };

  /// /=-Operator
  Mat3x3T<Scalar>& operator/=(const Scalar &s ) 
  {
    elem_[0] /= s;
    elem_[1] /= s;
    elem_[2] /= s;
    return *this;
  };

  /// *=-Operator : Matrix * Scalar
  Mat3x3T<Scalar>& operator*=( Scalar s ) 
  {
    elem_[0] *= s;  elem_[1] *= s;  elem_[2] *= s;
    return *this;
  };

  /// *=-Operator : Matrix * Matrix
  Mat3x3T<Scalar>& operator*=( const Mat3x3T<Scalar>& p )
  {
    return ( *this = *this * p );
  };

  Mat3x3T<Scalar>& left_mult( const Mat3x3T<Scalar>& p )
  {
    return ( *this = p * *this );
  };

  /// +-Operator 
  Mat3x3T<Scalar>  operator+( const Mat3x3T<Scalar>& p ) const
  {
    return Mat3x3T<Scalar>( elem_[0]+p.elem_[0],
			    elem_[1]+p.elem_[1],
			    elem_[2]+p.elem_[2] );
  };

  /// --Operator
  Mat3x3T<Scalar>  operator-( const Mat3x3T<Scalar>& p ) const 
  {
    return Mat3x3T<Scalar>( elem_[0]-p.elem_[0],
			    elem_[1]-p.elem_[1],
			    elem_[2]-p.elem_[2] );
  };
 
  ///-Operator 
  Mat3x3T<Scalar>  operator/(const Scalar &s ) const 
  {
    Mat3x3T<Scalar> m( elem_[0]/s, elem_[1]/s, elem_[2]/s );
    return m;
  };
  
  /// *-Operator : Matrix * Scalar
  Mat3x3T<Scalar>  operator*(const Scalar &s ) const
  {
    return Mat3x3T<Scalar>( elem_[0]*s, elem_[1]*s, elem_[2]*s );
  };

  /// friend *-Operator : Scalar * Matrix
  friend Mat3x3T operator*(const Scalar s, Mat3x3T& p )
  {
    return Mat3x3T<Scalar>( p.elem_[0]*s, p.elem_[1]*s, p.elem_[2]*s );
  };
  
  /// *-Operator : Matrix * Vector
  Vec3 operator*( const Vec3& vec ) const 
  {
    //return (elem_[0]*vec[0] + elem_[1] * vec[1] + elem_[2] * vec[2]);
    return Vec3(elem_[0]*vec[0] + elem_[1] * vec[1] + elem_[2] * vec[2]);
  };

  /// *-Operator : Matrix * Matrix
  Mat3x3T<Scalar>  operator*( const Mat3x3T<Scalar>& p ) const
  {
    Mat3x3T<Scalar> result;
    result.elem_[0] = ( elem_[0]*p.elem_[0][0] +
			elem_[1]*p.elem_[0][1] +
			elem_[2]*p.elem_[0][2] );
    result.elem_[1] = ( elem_[0]*p.elem_[1][0] +
			elem_[1]*p.elem_[1][1] +
			elem_[2]*p.elem_[1][2] );
    result.elem_[2] = ( elem_[0]*p.elem_[2][0] +
			elem_[1]*p.elem_[2][1] +
			elem_[2]*p.elem_[2][2] );
    return result;
  };


  /// read access for matrix elements
  Scalar operator() (int _row, int _col) const { return elem_[_col][_row]; };

  /// write access for matrix elements
  Scalar& operator() (int _row, int _col) { return elem_[_col][_row]; };

  /// read access for matrix elements
  Vec3 operator() (int _col) const { return elem_[_col]; };

  /// write access for matrix elements
  Vec3& operator() (int _col) { return elem_[_col]; };

  bool operator==( const Mat3x3T<Scalar>& p )
  {
    const Scalar epsilon = 0.0001;
    if ( ( fabs (elem_[0][0] - p.elem_[0][0]) > epsilon )  )return false;
    if ( ( fabs (elem_[0][1] - p.elem_[0][1]) > epsilon )  )return false;
    if ( ( fabs (elem_[0][2] - p.elem_[0][2]) > epsilon )  )return false;
    if ( ( fabs (elem_[1][0] - p.elem_[1][0]) > epsilon )  )return false;
    if ( ( fabs (elem_[1][1] - p.elem_[1][1]) > epsilon )  )return false;
    if ( ( fabs (elem_[1][2] - p.elem_[1][2]) > epsilon )  )return false;
    if ( ( fabs (elem_[2][0] - p.elem_[2][0]) > epsilon )  )return false;
    if ( ( fabs (elem_[2][1] - p.elem_[2][1]) > epsilon )  )return false;
    if ( ( fabs (elem_[2][2] - p.elem_[2][2]) > epsilon )  )return false;
    return true;
  };


  //
  // Methods
  //

   /// determinant of 3x3T Matrix
  Scalar  det()
  {
    return ( (elem_[0][1]*elem_[1][2] - elem_[0][2]*elem_[1][1]) * elem_[2][0]
	  +  (elem_[0][2]*elem_[1][0] - elem_[0][0]*elem_[1][2]) * elem_[2][1]
	  +  (elem_[0][0]*elem_[1][1] - elem_[0][1]*elem_[1][0]) * elem_[2][2]);
  };


  /// Transposed Matrix
  Mat3x3T<Scalar> transpose()
  {
    return( Mat3x3T( elem_[0][0], elem_[0][1], elem_[0][2],
		     elem_[1][0], elem_[1][1], elem_[1][2],
		     elem_[2][0], elem_[2][1], elem_[2][2] ) );
  }


  /// Eigenvalues and Eigenvectors
  int symm_eigenv(Vec3& eigenvals, 
		  Vec3& eigenvec1, 
		  Vec3& eigenvec2, 
		  Vec3& eigenvec3);

  
  bool invert()
  {
    Scalar dt = det();

    if ( fabs(dt) < std::numeric_limits<Scalar>::epsilon() ) return false;
    
    dt = 1.0/ dt;

    Scalar i_11 =  under_det(0,0);
    Scalar i_21 = -under_det(0,1);
    Scalar i_31 =  under_det(0,2);

    Scalar i_12 = -under_det(1,0);
    Scalar i_22 =  under_det(1,1);
    Scalar i_32 = -under_det(1,2);

    Scalar i_13 =  under_det(2,0);
    Scalar i_23 = -under_det(2,1);
    Scalar i_33 =  under_det(2,2);

    elem_[0] = Vec3 (i_11 * dt, i_21 * dt, i_31 * dt);
    elem_[1] = Vec3 (i_12 * dt, i_22 * dt, i_32 * dt);
    elem_[2] = Vec3 (i_13 * dt, i_23 * dt, i_33 * dt);

    return true;
  }


  /// assumes a symmetric matrix !
  bool invertSym()
  {
    Scalar dt = det();

    if ( fabs(dt) < std::numeric_limits<Scalar>::epsilon() ) return false;
    
    dt = 1.0 / dt;

    Scalar i_11 =  ( elem_[1][1] * elem_[2][2] - 
                     elem_[1][2] * elem_[2][1] ) * dt;

    Scalar i_21 = - ( elem_[0][1] * elem_[2][2] - 
                      elem_[0][2] * elem_[2][1] ) * dt;

    Scalar i_31 =  ( elem_[0][1] * elem_[1][2] - 
                     elem_[0][2] * elem_[1][1] ) * dt;

    Scalar i_22 =  ( elem_[0][0] * elem_[2][2] - 
                     elem_[0][2] * elem_[2][0] ) * dt;

    Scalar i_32 = -( elem_[0][0] * elem_[1][2] - 
                     elem_[0][2] * elem_[1][0] ) * dt;

    Scalar i_33 =  ( elem_[0][0] * elem_[1][1] - 
                     elem_[0][1] * elem_[1][0] ) * dt;

    elem_[0] = Vec3 (i_11, i_21, i_31);
    elem_[1] = Vec3 (i_21, i_22, i_32);
    elem_[2] = Vec3 (i_31, i_32, i_33);

    return true;
  }
  
  // TODO include the lud decomposition here see ludcmp3x3

  // frobenius norm
  Scalar sqr_norm()
  {
    return elem_[0].sqrnorm() + elem_[1].sqrnorm() + elem_[2].sqrnorm();
  }

  Scalar norm()
  {
    return sqrt (sqr_norm());
  }

  Scalar trace()
  {
    return (elem_[0][0] + elem_[1][1] + elem_[2][2]);
  }

  void identity()
  {
    elem_[0] = Vec3(1,0,0);
    elem_[1] = Vec3(0,1,0);
    elem_[2] = Vec3(0,0,1);
  }
  
  friend std::ostream& operator<<(std::ostream& os, Mat3x3T<Scalar> const &mat)
  {
    for(int j=0; j<3; ++j)
    {
      for(int i=0; i<3; ++i) 
        os << mat.elem_[i][j] << " ";
      os << " |  ";
    }
    return os;
  }


private:

  // row i and column j cleared off matrix -> 2x2 determint computed
  Scalar under_det(int i, int j)
  {
    int i_1 = (i + 1) % 3;
    int i_2 = (i + 2) % 3;
    int j_1 = (j + 1) % 3;
    int j_2 = (j + 2) % 3;

    if (i_1 > i_2) 
    {
      int swap = i_2;
      i_2 = i_1;
      i_1 = swap;
    }
    if (j_1 > j_2) 
    {
      int swap = j_2;
      j_2 = j_1;
      j_1 = swap;
    }

    return 
      elem_[j_1][i_1] * elem_[j_2][i_2] - 
      elem_[j_1][i_2] * elem_[j_2][i_1];
  }

  Vec3 elem_[3];
};


/// typedef
typedef Mat3x3T< float >  Mat3x3f;
/// typedef
typedef Mat3x3T< double > Mat3x3d;

//=============================================================================
} //namespace Math {
//=============================================================================
// include templates in header file (g++ only)
# if !defined(MAT3X3_CC)
#  include "Mat3x3T.cpp"
# endif

//=============================================================================

//=============================================================================
#endif // MAT3X3T_HH defined
