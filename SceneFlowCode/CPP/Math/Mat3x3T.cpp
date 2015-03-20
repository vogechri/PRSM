//=============================================================================
//
//  CLASS Mat3x3 - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//-***************************************************************************
// for g++ (template defs only compiled through header file)

#if !defined(MAT3X3_CC)
#define MAT3X3_CC
//-***************************************************************************

//== INCLUDES =================================================================

#include "Mat3x3T.h"
#include <math.h>
#include "eigenval.h"

#include <iostream>

//-----------------------------------------------------------------------------
namespace Math {

//== IMPLEMENTATION ========================================================== 


#define MAX_ITER 100
#define EPS 0.000001

template<class Scalar>
int
Mat3x3T<Scalar>::
symm_eigenv(Vec3& eigenvals,
	    Vec3& eigenvec1, Vec3& eigenvec2, Vec3& eigenvec3)
{
  unsigned int i, j;
  int num_iter=0;
  Scalar theta, t, c, s;
  Mat3x3T<Scalar>
    V(1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0);
  Mat3x3T<Scalar> R;
  Mat3x3T<Scalar> A=*this;
  
  while (num_iter < MAX_ITER) {
    
    // find largest off-diagonal m
    if ( fabs(A.elem_[0][1]) < fabs(A.elem_[0][2]))
    {
      if ( fabs(A.elem_[0][2]) < fabs(A.elem_[1][2]) )
        i = 1, j = 2;
      else
        i = 0, j = 2;
    }
    else 
    {
      if ( fabs(A.elem_[0][1]) < fabs(A.elem_[1][2]) ) 
        i = 1, j = 2;
      else
        i = 0, j = 1;
    }
    
    if ( fabs(A.elem_[i][j]) < EPS ) break;
    

    // compute Jacobi-Rotation
    theta = Scalar(0.5) * (A.elem_[j][j] - A.elem_[i][i]) / A.elem_[i][j];
    t = Scalar(1.0) / (fabs(theta) + sqrt( Scalar(1.0) + theta*theta ));    //  orig: t = 1.0 / (fabs(theta) + sqrt( 1.0 + theta*theta )); but it didnt work 
    if (theta > 0) t *= -1;

    c = Scalar(1.0) / sqrt(Scalar(1.0) + t*t);
    s = t*c;

    R = Mat3x3T<Scalar>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    R.elem_[i][i] = R.elem_[j][j] = c;
    R.elem_[i][j] = s;
    R.elem_[j][i] = -s;

    A = (R.transpose() * A) * R;
    V *= R;
    
    num_iter++;
  }


  if (num_iter < MAX_ITER) {

    // sort and return
    int sorted[3];
    Scalar d[3]={A.elem_[0][0], A.elem_[1][1], A.elem_[2][2]};

    if (d[0] > d[1])
    {
      if (d[1] > d[2])
      {
	sorted[0] = 0, sorted[1] = 1, sorted[2] = 2;
      }
      else // d2 > d1
      {
	if (d[0] > d[2])
	{
	  sorted[0] = 0, sorted[1] = 2, sorted[2] = 1;
	}
	else // d2 > d0
	{
	  sorted[0] = 2, sorted[1] = 0, sorted[2] = 1;
	}
      }
    }
    else // d1 > d0
    {
      if (d[0] > d[2])
      {
	sorted[0] = 1, sorted[1] = 0, sorted[2] = 2;
      }
      else // d2 > d0
      {
	if (d[1] > d[2]) 
	{
	  sorted[0] = 1, sorted[1] = 2, sorted[2] = 0;
	}
	else // d2 > d1
	{
	  sorted[0] = 2, sorted[1] = 1, sorted[2] = 0;
	}
      }
    }
                	
    eigenvals = Vec3( d[sorted[0]],
		      d[sorted[1]],
		      d[sorted[2]]);
		      
    eigenvec1 = V.elem_[sorted[0]];
    eigenvec2 = V.elem_[sorted[1]];

    eigenvec1.normalize();
    eigenvec2.normalize();

    eigenvec3 = eigenvec1 % eigenvec2;
    eigenvec3.normalize();
    
    return(1);
  }
  else return(0);
}

#undef MAX_ITER
#undef EPS

//=============================================================================

} //end namespace Math {

//=============================================================================
#endif
