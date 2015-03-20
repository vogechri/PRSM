//=============================================================================
//
//  CLASS Polar
//
//=============================================================================
 
#ifndef POLAR_HH
#define POLAR_HH
 
//== INCLUDES =================================================================
 

#include "../Mat3x3T.h"
#include "My_QuaternionT.h"

//== CLASS DEFINITION =========================================================

namespace Math
{

/*!
 * \class Polar
 * \brief computes the Polar decomposition of a 3x3 matrix.
 *
 * This class computes the Polar decomposition of a 3x3 matrix.
 * The matrix get s split into a scale and rotation part.
 * Such that M = R * S. (R : Rotation, S: Scale)
 *
 * Various output formats like quaternions or matrices are possible.
 *
 * \author C.Vogel
 * 
 */
template<class Scalar>
class Polar
{

public:

  typedef Mat3x3T<Scalar>         Matrix;
  typedef My_QuaternionT<Scalar>  MyQuaternion; 

  enum Interpolation_Mode
    {
      SLERP,
      LERP
    };


  Polar();

  Polar(Scalar w, Scalar x, Scalar y, Scalar z);

  Polar(MyQuaternion _q, Matrix _m);


  ~Polar();

  static void set_interpolation_mode (int _i)
  {
    ip_mode_ = _i;
  }

  //-----------------inits ----------------

  // set the matrix (init for decomposition)
  void set_matrix(Matrix _mat);

  // set the matrix (init for decomposition)
  void set_matrix(Scalar a00,Scalar a01,Scalar a02,
		  Scalar a10,Scalar a11,Scalar a12,
		  Scalar a20,Scalar a21,Scalar a22 );

  // no decomposition needed for give matrix is a rotation matrix
  void set_rot_matrix(Matrix _m);

  // set quaternion directly
  void set_quaternion(MyQuaternion _q);

  // set scale matrix directly
  void set_scale(Matrix _mat);

  // return the complete matrix
  Matrix get_matrix();

  // return complete matrix given by quaternion and scale part
  Matrix get_new_matrix();

  // compute intern matrix as composition of quaternion and scale part
  void update_matrix();

  // set polar to identity
  void identity();

  //------------computations---------------------

  // decomposes intern matrix in rotation-quaternion and scale matrix such that
  // Rot * Scale = _mat
  void decompose();

  // decomposes _mat in rotation-quaternion and scale matrix such that
  // Rot * Scale = _mat
  void decompose( Matrix _mat);

  // decomposes _mat in rotation-matrix and scale matrix such that
  // Rot * Scale = _mat
  void decompose2Matrix( Matrix _mat );

  //  returns composition of two Polar classes
  Polar operator*(const Polar& _p) const ;

  //  returns composition of two Polar classes
  Polar& operator*=(const Polar& _p);


  // return quaternion (rotation part)
  MyQuaternion get_quaternion();

  // interpolate and keep solution stored:
  // different interpolation of rotation and scale part
  void interpolate(const Scalar _interpol , Polar<Scalar> _p);


  // interpolate:
  // different interpolation of rotation and scale part
  // interpolate rotation part by the angle
  // sense is: if one quaternion is the identity then 
  // the new quaternion is the axis and a part of the angle of the
  // other quaternion
  // so a larger area can be interpolated
  void interpolate_angle(const Scalar _interpol , Polar<Scalar> _p);


  // interpolate and return a matrix:
  // different interpolation of rotation and scale part
  Matrix interpolate(Polar<Scalar> _p,
		     const Scalar _interpol );  //const


  // apply whole matrix on vector
  template<class Vec>
  Vec transform(Vec _vec);// const;

  // apply rotation part of matrix on vector
  template<class Vec>
  Vec transform_without_scale(Vec _vec) const;


  // return the rotation angle
  Scalar angle();

  // return the rotation axis
  typename MyQuaternion::Vec3 axis();

  // return just the scale part
  Matrix scale() const;

  // return just the rotation part
  Matrix rotation();

  // return the complete matrix
  Matrix trafo_matrix();

  // extern interpolation of two quaternions by SLERP
  MyQuaternion slerp(const MyQuaternion _q1,
		   const MyQuaternion _q2,
		   const Scalar _interpol) const;

  // extern interpolation of two quaternions by LERP
  MyQuaternion lerp(const MyQuaternion _q1,
		  const MyQuaternion _q2,
		  const Scalar _interpol) const;


private:

  // ! ATTENTION : matrix has to be pure rotation 
  MyQuaternion matrix_to_quaternion(Matrix _mat);


  static const Scalar epsilon;// = 0.00001;
  static int ip_mode_;

  Matrix matrix_;
  Matrix scale_;
  MyQuaternion quaternion_;

};

template<class Scalar>
int Polar< Scalar >::ip_mode_(0);

template<class Scalar>
const Scalar Polar< Scalar >::epsilon (Scalar(0.00001));
} // namespace Math

#if defined(INCLUDE_TEMPLATES) && !defined(Polar_C)
#define Polar_TEMPLATES
#include "Polar.cpp"
#endif

#endif