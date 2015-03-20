//=============================================================================
//
//  CLASS MYQuaternion
//
//=============================================================================
 
#ifndef MY_QUATERNION_HH
#define MY_QUATERNION_HH
 
//== INCLUDES =================================================================
 
#include "Math/VectorT.h"
#include "Math/QuaternionT.h"
#include "Math/Mat3x3T.h"
#include <math.h>
#include <float.h>
//== NAMESPACES  ==============================================================
 
 namespace Math {

//== CLASS DEFINITION =========================================================

template <class Scalar>
class My_QuaternionT : public Math::QuaternionT<Scalar>
{
public:

//#define my_debug

#define W values_[0]
#define X values_[1]
#define Y values_[2]
#define Z values_[3]

#define W_ this->values_[0]
#define X_ this->values_[1]
#define Y_ this->values_[2]
#define Z_ this->values_[3]

  typedef QuaternionT<Scalar>          Base;
  typedef My_QuaternionT<Scalar>       My_Quaternion;
  typedef Mat3x3T<Scalar>              My_Matrix;
  typedef VectorT<Scalar,3>            Vec3;
  typedef VectorT<Scalar,4>            Vec4;
  
  My_QuaternionT(Scalar _w=1.0, Scalar _x=0.0, Scalar _y=0.0, Scalar _z=0.0)
    :Base (_w, _x, _y, _z)  {}
  
  My_QuaternionT(const Vec3& _p)
    :Base (_p)  {}
  
  My_QuaternionT(const Vec4& _v)
    : Base(_v[0], _v[1], _v[2], _v[3]) {}
  
  My_QuaternionT(Vec3 _axis, Scalar _angle)
    : Base(_axis, _angle) {}

  
  // new implementations :
  
  My_Matrix to_rotation_matrix()
  {

    //    this->normalize();

    Scalar 
      //ww(W*W), 
      xx(X_*X_), yy(Y_*Y_), zz(Z_*Z_), wx(W_*X_),
      wy(W_*Y_), wz(W_*Z_), xy(X_*Y_), xz(X_*Z_), yz(Y_*Z_);
       
    My_Matrix m;

    m(0,0) = 1.0 - 2.0*(yy + zz);//ww + xx - yy - zz;
    m(1,0) = 2.0*(xy - wz);
    m(2,0) = 2.0*(xz + wy);
 
    m(0,1) = 2.0*(xy + wz);
    m(1,1) = 1.0 - 2.0*(zz + xx);//ww - xx + yy - zz;
    m(2,1) = 2.0*(yz - wx);
 
    m(0,2) = 2.0*(xz - wy);
    m(1,2) = 2.0*(yz + wx);
    m(2,2) = 1.0 - 2.0*(yy + xx);//ww - xx - yy + zz;

    return m;
  }

  My_Quaternion operator*(const Scalar& _s) 
  { return My_Quaternion(W_ *_s , X_ *_s , Y_ *_s , Z_ *_s ); }

  My_Quaternion& operator*=(const Scalar& _s)
  { return *this = *this * _s; }


  My_Quaternion operator*(const My_Quaternion& _q) const 
  { return My_Quaternion(W_*_q.W - X_*_q.X - Y_*_q.Y - Z_*_q.Z,
			 W_*_q.X + X_*_q.W + Y_*_q.Z - Z_*_q.Y,
			 W_*_q.Y - X_*_q.Z + Y_*_q.W + Z_*_q.X,
			 W_*_q.Z + X_*_q.Y - Y_*_q.X + Z_*_q.W); }


  /// quaternion *= quaternion
  My_Quaternion& operator*=(const My_Quaternion& _q)
  { return *this = *this * _q; }


  Scalar angle()
  {
    //    this->normalize();
    return Scalar(2) * acos(W_);
  }
  
  Scalar grad_angle()
  {
    //    this->normalize();
    return Scalar(360) * acos(W_) / M_PI;
  }

  Vec3 axis() const
  {
    //if (X==Y==Z == Scalar(0)) -> axis = 0,0,0
    Vec3 axis(X_, Y_, Z_);
    if ( axis.sqrnorm() > 0)
      axis.normalize();
    return axis;
  }

  Scalar dot(My_Quaternion _q)
  {
    return W_ * _q.W + X_ * _q.X + Y_ * _q.Y + Z_ * _q.Z ;
  }

  bool null_axis()
  {
    return ((fabs (X_) < 0.001) && (fabs (Y_) < 0.001) && (fabs (Z_) < 0.001));
  }

  My_Quaternion slerp(My_Quaternion _q , Scalar _interpol)
  {

#ifdef my_debug
    if (_interpol > Scalar(1.0001)) 
      std::cerr << "error in slerp: " << _interpol << " not in [0,1]\n";
    if (_interpol < Scalar(-0.0001))
      std::cerr << "error in slerp: " << _interpol << " not in [0,1]\n";
#endif    

    // winkel immer kleiner 180° -> nicht von orientierung der achsen abhängig
    // winkel zwischen Achsen


    if (dot(_q) < Scalar(0)) 
 	_q = _q * Scalar(-1);

    Scalar alpha;
    Scalar beta;

    Scalar skp ( axis() | _q.axis() );
    skp = std::max (skp, -skp);
    Scalar theta ( acos( skp ) );

    Scalar sin_theta( sin (theta) );

    if (fabs(sin_theta) >  0.05)
    {
      Scalar denom = Scalar(1) / sin_theta;
      alpha = sin( (Scalar (1) - _interpol) * theta) * denom;
      beta  = sin( _interpol * theta) * denom;
    }
    else // bei 0 ist sin ~ identitaet , schutz vor 0/0:
    {
      alpha = (Scalar (1) - _interpol);
      beta  = _interpol;
    }

    // normalize because of interpolation of several handles
    return (*this * alpha + _q * beta).normalize();
  }


  My_Quaternion lerp(My_Quaternion _q , Scalar _interpol)
  {

#ifdef my_debug
    if (_interpol > Scalar(1.0001)) 
      std::cerr << "error in lerp: " << _interpol << " not in [0,1]\n";
    if (_interpol < Scalar(-0.0001)) 
      std::cerr << "error in lerp: " << _interpol << " not in [0,1]\n";
#endif

    Scalar alpha;
    Scalar beta;

    if (dot(_q) < Scalar(0)) 
      _q = _q * Scalar(-1);

    alpha = (Scalar (1) - _interpol);
    beta  = _interpol;

    return (*this * alpha + _q * beta).normalize();
  }


  My_Quaternion interpolate_angle(My_Quaternion _q , Scalar _interpol)
  {

#ifdef my_debug
    if (_interpol > Scalar(1.0001)) 
      std::cerr << "error in interpolate_angle: " 
		<< _interpol << " not in [0,1]\n";
    if (_interpol < Scalar(-0.0001)) 
      std::cerr << "error in interpolate_angle: " 
		<< _interpol << " not in [0,1]\n";
#endif

    Scalar alpha = (Scalar (1) - _interpol);
    Scalar beta  = _interpol;
    
    if (null_axis()) 
      {
	if (_q.null_axis()) return  My_Quaternion(); // id
	return  My_Quaternion( _q.axis(), _q.angle() * beta);
      }
      if (_q.null_axis()) return  My_Quaternion(axis (),  alpha  * angle() );

#ifdef my_debug
      std::cerr << "error: interpolate_angle: " 
		<< "no quaternion is identity\n";
#endif

    return slerp(  _q , _interpol);

  }

//#undef my_debug

#undef W
#undef X
#undef Y
#undef Z

#undef W_
#undef X_
#undef Y_
#undef Z_

};

/// typedef
typedef My_QuaternionT<float>  MyQuaternionf;
/// typedef
typedef My_QuaternionT<double> MyQuaterniond;
//=============================================================================
} // namespace Math
//=============================================================================

#endif
