//#if !defined(Polar_TEMPLATES)
#if !defined(Polar_C)
#define Polar_C

#include "Polar.h"

//-----------------------------------------------------------------------------
namespace Math
{

//== IMPLEMENTATION ========================================================== 


template<class Scalar>
Polar<Scalar>::
Polar() : quaternion_()
{
  scale_ =  Matrix(1 ,0 ,0,
		   0, 1, 0, 
		   0, 0, 1);
}


template<class Scalar>
Polar<Scalar>::
Polar(Scalar w, Scalar x, Scalar y, Scalar z) : quaternion_(w,x,y,z)
{
  scale_ =  Matrix(1 ,0 ,0,
		   0, 1, 0, 
		   0, 0, 1);
}


template<class Scalar>
Polar<Scalar>::
Polar( MyQuaternion _q, Matrix _m )
{
  scale_      = _m;
  quaternion_ = _q;
  matrix_ = get_new_matrix();
}


template<class Scalar>
Polar<Scalar>::
~Polar() {} ;


//-----------------inits ----------------

template<class Scalar>
void
Polar<Scalar>::
set_matrix( Matrix _mat )
{
  matrix_ =  _mat;
}


template<class Scalar>
void
Polar<Scalar>::
set_matrix( Scalar a00,Scalar a01,Scalar a02,
	    Scalar a10,Scalar a11,Scalar a12,
	    Scalar a20,Scalar a21,Scalar a22 )
{
  matrix_ =  Matrix(a00 ,a01 ,a02,
		    a10, a11, a12, 
		    a20, a21, a22);
}


template<class Scalar>
void
Polar<Scalar>::
set_quaternion( MyQuaternion _q )
{
  quaternion_ = _q;
}


template<class Scalar>
void
Polar<Scalar>::
set_scale( Matrix _mat )
{
  scale_ = _mat;
}


template<class Scalar>
typename Polar<Scalar>::Matrix
Polar<Scalar>::
get_matrix()
{
  return matrix_;
}


template<class Scalar>
typename Polar<Scalar>::Matrix 
Polar<Scalar>::
get_new_matrix()
{
  return quaternion_.to_rotation_matrix() * scale_;
}


template<class Scalar>
void
Polar<Scalar>::
update_matrix()
{
  matrix_ = quaternion_.to_rotation_matrix() * scale_;
}


template<class Scalar>
void
Polar<Scalar>::
identity()
{
  quaternion_.identity();
  scale_.identity();
  matrix_.identity();
}

//------------computations---------------------

template<class Scalar>
Polar<Scalar>
Polar<Scalar>::
operator*(const Polar<Scalar>& _p) const 
{ 
  return Polar( this->quaternion_ * _p.quaternion_, 
		this->scale_ * _p.scale_         );
}


template<class Scalar>
Polar<Scalar>&
Polar<Scalar>::
operator*=(const Polar<Scalar>& _p)
{ return *this = *this * _p; }


// decomposes _mat in rotation-quaternion and scale matrix such that
// Rot * Scale = _mat
template<class Scalar>
void
Polar<Scalar>::
decompose()
{
  decompose( matrix_);
}


template<class Scalar>
void
Polar<Scalar>::
set_rot_matrix( Matrix _m )
{
  quaternion_ = matrix_to_quaternion(_m);
  scale_.identity();
};


template<class Scalar>
typename Polar<Scalar>::MyQuaternion 
Polar<Scalar>::
get_quaternion()
{
  return quaternion_;
}


template<class Scalar>
void
Polar<Scalar>::
interpolate(const Scalar _interpol , Polar<Scalar> _p)
{
  switch (Interpolation_Mode(ip_mode_))
    {
    case SLERP:
      quaternion_ = quaternion_.slerp(_p.quaternion_, _interpol);
      break;
    case LERP:
      quaternion_ = quaternion_.lerp(_p.quaternion_, _interpol);
      break;
    default:
      quaternion_ = quaternion_.slerp(_p.quaternion_, _interpol);
    }
  scale_ = (scale_ * (Scalar (1) - _interpol)) + (_p.scale_ * _interpol);
}


template<class Scalar>
void
Polar<Scalar>::
interpolate_angle(const Scalar _interpol , Polar<Scalar> _p)
{
  //    quaternion_ = quaternion_.interpolate_angle(_p.quaternion_, _interpol);

  switch (Interpolation_Mode(ip_mode_))
    {
    case SLERP:
      quaternion_ = quaternion_.interpolate_angle(_p.quaternion_, _interpol);
      break;
    case LERP:
      quaternion_ = quaternion_.lerp(_p.quaternion_, _interpol);
      break;
    default:
      quaternion_ = quaternion_.interpolate_angle(_p.quaternion_, _interpol);
    }
  scale_ = (scale_ * (Scalar (1) - _interpol)) + (_p.scale_ * _interpol);
}


template<class Scalar>
typename Polar<Scalar>::Matrix
Polar<Scalar>::  
interpolate(Polar<Scalar> _p,
	    const Scalar _interpol )  //const
{
  Matrix scale;
  MyQuaternion quat;

  //    quat  = quaternion_.slerp(_p.quaternion_, _interpol);
  switch (Interpolation_Mode(ip_mode_))
    {
    case SLERP:
      quat = quaternion_.slerp(_p.quaternion_, _interpol);
      break;
    case LERP:
      quat = quaternion_.lerp(_p.quaternion_, _interpol);
      break;
    default:
      quat = quaternion_.slerp(_p.quaternion_, _interpol);
    }

  scale = (scale_ * (Scalar (1) - _interpol)) + (_p.scale_ * _interpol);

  return quat.to_rotation_matrix() * scale;
}


//template<class Vec>
template<class Scalar>
template<class Vec>
Vec
Polar<Scalar>::
transform(Vec _vec) //const
{
  Vec vec (scale_ * _vec);
  // normalize ??
  return quaternion_.rotate( vec );
}


//template<class Vec>
template<class Scalar>
template<class Vec>
Vec
Polar<Scalar>::
transform_without_scale(Vec _vec) const
{
  return  quaternion_.rotate( _vec );
}


template<class Scalar>
//typename Polar<Scalar>::Scalar
Scalar
Polar<Scalar>::  
angle()
{
  return quaternion_.angle();
}
  

template<class Scalar>
typename My_QuaternionT<Scalar>::Vec3
Polar<Scalar>::  
axis()
{
  return quaternion_.axis();
}


template<class Scalar>
typename Polar<Scalar>::Matrix
Polar<Scalar>::  
scale() const
{
  return scale_;
}


template<class Scalar>
typename Polar<Scalar>::Matrix
Polar<Scalar>::  
rotation()
{
  return quaternion_.to_rotation_matrix();
}


template<class Scalar>
typename Polar<Scalar>::Matrix
Polar<Scalar>::  
trafo_matrix()
{
  return rotation() * scale();
}


template<class Scalar>
typename Polar<Scalar>::MyQuaternion 
Polar<Scalar>::
slerp( const MyQuaternion _q1, const MyQuaternion _q2, const Scalar _interpol) const
{
  MyQuaternion q = _q1;
  return (q.slerp(_q2,_interpol));
}


template<class Scalar>
typename Polar<Scalar>::MyQuaternion 
Polar<Scalar>::
lerp( const MyQuaternion _q1, const MyQuaternion _q2, const Scalar _interpol) const
{
  MyQuaternion q = _q1;
  return (q.lerp(_q2,_interpol));
}


template<class Scalar>
typename Polar<Scalar>::MyQuaternion 
Polar<Scalar>::
matrix_to_quaternion(Matrix _mat)
{
  Scalar w,x,y,z;
  // now heading for quaternion as in : 
  // Shoemake Animating rotation  with quaternion curves
  Scalar squared_w = 0.25 * (_mat.trace() + 1.0);

  if (squared_w  > 0.0001)
    {
      // immer positiv -> winkel zwischen -180 und 180 grad
      // da cos (phi/2) == 0 bei phi = +- 180 grad

      w = sqrt(squared_w);
      
      Scalar denom = 1.0 / (4.0 * w);

      x = (_mat(1,2) - _mat(2,1)) * denom;
      y = (_mat(2,0) - _mat(0,2)) * denom;
      z = (_mat(0,1) - _mat(1,0)) * denom;
      
    }
  else // w == 0 -> angle = 0 
    {
      w = 0.0;
      Scalar squared_x = -0.5 * (_mat(1,1) + _mat(2,2));
      if (squared_x > 0.001)
	{
	  x = sqrt(squared_x);
	  y = _mat(0,1) / (2.0 * x);
	  z = _mat(0,2) / (2.0 * x);
	}
      else // w == x == 0
	{
	  x = 0.0;
	  Scalar squared_y = 0.5 * ( 1.0 - _mat(2,2));
	  
	  if (squared_y > 0.001)
	    {
	      y = sqrt(squared_y);
	      z = _mat(1,2) / (2.0 * y);
	    }
	  else  // w==x==y==0
	    {
	      y = 0.0;
	      z = 1.0;
	    }

	}
    }
  return MyQuaternion (w,x,y,z);
}


template<class Scalar>
void
Polar<Scalar>::
decompose( Matrix _mat )
{
  Matrix oldmat;

  Matrix newmat = _mat;
    
  do
    {
      oldmat = newmat;

      Matrix help (newmat);
      help.invert();

      newmat = (oldmat + help.transpose()) * 0.5;
      
    } while ((oldmat-newmat).sqr_norm() > epsilon);


  // if not pure rotation -> flip (180° = flip!)
  if (newmat.det() < 0.0)
    newmat *= -1.0;

  {
    Matrix help (newmat);
    help.invert();
    scale_ = help * _mat;
  }

  quaternion_ = matrix_to_quaternion( newmat );
}

template<class Scalar>
void
Polar<Scalar>::
decompose2Matrix( Matrix _mat )
{
  Matrix oldmat;

  Matrix newmat = _mat;
    
  do
    {
      oldmat = newmat;

      Matrix help (newmat);
      help.invert();

      newmat = (oldmat + help.transpose()) * 0.5;
      
    } while ((oldmat-newmat).sqr_norm() > epsilon);


  // if not pure rotation -> flip (180° = flip!)
  if (newmat.det() < 0.0)
    newmat *= -1.0;

  {
    Matrix help (newmat);
    help.invert();
    scale_ = help * _mat;
  }

  matrix_ = newmat;
}


} // namespace Math
#endif // defined