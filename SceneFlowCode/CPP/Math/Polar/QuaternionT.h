//=============================================================================
//
//  CLASS Quaternion
//
//=============================================================================
 
#ifndef _QUATERNION_HH
#define _QUATERNION_HH
 
//== INCLUDES =================================================================
 
#include "../VectorT.h"
#include "../Matrix4x4T.h"
 
 
//== NAMESPACES  ==============================================================
 
 namespace Math {
 
 
//== CLASS DEFINITION =========================================================
 
 
   template <class Scalar>
   class QuaternionT : public VectorT<Scalar,4>
   {
   public:
     
#define W values_[0]
#define X values_[1]
#define Y values_[2]
#define Z values_[3]
 
     
     typedef QuaternionT<Scalar>  Quaternion;
     typedef VectorT<Scalar,3>    Vec3;
     typedef VectorT<Scalar,4>    Vec4;
     typedef Matrix4x4T<Scalar>   Matrix;
     
 
     QuaternionT(Scalar _w=1.0, Scalar _x=0.0, Scalar _y=0.0, Scalar _z=0.0)
       : Vec4(_w, _x, _y, _z) {}
  
     QuaternionT(const Vec3& _p)
       : Vec4(0, _p[0], _p[1], _p[2]) {}
 
     QuaternionT(const Vec4& _v)
       : Vec4(_p[0], _p[1], _p[2], _p[3]) {}
   
     QuaternionT(Vec3 _axis, Scalar _angle)
     {
       _axis.normalize();
       Scalar theta = 0.5 * _angle;
       Scalar sin_theta = sin(theta);
       W = cos(theta);
       X = sin_theta * _axis[0];
       Y = sin_theta * _axis[1];
       Z = sin_theta * _axis[2];
     }
  
     
     void identity() { W=1.0; X=Y=Z=0.0; }


     Quaternion conjugate() const
     { return Quaternion(W, -X, -Y, -Z); }

     
     Quaternion invert() const
     { return conjugate()/sqrnorm(); }

  
     Quaternion operator*(const Quaternion& _q) const 
     { return Quaternion(W*_q.W - X*_q.X - Y*_q.Y - Z*_q.Z,
			 W*_q.X + X*_q.W + Y*_q.Z - Z*_q.Y,
			 W*_q.Y - X*_q.Z + Y*_q.W + Z*_q.X,
			 W*_q.Z + X*_q.Y - Y*_q.X + Z*_q.W); }
     
 
     Quaternion& operator*=(const Quaternion& _q)
     { return *this = *this * _q; }
   

     Vec3 rotate(const Vec3& _v)
     { return *this * _v * conjugate(); }
 
     
     void axis_angle(Vec3& _axis, Scalar& _angle) const
     {
       if (fabs(W) > 0.999999)
       {
	 _axis  = Vec3(1,0,0);
	 _angle = 0;
       }
       else
       {
	 _angle = 2.0 * acos(W);
	 _axis  = Vec3(X, Y, Z).normalize();
       }
     }
 
  
 
     Matrix rotation_matrix() const
     {
       Scalar 
	 ww(W*W), xx(X*X), yy(Y*Y), zz(Z*Z), wx(W*X),
	 wy(W*Y), wz(W*Z), xy(X*Y), xz(X*Z), yz(Y*Z);
       
       Matrix m;
       
       m(0,0) = ww + xx - yy - zz;
       m(1,0) = 2.0*(xy + wz);
       m(2,0) = 2.0*(xz - wy);
 
       m(0,1) = 2.0*(xy - wz);
       m(1,1) = ww - xx + yy - zz;
       m(2,1) = 2.0*(yz + wx);
       m(0,2) = 2.0*(xz + wy);
       m(1,2) = 2.0*(yz - wx);
       m(2,2) = ww - xx - yy + zz;
       
       m(0,3) = m(1,3) = m(2,3) = m(3,0) = m(3,1) = m(3,2) = 0.0;
       m(3,3) = 1.0;
     
       return m;
  }

  
     Matrix right_mult_matrix() const 
     {
       Matrix m;
       m(0,0) =  W; m(0,1) = -X; m(0,2) = -Y; m(0,3) = -Z;
       m(1,0) =  X; m(1,1) =  W; m(1,2) = -Z; m(1,3) =  Y;
       m(2,0) =  Y; m(2,1) =  Z; m(2,2) =  W; m(2,3) = -X;
       m(3,0) =  Z; m(3,1) = -Y; m(3,2) =  X; m(3,3) =  W;
       return m;
     }


     Matrix left_mult_matrix() const 
     {
       Matrix m;
       m(0,0) =  W; m(0,1) = -X; m(0,2) = -Y; m(0,3) = -Z;
       m(1,0) =  X; m(1,1) =  W; m(1,2) =  Z; m(1,3) = -Y;
       m(2,0) =  Y; m(2,1) = -Z; m(2,2) =  W; m(2,3) =  X;
       m(3,0) =  Z; m(3,1) =  Y; m(3,2) = -X; m(3,3) =  W;
       return m;
     }


#undef W
#undef X
#undef Y
#undef Z
};
 
 
   typedef QuaternionT<float>  Quaternionf;
   typedef QuaternionT<double> Quaterniond;
 
  
//=============================================================================
} // namespace ACG
//=============================================================================
#endif // _QUATERNION_HH defined
//=============================================================================
