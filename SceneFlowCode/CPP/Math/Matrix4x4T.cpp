//=============================================================================
//
//  CLASS Matrix4x4T - IMPLEMENTATION
//
//=============================================================================

#if !defined(_MATRIX4X4_C)
#define _MATRIX4X4_C


//== INCLUDES =================================================================


#include "Matrix4x4T.h"
//#include <ACG/Utils/NumLimitsT.hh>


//== IMPLEMENTATION ========================================================== 


namespace Math {


#define MAT(m,r,c)  ((m)[(r)+((c)<<2)])
#define M(r,w)      (MAT(mat_,r,w))


//-----------------------------------------------------------------------------


template <typename Scalar> 
Matrix4x4T<Scalar>
Matrix4x4T<Scalar>::
operator* (const Matrix4x4T<Scalar>& _rhs) const
{
#define RHS(row,col)  MAT(_rhs.mat_, row,col)
#define TMP(row,col)  MAT(tmp.mat_, row,col)

  Matrix4x4T<Scalar> tmp;
  Scalar mi0, mi1, mi2, mi3;
  int i;

  for (i = 0; i < 4; i++) {
    mi0=M(i,0);  mi1=M(i,1);  mi2=M(i,2);  mi3=M(i,3);
    TMP(i,0) = mi0*RHS(0,0) + mi1*RHS(1,0) + mi2*RHS(2,0) + mi3*RHS(3,0);
    TMP(i,1) = mi0*RHS(0,1) + mi1*RHS(1,1) + mi2*RHS(2,1) + mi3*RHS(3,1);
    TMP(i,2) = mi0*RHS(0,2) + mi1*RHS(1,2) + mi2*RHS(2,2) + mi3*RHS(3,2);
    TMP(i,3) = mi0*RHS(0,3) + mi1*RHS(1,3) + mi2*RHS(2,3) + mi3*RHS(3,3);
  }

  return tmp;

#undef RHS
#undef TMP
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
Matrix4x4T<Scalar>&
Matrix4x4T<Scalar>::
operator*= (const Matrix4x4T<Scalar>& _rhs)
{
#define RHS(row,col)  MAT(_rhs.mat_, row,col)

  int i;
  Scalar mi0, mi1, mi2, mi3;

  for (i = 0; i < 4; i++) 
  {
    mi0=M(i,0);  mi1=M(i,1);  mi2=M(i,2);  mi3=M(i,3);
    M(i,0) = mi0 * RHS(0,0) + mi1 * RHS(1,0) + mi2 * RHS(2,0) + mi3 * RHS(3,0);
    M(i,1) = mi0 * RHS(0,1) + mi1 * RHS(1,1) + mi2 * RHS(2,1) + mi3 * RHS(3,1);
    M(i,2) = mi0 * RHS(0,2) + mi1 * RHS(1,2) + mi2 * RHS(2,2) + mi3 * RHS(3,2);
    M(i,3) = mi0 * RHS(0,3) + mi1 * RHS(1,3) + mi2 * RHS(2,3) + mi3 * RHS(3,3);
  }

  return *this;

#undef RHS
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
Matrix4x4T<Scalar>&
Matrix4x4T<Scalar>::
leftMult(const Matrix4x4T<Scalar>& _rhs)
{
#define RHS(row,col)  MAT(_rhs.mat_, row,col)
  int i;
  Scalar m0i, m1i, m2i, m3i; 
  for(i=0;i<4;i++)
  {
    m0i = M(0,i);  m1i = M(1,i);  m2i = M(2,i);  m3i = M(3,i); 
    M(0,i) = RHS(0,0)*m0i + RHS(0,1)*m1i + RHS(0,2)*m2i + RHS(0,3)*m3i;
    M(1,i) = RHS(1,0)*m0i + RHS(1,1)*m1i + RHS(1,2)*m2i + RHS(1,3)*m3i;
    M(2,i) = RHS(2,0)*m0i + RHS(2,1)*m1i + RHS(2,2)*m2i + RHS(2,3)*m3i;
    M(3,i) = RHS(3,0)*m0i + RHS(3,1)*m1i + RHS(3,2)*m2i + RHS(3,3)*m3i;
  }
  return *this;
#undef RHS
}  


//-----------------------------------------------------------------------------


template <typename Scalar> 
template <typename T>
VectorT<T,4>
Matrix4x4T<Scalar>::
operator*(const VectorT<T,4>& _v) const
{
  return VectorT<T,4> (
    M(0,0)*_v[0] + M(0,1)*_v[1] + M(0,2)*_v[2] + M(0,3)*_v[3],
    M(1,0)*_v[0] + M(1,1)*_v[1] + M(1,2)*_v[2] + M(1,3)*_v[3],
    M(2,0)*_v[0] + M(2,1)*_v[1] + M(2,2)*_v[2] + M(2,3)*_v[3],
    M(3,0)*_v[0] + M(3,1)*_v[1] + M(3,2)*_v[2] + M(3,3)*_v[3]);
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
template <typename T>
VectorT<T,3>
Matrix4x4T<Scalar>::
transform_point(const VectorT<T,3>& _v) const
{
  Scalar x = M(0,0)*_v[0] + M(0,1)*_v[1] + M(0,2)*_v[2] + M(0,3);
  Scalar y = M(1,0)*_v[0] + M(1,1)*_v[1] + M(1,2)*_v[2] + M(1,3);
  Scalar z = M(2,0)*_v[0] + M(2,1)*_v[1] + M(2,2)*_v[2] + M(2,3);
  Scalar w = M(3,0)*_v[0] + M(3,1)*_v[1] + M(3,2)*_v[2] + M(3,3);

  if (w)
  {
    w = 1.0 / w;
    return VectorT<T,3>(x*w, y*w, z*w);
  } 
  else return VectorT<T,3>(x, y, z);
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
template <typename T>
VectorT<T,3>
Matrix4x4T<Scalar>::
transform_vector(const VectorT<T,3>& _v) const
{
  Scalar x = M(0,0)*_v[0] + M(0,1)*_v[1] + M(0,2)*_v[2];
  Scalar y = M(1,0)*_v[0] + M(1,1)*_v[1] + M(1,2)*_v[2];
  Scalar z = M(2,0)*_v[0] + M(2,1)*_v[1] + M(2,2)*_v[2];
  return VectorT<T,3>(x, y, z);
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
void
Matrix4x4T<Scalar>::
clear()
{
  Scalar* m = mat_;
  *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; *m   = 0.0; 
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
void
Matrix4x4T<Scalar>::
identity()
{
  Scalar* m = mat_;
  *m++ = 1.0; *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 1.0; *m++ = 0.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 0.0; *m++ = 1.0; *m++ = 0.0; 
  *m++ = 0.0; *m++ = 0.0; *m++ = 0.0; *m   = 1.0; 
}


//-----------------------------------------------------------------------------


template <typename Scalar> 
void
Matrix4x4T<Scalar>::
transpose()
{
  Scalar tmp;
  for( int i=0; i<4; i++ ) 
  {
    for( int j=i+1; j<4; j++ ) 
    {
      tmp           = MAT(mat_,i,j);
      MAT(mat_,i,j) = MAT(mat_,j,i);
      MAT(mat_,j,i) = tmp;
    }
  }
}
  

//-----------------------------------------------------------------------------


/*
 * Compute inverse of 4x4 transformation matrix.
 * Taken from Mesa3.1
 * Code contributed by Jacques Leroy jle@star.be */
template <typename Scalar> 
bool
Matrix4x4T<Scalar>::
invert()
{
#define SWAP_ROWS(a, b) { Scalar *_tmp = a; (a)=(b); (b)=_tmp; }

  Scalar wtmp[4][8];
  Scalar m0, m1, m2, m3, s;
  Scalar *r0, *r1, *r2, *r3;
  
  r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];
  
  r0[0] = M(0,0); r0[1] = M(0,1);
  r0[2] = M(0,2); r0[3] = M(0,3);
  r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0;
    
  r1[0] = M(1,0); r1[1] = M(1,1);
  r1[2] = M(1,2); r1[3] = M(1,3);
  r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0;
  
  r2[0] = M(2,0); r2[1] = M(2,1);
  r2[2] = M(2,2); r2[3] = M(2,3);
  r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0;
  
  r3[0] = M(3,0); r3[1] = M(3,1);
  r3[2] = M(3,2); r3[3] = M(3,3);
  r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;
  

  /* choose pivot - or die */
  if (fabs(r3[0])>fabs(r2[0])) SWAP_ROWS(r3, r2);
  if (fabs(r2[0])>fabs(r1[0])) SWAP_ROWS(r2, r1);
  if (fabs(r1[0])>fabs(r0[0])) SWAP_ROWS(r1, r0);
  if (0.0 == r0[0])  return false;

  
  /* eliminate first variable     */
  m1 = r1[0]/r0[0]; m2 = r2[0]/r0[0]; m3 = r3[0]/r0[0];
  s = r0[1]; r1[1] -= m1 * s; r2[1] -= m2 * s; r3[1] -= m3 * s;
  s = r0[2]; r1[2] -= m1 * s; r2[2] -= m2 * s; r3[2] -= m3 * s;
  s = r0[3]; r1[3] -= m1 * s; r2[3] -= m2 * s; r3[3] -= m3 * s;
  s = r0[4];
  if (s != 0.0) { r1[4] -= m1 * s; r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r0[5];
  if (s != 0.0) { r1[5] -= m1 * s; r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r0[6];
  if (s != 0.0) { r1[6] -= m1 * s; r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r0[7];
  if (s != 0.0) { r1[7] -= m1 * s; r2[7] -= m2 * s; r3[7] -= m3 * s; }
  

  /* choose pivot - or die */
  if (fabs(r3[1])>fabs(r2[1])) SWAP_ROWS(r3, r2);
  if (fabs(r2[1])>fabs(r1[1])) SWAP_ROWS(r2, r1);
  if (0.0 == r1[1])  return false;
  

  /* eliminate second variable */
  m2 = r2[1]/r1[1]; m3 = r3[1]/r1[1];
  r2[2] -= m2 * r1[2]; r3[2] -= m3 * r1[2];
  r2[3] -= m2 * r1[3]; r3[3] -= m3 * r1[3];
  s = r1[4]; if (0.0 != s) { r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r1[5]; if (0.0 != s) { r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r1[6]; if (0.0 != s) { r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r1[7]; if (0.0 != s) { r2[7] -= m2 * s; r3[7] -= m3 * s; }
  

  /* choose pivot - or die */
  if (fabs(r3[2])>fabs(r2[2])) SWAP_ROWS(r3, r2);
  if (0.0 == r2[2])  return false;
  
  /* eliminate third variable */
  m3 = r3[2]/r2[2];
  r3[3] -= m3 * r2[3]; 
  r3[4] -= m3 * r2[4];
  r3[5] -= m3 * r2[5]; 
  r3[6] -= m3 * r2[6];
  r3[7] -= m3 * r2[7];
  
  /* last check */
  if (0.0 == r3[3]) return false;
  
  s = 1.0/r3[3];              /* now back substitute row 3 */
  r3[4] *= s; r3[5] *= s; r3[6] *= s; r3[7] *= s;
  
  m2 = r2[3];                 /* now back substitute row 2 */
  s  = 1.0/r2[2];
  r2[4] = s * (r2[4] - r3[4] * m2); r2[5] = s * (r2[5] - r3[5] * m2);
  r2[6] = s * (r2[6] - r3[6] * m2); r2[7] = s * (r2[7] - r3[7] * m2);
  m1 = r1[3];
  r1[4] -= r3[4] * m1; r1[5] -= r3[5] * m1;
  r1[6] -= r3[6] * m1; r1[7] -= r3[7] * m1;
  m0 = r0[3];
  r0[4] -= r3[4] * m0; r0[5] -= r3[5] * m0;
  r0[6] -= r3[6] * m0; r0[7] -= r3[7] * m0;
  
  m1 = r1[2];                 /* now back substitute row 1 */
  s  = 1.0/r1[1];
  r1[4] = s * (r1[4] - r2[4] * m1); r1[5] = s * (r1[5] - r2[5] * m1);
  r1[6] = s * (r1[6] - r2[6] * m1); r1[7] = s * (r1[7] - r2[7] * m1);
  m0 = r0[2];
  r0[4] -= r2[4] * m0; r0[5] -= r2[5] * m0;
  r0[6] -= r2[6] * m0; r0[7] -= r2[7] * m0;
  
  m0 = r0[1];                 /* now back substitute row 0 */
  s  = 1.0/r0[0];
  r0[4] = s * (r0[4] - r1[4] * m0); r0[5] = s * (r0[5] - r1[5] * m0);
  r0[6] = s * (r0[6] - r1[6] * m0); r0[7] = s * (r0[7] - r1[7] * m0);
  
  M(0,0) = r0[4]; M(0,1) = r0[5];
  M(0,2) = r0[6]; M(0,3) = r0[7];
  M(1,0) = r1[4]; M(1,1) = r1[5];
  M(1,2) = r1[6]; M(1,3) = r1[7];
  M(2,0) = r2[4]; M(2,1) = r2[5];
  M(2,2) = r2[6]; M(2,3) = r2[7];
  M(3,0) = r3[4]; M(3,1) = r3[5];
  M(3,2) = r3[6]; M(3,3) = r3[7]; 
  
  return true;
#undef SWAP_ROWS
}


//-----------------------------------------------------------------------------


#undef MAT
#undef M


//=============================================================================
} // namespace Math
//=============================================================================
#endif // !defined(_MATRIX4X4_C)
