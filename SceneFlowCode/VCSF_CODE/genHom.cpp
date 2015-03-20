/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

//////////////////////////////////////////////////////////////////////////////
// deliver homograpies and normals - wrt to a certain view and timestep     //
//////////////////////////////////////////////////////////////////////////////

// with _flipMC_ still both possible - =-O
//#define _flipCorrect_
//#define _flipMC_
// all the same ; either nothing or _flipCorrect_ or _flipCorrect_ and _flipMC_

#ifndef __genHom__cpp
#define __genHom__cpp

#include "genHom.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;
///////////////////////////////////////


/// does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled: p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::getView_N( int index )
{
    P3 n ( (*normals)[index] );
    P3 rn = MC*n;
    Scalar rnt = 1. - (mc|rn);rn/=rnt;
    return rn;
}
 
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::getView_NC0_paper( int index )
{
    P3 n ( (*normals)[index] );
    P3 rn = n;
    // well should be 1 - m|n, since n|p_3D =-1 -> n_r * (p+m) =  -1+ m|n/ (1- m|n) = -1
    Scalar rnt = 1. + (mc|rn);rn/=rnt; 
    return rn;
}

  // has some problems with general cameras ?! mc positive or negative: does not matter here
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::
getView_NC0_new( int index )
  {
#ifndef _flipCorrect_

    P3 n ( -(*normals)[index] );
    P3 q1(0.0,0.0,1.0);
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n))
      q2 = q22;
    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n))
      q3=q33;

    // qs are some point on image plane and (by division in 3D) + mc -> 3D point in second cam
    q1 /= (q1|n); q1+= mc;
    q2 /= (q2|n); q2+= mc;
    q3 /= (q3|n); q3+= mc;
    M3 mat(q1, q2, q3);mat=mat.transpose();mat.invert();
    P3 rn = mat*P3(1.,1.,1.);
    return -rn;
#else
    // WELL it is exactly the same but leads to completely different results:
    // assume p on plane : p|n_l =1
    // now n_r = n_l / (1+ m|n_l)
    // then (p+m) | n_r = 1
    // this is CORRECT
    P3 n ( (*normals)[index] );
    P3 q1(0.0,0.0,1.0);
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n))
      q2 = q22;
    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n))
      q3=q33;

    // qs are some point on image plane and (by division in 3D) + mc -> 3D point in second cam
#ifndef _flipMC_
    q1 /= (q1|n); q1+= mc;
    q2 /= (q2|n); q2+= mc;
    q3 /= (q3|n); q3+= mc;
#else
    q1 /= (q1|n); q1-= mc;
    q2 /= (q2|n); q2-= mc;
    q3 /= (q3|n); q3-= mc;
#endif
    M3 mat(q1, q2, q3);mat=mat.transpose();mat.invert();
    P3 rn = mat*P3(1.,1.,1.);
    return rn;
#endif

//    $(x^t \overline{\mathbf{n}}^1 = 1 \leftrightarrow $
//$(R^{-1}x-R^{-1}t)^t \overline{\mathbf{n}} = 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} - R^{-1}t)^t \overline{\mathbf{n}}= 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} = 1 + R^{-1}t)^t \overline{\mathbf{n}} \leftrightarrow $
//$ {\mathbf{n}}^1 = (R^{-1}^t \overline{\mathbf{n}} / (1 + R^{-1}t)^t \overline{\mathbf{n}})$.
  }
    
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::
getView_rNC0_new( const int index )
  {
#ifndef _flipCorrect_
    P3 n( -(*normals)[index] );	// guess : after flip this delivers 3d coords
#else
    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
#endif
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );
    ////
    P3 q1(0.0,0.0,1.0);

    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n)) q2 = q22;

    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n)) q3 = q33;

#ifndef _flipMC_
    q1 /= (q1|n);
    q1 = MC*(R*q1+t)+mc;
    q2 /= (q2|n);
    q2 = MC*(R*q2+t)+mc;
    q3 /= (q3|n);
    q3 = MC*(R*q3+t)+mc;
#else
    q1 /= (q1|n);
    q1 = MC*(R*q1+t)-mc;
    q2 /= (q2|n);
    q2 = MC*(R*q2+t)-mc;
    q3 /= (q3|n);
    q3 = MC*(R*q3+t)-mc;
#endif


    M3 mat(q1, q2, q3);mat=mat.transpose();mat.invert();
    P3 rn = mat*P3(1.,1.,1.);
#ifdef _DEBUG
    // probe: should be one
    Scalar diff1 = (rn|q1);
    Scalar diff2 = (rn|q2);
    Scalar diff3 = (rn|q3);
#endif
#ifndef _flipCorrect_
    return -rn;
#else
    return rn;
#endif
  };

template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::
getView_rNC0_new_x2( const int index )
  {
#ifndef _flipCorrect_
    P3 n( -(*normals)[index] );	// guess : after flip this delivers 3d coords
#else
    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
#endif
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    // twice and global adjustment:
    M3 R2  = R*R;
    P3 t2  = t + (R*t);

    t2  = R_Glob * t2 + T_Glob; 
    R2  = R_Glob * R2;

    ////
    P3 q1(0.0,0.0,1.0);
    
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n))
      q2 = q22;
    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n))
      q3=q33;

#ifndef _flipMC_
    q1 /= (q1|n);q1 = MC*(R2*q1+t2)+mc;
    q2 /= (q2|n);q2 = MC*(R2*q2+t2)+mc;
    q3 /= (q3|n);q3 = MC*(R2*q3+t2)+mc;
#else
    q1 /= (q1|n);q1 = MC*(R2*q1+t2)-mc;
    q2 /= (q2|n);q2 = MC*(R2*q2+t2)-mc;
    q3 /= (q3|n);q3 = MC*(R2*q3+t2)-mc;
#endif

    M3 mat(q1, q2, q3);mat=mat.transpose();mat.invert();
    P3 rn = mat*P3(1.,1.,1.);

#ifndef _flipCorrect_
    return -rn;
#else
    return rn;
#endif
  };

#ifdef _outCommented
  P3 getView_rNC0_new_x2( const int index )
  {
    P3 n = - getView_rNC0_new( index ); // or + ????

    M3 R( (*rot)[index] );
    P3 t( (*tra)[index] );

    // twice and global adjustment:
    t   = R_Glob * t + T_Glob; 
    R   = R_Glob * R;

    ////
    P3 q1(0.0,0.0,1.0);q1 /= (q1|n);q1 = MC*(R*q1+t)+mc;
    P3 q2(1.0,0.0,1.0);q2 /= (q2|n);q2 = MC*(R*q2+t)+mc;
    P3 q3(0.0,1.0,1.0);q3 /= (q3|n);q3 = MC*(R*q3+t)+mc;
    M3 mat(q1, q2, q3);mat=mat.transpose();mat.invert();
    P3 rn = mat*P3(1.,1.,1.);
    return -rn;
  };
#endif

template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::getView_rNC0_new_inv( const int index )
  {
#ifndef _flipCorrect_
    P3 n( -(*normals)[index] );	// guess : after flip this delivers 3d coords
#else
    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
#endif
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    // inverse:
    M3 Ri = R;Ri = Ri.transpose();
    P3 ti = -(Ri*t);

    ti += Ri * T_Glob;
    Ri  = Ri * R_Glob;

    // or: test by experiment
//    ti  = T_Glob + R_Glob * ti; 
//    Ri  = R_Glob * Ri;

    // Rt_inv(:,:,i) =  inv(Rt_lin(:,:,i)) * Rt_cam;
    // or     Rt_inv(:,:,i) =  Rt_cam * inv(Rt_lin(:,:,i)) ; % both are fine from test - one is true though
    ////
    P3 q1(0.0,0.0,1.0);
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n))
      q2=q22;
    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n))
      q3=q33;
#ifndef _flipMC_
    q1 /= (q1|n);q1 = MC*(Ri*q1+ti)+mc;
    q2 /= (q2|n);q2 = MC*(Ri*q2+ti)+mc;
    q3 /= (q3|n);q3 = MC*(Ri*q3+ti)+mc;
#else
    q1 /= (q1|n);q1 = MC*(Ri*q1+ti)-mc;
    q2 /= (q2|n);q2 = MC*(Ri*q2+ti)-mc;
    q3 /= (q3|n);q3 = MC*(Ri*q3+ti)-mc;
#endif
    M3 mat(q1, q2, q3);mat=mat.transpose();
    
    if ( !mat.invert() )
    {
      printf("Fail inversion ");
    }

    P3 rn = mat*P3(1.,1.,1.);
#ifndef _flipCorrect_
    return -rn;
#else
    return rn;
#endif
  };

template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::getView_rN(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    center = center * (1./(n|center));
    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...

    ////
    P3 rn = R*n;
    Scalar rnt = 1. + (n|center) - (rn|(t+center));
    rn /= rnt;
    //
    rn = MC*rn;
    rnt = 1. - (mc|rn);rn/=rnt;
    return rn;
  };

  /// that is fine - ok
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::
getView_rNC0( const int index )
  {
    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );
    ////
    P3 rn = R*n;
    Scalar rnt = 1. - (rn|(t));// should be +  no eq is n^t p = -1 for me - LOL
    rn /= rnt;
    //
    rn = MC*rn;
    rnt = 1. - (mc|rn);rn/=rnt;// and this ?

#ifdef _DEBUG
// probe
    P3 q1(0.0,0.0,1.0);
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n)) q2 = q22;

    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n)) q3 = q33;

    q1 /= fabs(q1|n);
    q1 = MC*(R*q1+t)+mc;
    q2 /= fabs(q2|n);q2 = MC*(R*q2+t)+mc;
    q3 /= fabs(q3|n);q3 = MC*(R*q3+t)+mc;

    Scalar diff1 = (rn|q1);
    Scalar diff2 = (rn|q2);
    Scalar diff3 = (rn|q3);
#endif
    return rn;
  };

#ifdef _DEBUG
  /// wrong
template<class Scalar>
typename Math::VectorT<Scalar, 3>
genHomoG<Scalar>::
getView_rNC0P( const int index )
{
    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );
    ////
    P3 rn = R*n;
    Scalar rnt = -1. + (rn|(t));// should be + ?
    rn /= rnt;
    //
    rn = MC*rn;
    rnt = -1. + (mc|rn);rn/=rnt;// and this ?

    // probe
    P3 q1(0.0,0.0,1.0);
    P3 q2(1.0,0.0,1.0);P3 q22(1.0,0.0,0.0);
    if ( fabs(q2|n) < fabs(q22|n)) q2 = q22;

    P3 q3(0.0,1.0,1.0);P3 q33(0.0,1.0,0.0);
    if ( fabs(q3|n) < fabs(q33|n)) q3 = q33;

    q1 /= fabs(q1|n);
    q1 = MC*(R*q1+t)+mc;
    q2 /= fabs(q2|n);q2 = MC*(R*q2+t)+mc;
    q3 /= fabs(q3|n);q3 = MC*(R*q3+t)+mc;

    Scalar diff1 = (rn|q1);
    Scalar diff2 = (rn|q2);
    Scalar diff3 = (rn|q3);

    return rn;
  };
#endif

//$(x^t \overline{\mathbf{n}}^1 = 1 \leftrightarrow $
//$(R^{-1}x-R^{-1}t)^t \overline{\mathbf{n}} = 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} - R^{-1}t)^t \overline{\mathbf{n}}= 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} = 1 + R^{-1}t)^t \overline{\mathbf{n}} \leftrightarrow $
//$ {\mathbf{n}}^1 = (R^{-1}^t \overline{\mathbf{n}} / (1 + R^{-1}t)^t \overline{\mathbf{n}})$.
//$ {\mathbf{n}}^1 = (R \overline{\mathbf{n}} / (1 + R^{-1}t)^t \overline{\mathbf{n}})$.
//$ {\mathbf{n}}^1 = (R \overline{\mathbf{n}} / (1 + t^t  R\overline{\mathbf{n}})$.

/*
template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rot_n(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    P3 n( -(*normals)[index] );	//-: after normalflip() - delivers usual/correct normals for homography
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    center = center * (1./(n|center));
    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...
    P3 a2 = a1 + A1*t;
    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_N_n(const int index)
  {
    P3 n (-(*normals)[index] );//-: after normalflip() - delivers usual/correct normals for homography
    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;// from pixel to camera, but already in camera coords
  };
*/

template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rot(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    //    printf("center2d %d: (%f,%f,%f)\n", index, (center)[0], (center)[1], (center)[2] );

    P3 n( (*normals)[index] );	
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    center = center * (1./(n|center));

    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...

    P3 a2 = a1 + A1*t;

    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rotC0( const int index )
  {
    P3 n( (*normals)[index] );	
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    P3 a2 = a1 + A1*t;

    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

  // TODO check this
template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rotC0_x2( const int index )
  {
    P3 n( (*normals)[index] );	
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] ); 

    // twice and global adjustment: REALLY TWICE
    M3 R2  = R*R;
    P3 t2  = t + (R*t);

    t2  = R_Glob * t2 + T_Glob; 
    R2  = R_Glob * R2;

    P3 a2 = a1 + A1*t2;

    return (A1*R2 - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

  /*
  // Should work (unless n1 is negative or so)
template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rotC0_x2( const int index )
  {
    M3 R( (*rot)[index] );
    P3 t( (*tra)[index] );

    // twice and global adjustment: REALLY TWICE
    M3 H1 = getHom_rotC0( index );
    P3 n  = getView_rNC0_new( index );

    t   = R_Glob * t + T_Glob; 
    R   = R_Glob * R;

    // just multiply both homographies :)

    P3 a2 = a1 + A1*t;

    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK*H1;
  };
  */

template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_rotC0_inv( const int index )
  {
    P3 n( (*normals)[index] );	
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

    // inverse:
    M3 Ri = R;Ri = Ri.transpose();
    P3 ti = -(Ri*t);

    ti += Ri * T_Glob; 
    Ri  = Ri * R_Glob;

    // or: test by experiment
//    ti  = T_Glob + R_Glob * ti; 
//    Ri  = R_Glob * Ri;

    P3 a2 = a1 + A1*ti;

    return (A1*Ri - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_N(const int index)
  {
    P3 n ( (*normals)[index] );
    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;
  };

// one for all
template<class Scalar>
Math::Mat3x3T<Scalar>
genHomoG<Scalar>::
getHom_Normal( const P3& n )
  {
    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;
  };

#endif