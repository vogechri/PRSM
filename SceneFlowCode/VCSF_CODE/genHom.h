/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

This software is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this software.
*/

#ifndef __GEN_HOM_HH
#define __GEN_HOM_HH

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;

///////////////////////////////////////
/*!
Class produces homographies and normals as observed from the given view and timestep.
Homographies are always given from canonical view to others.
The other view is specified by:
- inverseMotion: one timestep back
- doubleMotion: two timesteps in the future
- no MC_, mc_ given: 2nd view is wrt the same camera, else into cam specified (by MC,mc)
- no rotation/translation given : transformation specified by normal only
*/
template<typename Scalar> class genHomoG
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;

  // default constructor to allow a vector
  genHomoG() : centers(NULL), normals(NULL), tra(NULL), rot(NULL), inverseMotion(false), doubleMotion(false)
  {};

  genHomoG(Scalar* K_)
    : Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    Kr(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    MC(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0), mc(0.0,0.0,0.0), 
    centers(NULL), normals(NULL), tra(NULL), rot(NULL), inverseMotion(false), doubleMotion(false)
  {init();};

// assume Kr == Kl
  genHomoG(Scalar* K_, Scalar* MC_, Scalar* mc_)
    : Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    Kr(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    MC(MC_[0], MC_[3], MC_[6], MC_[1], MC_[4], MC_[7], MC_[2], MC_[5], MC_[8]), mc(mc_[0], mc_[1], mc_[2]),
    centers(NULL), normals(NULL), tra(NULL), rot(NULL), inverseMotion(false), doubleMotion(false)
  {init();};

  ///when going from left to right camera Kr has to be set if Kr != Kl that is to say
  genHomoG(Scalar* K_, Scalar* MC_, Scalar* mc_, Scalar* Kr_)
    : Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    Kr(Kr_[0], Kr_[3], Kr_[6], Kr_[1], Kr_[4], Kr_[7], Kr_[2], Kr_[5], Kr_[8]),
    MC(MC_[0], MC_[3], MC_[6], MC_[1], MC_[4], MC_[7], MC_[2], MC_[5], MC_[8]), mc(mc_[0], mc_[1], mc_[2]),
    centers(NULL), normals(NULL), tra(NULL), rot(NULL), inverseMotion(false), doubleMotion(false)
  {init();};

  ~genHomoG(){};

  size_t getSize() {return normals->size();};

  void init()
  {
    iK = Kl; iK.invert();
    A1 = Kr*MC;
    a1 = Kr*mc;
  }

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<P3>* tra_, const std::vector<M3>* rot_ ) {tra = tra_;rot = rot_;}

  void setCenters( Scalar* centers_) {centers = centers_;}

  /// get homos, normals for 3rd frame: one timestep into the past
  void do_inverse_motion(bool do_it, M3 globMat, P3 globTra )
  {inverseMotion = do_it; T_Glob = globTra; R_Glob = globMat;}

  /// 4 frame version: to get the mappings and normals two frames in the future, so from t to t+2
  void do_double_motion(bool do_it, M3 globMat, P3 globTra )
  {
    doubleMotion = do_it; 
    T_Glob = globTra; 
    R_Glob = globMat;

    R_Glob = R_Glob.transpose();
    T_Glob = -( R_Glob * T_Glob );
  }

  /*! return the normal as seen from the view, moving plane proposal given by index
  Normal assuming rotation centered at central pixel of patch (representation not used in view-consistency)
  */
  P3 getViewNormal (const int index, const int pIndex)
  {
    if (rot==NULL)  
      return getView_N(index);
    else
      return getView_rN(index, pIndex);
  }

  /*! return the normal as seen from the view, moving plane proposal given by index.
  Assuming camera centered rotations. (Default for view consistency)
  */
  P3 getViewNormalC0 (const int index)
  {
    P3 norm(0,0,1);
    if (rot==NULL)
//      norm = getView_NC0_paper(index);
      norm = getView_NC0_new(index);
    else
      if (!inverseMotion)
      {
        if (!doubleMotion)
          norm = getView_rNC0_new(index);
        else
          norm = getView_rNC0_new_x2(index);
      }
      else
        norm = getView_rNC0_new_inv(index);

#ifdef _DEBUG
        Scalar dn_ =  1. / sqrt(norm|norm);
        if ( dn_ != dn_ || std::numeric_limits<Scalar>::infinity()==dn_ )
          printf("Error in viewnormal too small : %f  id: %d \n", dn_, index );
#endif
        return norm;
  }

  /*! does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled
      p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
  */
  P3 getView_N( int index );

  P3 getView_NC0_paper( int index );

  // has some problems with general cameras ?! mc positive or negative: does not matter here
  P3 getView_NC0_new( int index );

  P3 getView_rNC0_new( const int index );

  /// two timesteps into the future
  P3 getView_rNC0_new_x2( const int index );

  /// one timestep into the past
  P3 getView_rNC0_new_inv( const int index );

  P3 getView_rN(const int index, const int pIndex);

  /// that is fine - ok
  P3 getView_rNC0( const int index );

#ifdef _DEBUG
  P3 getView_rNC0P( const int index );
#endif

//$(x^t \overline{\mathbf{n}}^1 = 1 \leftrightarrow $
//$(R^{-1}x-R^{-1}t)^t \overline{\mathbf{n}} = 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} - R^{-1}t)^t \overline{\mathbf{n}}= 1 \leftrightarrow $
//$x(R^{-1}^t \overline{\mathbf{n}} = 1 + R^{-1}t)^t \overline{\mathbf{n}} \leftrightarrow $
//$ {\mathbf{n}}^1 = (R^{-1}^t \overline{\mathbf{n}} / (1 + R^{-1}t)^t \overline{\mathbf{n}})$.
//$ {\mathbf{n}}^1 = (R \overline{\mathbf{n}} / (1 + R^{-1}t)^t \overline{\mathbf{n}})$.
//$ {\mathbf{n}}^1 = (R \overline{\mathbf{n}} / (1 + t^t  R\overline{\mathbf{n}})$.


  /*! return the homography from the canonical view to the one encoded.
  Moving plane proposal given by index.
  Assuming segment centered rotations. (not used for view consistency)
  */
  M3 getHom(const int index, const int pIndex = -1)
  {
    if (rot == NULL)
      return getHom_N(index);
    else
      return getHom_rot(index, pIndex);
  }

  /*! return the homography from the canonical view to the one encoded.
  Moving plane proposal given by index.
  Assuming camera centered rotations. (used for view consistency)
  */
  M3 getHomC0( const int index )
  {
    if (rot == NULL)
        return getHom_N(index);
    else
      if(inverseMotion)
        return getHom_rotC0_inv(index);
      else
      {
        if (!doubleMotion)
          return getHom_rotC0(index);
        else
          return getHom_rotC0_x2(index);
      }
  }

  M3 getI(const int index)  { return A1; }

  /// one timestep into the future, segment centered rotations
  M3 getHom_rot(const int index, const int pIndex);

  /// one timestep into the future, view centered rotations
  M3 getHom_rotC0( const int index );

  /// two time-steps into the future
  M3 getHom_rotC0_x2( const int index );

  /// one timestep into the past
  M3 getHom_rotC0_inv( const int index );

  M3 getHom_N(const int index);
  // one for all
  M3 getHom_Normal( const P3& n );

/*
  M3 getHom_n(const int index, const int pIndex = -1)
  {
    if (rot == NULL)
      return getHom_N_n(index);
    else
      return getHom_rot_n(index, pIndex);
  }

  M3 getHom_rot_n(const int index, const int pIndex);

  M3 getHom_N_n(const int index);
*/

private:

  M3 iK;
  M3 Kl;
  M3 Kr;
  M3 MC;
  P3 mc;
  M3 A1;
  P3 a1;

  /// centers of the pathces
  Scalar* centers;

  const std::vector<P3>* normals;
  /// translations
  const std::vector<P3>* tra;
  /// rotations
  const std::vector<M3>* rot;

  bool inverseMotion;

  bool doubleMotion;

  P3 T_Glob;
  M3 R_Glob;
};

#if defined(INCLUDE_TEMPLATES) && !defined(__genHom__cpp)
#include "genHom.cpp"
#endif

#endif // __GEN_HOM_HH