/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
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
License along with this software.*/

#ifndef _GEN_HOM_H
#define _GEN_HOM_H

#include "DataDefinitions.h"

///////////////////////////////////////
template<typename Scalar> class genHomography
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  genHomography(Scalar* K_)
    :
  Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
  Kr(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    MC(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0), mc(0.0,0.0,0.0), 
    centers(NULL), normals(NULL), tra(NULL), rot(NULL)
  {init();};

// assume Kr == Kl
  genHomography(Scalar* K_, Scalar* MC_, Scalar* mc_)
    :
  Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
  Kr(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
    MC(MC_[0], MC_[3], MC_[6], MC_[1], MC_[4], MC_[7], MC_[2], MC_[5], MC_[8]), mc(mc_[0], mc_[1], mc_[2]),
    centers(NULL), normals(NULL), tra(NULL), rot(NULL)
  {init();};

  ///when going from left to right camera Kr has to be set if Kr != Kl that is to say
  genHomography(Scalar* K_, Scalar* MC_, Scalar* mc_, Scalar* Kr_)
    : 
  Kl(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]),
  Kr(Kr_[0], Kr_[3], Kr_[6], Kr_[1], Kr_[4], Kr_[7], Kr_[2], Kr_[5], Kr_[8]),
    MC(MC_[0], MC_[3], MC_[6], MC_[1], MC_[4], MC_[7], MC_[2], MC_[5], MC_[8]), mc(mc_[0], mc_[1], mc_[2]),
    centers(NULL), normals(NULL), tra(NULL), rot(NULL)
  {init();};

  ~genHomography(){};

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

  M3 getHom(const int index, const int pIndex = -1)
  {
    if (rot == NULL)
      return getHom_N(index);
    else
      return getHom_rot(index, pIndex);
  }

  M3 getI(const int index)
  {
    return A1;
  }

  M3 getHom_rot(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    P3 n( (*normals)[index] );	
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

#ifdef _use_patchCenters_
    center = center * (1./(n|center));
    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...
#endif
    P3 a2 = a1 + A1*t;

    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

  M3 getHom_N(const int index)
  {
    P3 n ( (*normals)[index] );

    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;
  };

  // one for all
  M3 getHom_Normal( P3 n )
  {
    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;
  };

  /// wrt the second view:
  M3 getHom_n(const int index, const int pIndex = -1)
  {
    if (rot == NULL)
      return getHom_N_n(index);
    else
      return getHom_rot_n(index, pIndex);
  }

  /// wrt the second view:
  M3 getHom_rot_n(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    P3 n( -(*normals)[index] );	//-: after normalflip() - delivers usual/correct normals for homography
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );

#ifdef _use_patchCenters_
    center = center * (1./(n|center));
    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...
#endif

    P3 a2 = a1 + A1*t;
    return (A1*R - M3( a2[0]*n[0], a2[0]*n[1], a2[0]*n[2],  a2[1]*n[0], a2[1]*n[1], a2[1]*n[2],  a2[2]*n[0], a2[2]*n[1], a2[2]*n[2]))*iK;
  };

  /// wrt the second view:
  M3 getHom_N_n(const int index)
  {
    P3 n (-(*normals)[index] );//-: after normalflip() - delivers usual/correct normals for homography
    return (A1 - M3( a1[0]*n[0], a1[0]*n[1], a1[0]*n[2],  a1[1]*n[0], a1[1]*n[1], a1[1]*n[2],  a1[2]*n[0], a1[2]*n[1], a1[2]*n[2]))*iK;// from pixel to camera, but already in camera coords
  };

  /// normals in the second view:
  P3 getViewNormal (const int index, const int pIndex)
  {
    if (rot==NULL)  
      return getView_N(index);
    else
      return getView_rN(index, pIndex);
  }

  P3 getView_N( int index )
  {
    P3 n ( (*normals)[index] );
    P3 rn = MC*n;
    Scalar rnt = 1. - (mc|rn);rn/=rnt;
    return rn;
  }
  
  P3 getView_rN(const int index, const int pIndex)
  {
    P3 center ( &centers[3 * pIndex] ); // the center depends on the patch not the label ARHG

    P3 n( (*normals)[index] );	// guess : after flip this delivers 3d coords
    const M3& R( (*rot)[index] );
    P3 t( (*tra)[index] );
#ifdef _use_patchCenters_
    center = center * (1./(n|center));
    t -= center - R*center;	// changed to react on changed situation -> also input is -center instead of t -= ...
#endif
    ////
    P3 rn = R*n;
    Scalar rnt = 1. + (n|center) - (rn|(t+center));
    rn /= rnt;
    //

    rn = MC*rn;
    rnt = 1. - (mc|rn);rn/=rnt;
    return rn;
  };

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
};

///////////////////////////////////////


#endif