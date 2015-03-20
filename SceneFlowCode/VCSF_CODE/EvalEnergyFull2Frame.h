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

///////////////////////////////////////////////////////
//////////      Per segment smoothness      ///////////
///////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANOR__
#define __ENERGY_ROTTRANOR__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
// e.g. defines if 3D smoothing is enabled
#include "DataDefinitionsVC.h"

#include <map>
#include <vector>

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

/// provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
template<typename Scalar> class EvalEnergyFullFrame
{
public:

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 4> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::Mat3x3T<Scalar>    M3; 

  EvalEnergyFullFrame(): gamma(1.0), epsilon(0.00), normals(NULL), tra(), rot(), rotJump(20.), depthJump(20.), nSegments(0),
    Pl(1.,0.,0., 0.,1.,0., 0.,0.,1.), Pr(1.,0.,0., 0.,1.,0., 0.,0.,1.), pr(0.,0.,0.)
  {};

  ~EvalEnergyFullFrame(){};

  void set2dMotionMatrix ( Scalar* Kl, Scalar* Kr, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(Kl);
    Pr = M3(Kr) * M3(Rot);
    pr = M3(Kr) * P3(mc);
  };

  void set2dMotionMatrix_inv ( Scalar* K, Scalar* Kr, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(K);
    Pr = M3(Kr) * M3(Rot).invert();
    pr = Pr * (-P3(mc));
  }

  void setGamma (Scalar gamma_)               {gamma = gamma_;}

  void setRotWeight (Scalar rotWeight_)       { rotWeight = rotWeight_;};

  void setRotJump   (Scalar rotJump_)         {rotJump = rotJump_;}
  void setDepthJump (Scalar depthJump_)       {depthJump = depthJump_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  const std::vector<int>& getIdl()  { return Idl; }
  const std::vector<int>& getIdk()  { return Idk; }

  const std::vector<Scalar>& getF11()  { return F11; }
  const std::vector<Scalar>& getF00()  { return F00; }
  const std::vector<Scalar>& getF10()  { return F10; }
  const std::vector<Scalar>& getF01()  { return F01; }

  void setHalfEdges( Scalar* _halfEdgeX, Scalar* _halfEdgeY) {halfEdgeX = _halfEdgeX; halfEdgeY = _halfEdgeY;}
  /// weights on the cross egges in a 8 neighbourhood
  void setCrossHalfEdges( Scalar* _halfEdgeXY, Scalar* _halfEdgeiXY) {halfEdgeXY = _halfEdgeXY; halfEdgeiXY = _halfEdgeiXY;}

  /// segimg must be flipped
  void buildWeights(int* segImg, int width, int height, int nSegments);

  /// segimg must not! be flipped
  void buildWeightsVC(int* segImg, int width, int height, int nSegments);

  void prepare( int nSegments_, const mxArray* edges_, const mxArray* weights_, Scalar* centers_ );

  /// segimage flipped so also width and height
  void prepareOwnWeights( int nSegments_, const mxArray* edges_, Scalar* centers_, int* segImg, int width, int height);

  size_t edges_num() {return Idk.size();};

#ifdef __USE3D__
  Scalar compute_score_Fuse( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 );
#else
  Scalar compute_score_Fuse3D( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 );
#endif

#ifdef __USE3D__
  Scalar compute_score_Fuse2D( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 );
#else
  Scalar compute_score_Fuse( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 );
#endif

  /// faster version if only the current solution is to be evaluated - called from getEnergy , etc.
#ifdef __USE3D__
  Scalar compute_score_Fuse2D00( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 );
#endif
  Scalar compute_score_Fuse00( std::vector<int>& curr_ids );

  /// view consistent function
  Scalar compute_score_Fuse_local( const std::vector<int>& curr_ids, const std::vector<int>& curr_ids2, const std::vector<int>& localIds, 
    std::vector<Scalar>& _F00, std::vector<Scalar>& _F01, std::vector<Scalar>& _F10, std::vector<Scalar>& _F11,
    std::vector<int>& _idk, std::vector<int>& _idl, bool fillVecs=true );

private:
  
  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2dProjections( const P3& p_t0, P2& pixP_l, P2& pixP_r) const 
  {
    P3 pixP_l_3 = Pl * p_t0;     pixP_l = P2(pixP_l_3[0], pixP_l_3[1]); pixP_l  /= pixP_l_3[2];
    P3 pixP_r_3 = Pr * p_t0 +pr; pixP_r = P2(pixP_r_3[0], pixP_r_3[1]); pixP_r  /= pixP_r_3[2];
  }

  inline void get2dDiff_Fast( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                          const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                          P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const 
  {
    P2 difft0   = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P2 difft1   = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
       flowDiff = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }

  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2dDiff( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff )
  {
    // that is constant
    P3 pixP_t0_l = Pl * p_t0;pixP_t0_l /= pixP_t0_l[2];
    P3 pixQ_t0_l = Pl * q_t0;pixQ_t0_l /= pixQ_t0_l[2];

    P3 pixP_t1_l = Pl * p_t1;pixP_t1_l /= pixP_t1_l[2];
    P3 pixQ_t1_l = Pl * q_t1;pixQ_t1_l /= pixQ_t1_l[2];

    // if Pr =Pl that is also const, if pr = x,0,0 even more
    P3 pixP_t0_r = Pr * p_t0 +pr;pixP_t0_r /= pixP_t0_r[2];
    P3 pixQ_t0_r = Pr * q_t0 +pr;pixQ_t0_r /= pixQ_t0_r[2];

    P3 pixP_t1_r = Pr * p_t1 +pr; pixP_t1_r /= pixP_t1_r[2];
    P3 pixQ_t1_r = Pr * q_t1 +pr; pixQ_t1_r /= pixQ_t1_r[2];

    P3 difft0 = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P3 difft1 = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
    P3 diffF  = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    flowDiff   = P2(diffF[0],  diffF[1]);
    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }

   inline void get2dDiff_vc( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                             const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                             P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const 
  {
    P2 difft0   = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P2 difft1   = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
       flowDiff = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

//    flowDiff   = P2(diffF[0],  diffF[1]);
    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }
  //////////////////////////////////////////////////

  M3 Pl;//R=I t=0

  M3 Pr;//K*R
  P3 pr;//K*t

  /// helper
  std::vector<std::map<int, Scalar> > edgeWeightMap;

  /// the edges 0,1,2 are 2d coordinates of point p 
  std::vector<P3> edgesP;
  /// the edges 0,1,2 are 2d coordinates of point q
  std::vector<P3> edgesQ;

  /// the centers of the pathces in 2d
  std::vector<P3> centersA;
  std::vector<P3> centersB;

  /// the weights of the edges
  std::vector<Scalar> weights;

  /// ids of the patches id k has a corresponding neighbour l in Idl
  std::vector<int> Idk;

  /// ids of the patches id l has a corresponding neighbour k in Idl
  std::vector<int> Idl;

  const std::vector<P3>* normals;

  /// translations
  const std::vector<P3>* tra;
  /// rotations
  const std::vector<M3>* rot;


  std::vector<Scalar> F00;
  std::vector<Scalar> F01;
  std::vector<Scalar> F10;
  std::vector<Scalar> F11;
  
  /// weights on the edges per pixel (4/8-neighbourhood)
  Scalar* halfEdgeX;
  Scalar* halfEdgeY;

  /// weight on the cross edges (8-neighbourhood)
  Scalar* halfEdgeXY;
  Scalar* halfEdgeiXY;

  Scalar gamma;
  double epsilon;
  Scalar rotJump;
  Scalar depthJump;
  Scalar rotWeight;

  int nSegments;
};
#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANOR__cpp)
#include "EvalEnergyFull2Frame.cpp"
#endif

#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANORVC__cpp)
#include "EvalEnergyFull2FrameVC.cpp"
#endif
#endif
