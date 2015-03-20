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
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANORPIXEL__
#define __ENERGY_ROTTRANORPIXEL__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

//#define __use_m128__

#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include <map>
#include <vector>

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

// uses a field which maps from pixelids to ids in the mrf

/// provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
template<typename Scalar> class EvalEnergyFullFramePixel
{
public:

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 4> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::VectorT<int, 4>    P4i;

  typedef Math::Mat3x3T<Scalar>    M3; 

  EvalEnergyFullFramePixel(int _width, int _height)
    : gamma(1.0), temporalWeight(0.03), epsilon(0.00), normals(NULL), tra(), rot(), 
    rotJump(15.), tempJump(100.), segImg(NULL), depthJump(sqrt(1300.)), nSegments(0), weightMap(), 
    width(_width), height(_height), centers(NULL), iK(1,0,0,0,1,0,0,0,1), iiKt(1,0,0,0,1,0,0,0,1),
    halfEdgeX(NULL), halfEdgeY(NULL), halfEdgeXY(NULL), halfEdgeiXY(NULL), edgeScale(0.1) //, epix(0)
    //, p2dSet(NULL)
  {};

  ~EvalEnergyFullFramePixel(){};

  void setEdgeScale (Scalar _edgeScale) {edgeScale = _edgeScale;};
  
  void set2dMotionMatrix ( Scalar* K, Scalar* Rot, Scalar* mc, Scalar pixJump_ );
  void set2dMotionMatrix_inv ( Scalar* K, Scalar* Rot, Scalar* mc, Scalar pixJump_ );
  void set2dMotionMatrix ( Scalar* Kl, Scalar* Kr, Scalar* Rot, Scalar* mc, Scalar pixJump_ );

  void setRotWeight (Scalar rotWeight_)         {rotWeight = rotWeight_;};

  void setWidthHeight (int _width, int _height) {width=_width; height=_height;}

  void setEpsilon (Scalar epsilon_)             {epsilon = epsilon_;}

  void setK (Scalar* K_);

  /// weights on the edges, based on the input image
  void setHalfEdges( Scalar* _halfEdgeX, Scalar* _halfEdgeY) {halfEdgeX = _halfEdgeX; halfEdgeY = _halfEdgeY;}
  /// weights on the cross egges in a 8 neighbourhood
  void setCrossHalfEdges( Scalar* _halfEdgeXY, Scalar* _halfEdgeiXY) {halfEdgeXY = _halfEdgeXY; halfEdgeiXY = _halfEdgeiXY;}

  void setGamma (Scalar gamma_)                 {gamma = gamma_;}
  void setTemporalWeight(Scalar ww)             {temporalWeight = ww;};

  void setTempJump  (Scalar tempJump_)          {tempJump  = tempJump_;}
  void setRotJump   (Scalar rotJump_)           {rotJump   = rotJump_;}
  void setDepthJump (Scalar depthJump_)         {depthJump = depthJump_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}
//  void setP2d (Scalar* p2d_) { p2dSet = p2d_;}
  void setCenters (Scalar* centers_) { centers = centers_;}

  const std::vector<int>& getIdl()  { return Idl; }
  const std::vector<int>& getIdk()  { return Idk; }

  const std::vector<Scalar>& getF11()  { return F11; }
  const std::vector<Scalar>& getF00()  { return F00; }
  const std::vector<Scalar>& getF10()  { return F10; }
  const std::vector<Scalar>& getF01()  { return F01; }

  const std::vector<Scalar>& getF0()  { return F0; }
  const std::vector<Scalar>& getF1()  { return F1; }

  const std::vector<P3>& getTra()      { return tra; }
  const std::vector<M3>& getRot()      { return rot; }
  const std::vector<P3>& getCenters()  { return centers; }

  void prepare( int Npixel, int nSegments_);//, const mxArray* edges_ )//, Scalar* centers_

  size_t edges_num() {return Idk.size();};

  /// returns the map describing the submodularity of the edge between a segment and its neighbour
  std::vector< std::map< int, Scalar> >&  getWeightMap() { return weightMap; };

  void setSegImg    ( int _w, int _h, int* _segImg, int _nSegments );
  /*!
  compute the energy - debug info no need to fill any containers etc. just return energy
  */
  Scalar compute_score_combiDepth_consistentScore( 
    const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt );


// that is not the problem - rather that the graph constructed becomes too large
#define _use_safe_version_to_test_against_
#ifdef _use_safe_version_to_test_against_
  ///////////////////// new version without local rotation stuff ////////////////////
  //compute_score_2D_homs
  /// PURE 2d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution
  Scalar compute_score_combiDepth_consistent(const int id, const P4i& box, const std::vector<int>& loc2glob,// const int localSize, 
    const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt, //)//, std::vector<int>& seg2Prop )
    std::vector<Scalar>& F00, std::vector<Scalar>& F01, std::vector<Scalar>& F10, std::vector<Scalar>& F11,
    std::vector<Scalar>& F0, std::vector<Scalar>& F1, std::vector<int>& Idk, std::vector<int>& Idl );

//  void precompute_Displacements();
#else
  void precompute_Displacements( const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt );

  Scalar compute_score_combiDepth_consistent(const int id, P4i box, std::vector<int>& loc2glob, int localSize, 
    const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt, std::vector<int>& seg2Prop );
#endif

private:

  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2dDiff( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const; 

  inline void get2dDiff_x2disp( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                                const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                                P4& dispDiff, P2& flowDiff ) const;

  inline void get2dDiff_x1disp( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                                const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                                P4& dispDiff, P2& flowDiff ) const;

  inline void get2dProjections( const P3& p_t0, P2& pixP_l, P2& pixP_r) const;

  inline void get2dDiff_( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                          const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
    P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const ;


  //////////////////////////////////////////////////

  M3 Pl;//R=I t=0
  M3 Pr;//K*R
  P3 pr;//K*t

  /// the segment ids for querying the assigned segment of a pixel
  int* segImg;

  int nSegments;

  /// the width of the image region
  int width;
  /// the height of the image region
  int height;

  /// the edges 0,1,2 are 2d coordinates of point p 
  std::vector<P3> edgesP;
  /// the edges 0,1,2 are 2d coordinates of point q
  std::vector<P3> edgesQ;

  /// the centers of the pathces in 2d
//  std::vector<P3> centers;
  M3 iK;

  /// transfering per pixel normals (image base) into 3d normals (camear space)
  M3 iiKt;

  std::vector<P3> centersA;
  std::vector<P3> centersB;
  Scalar* centers;

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

#ifdef __use_m128__
  std::vector<__m128d> homs_l_sse;  // homographies
  std::vector<__m128d> homs_r_sse;  // homographies
  std::vector<__m128d> homs_rt_sse;  // homographies
  std::vector<__m128d> homs_li_sse;  // homographies
  std::vector<__m128d> homs_rti_sse;  // homographies
  std::vector<__m128d> homs_ti_sse;  // homographies, special case backward motion inversely applied on current frame
#endif

  std::vector<M3> homs_l;  // homographies
  std::vector<M3> homs_r;  // homographies
  std::vector<M3> homs_rt; // homographies
  std::vector<M3> homs_li;  // homographies
  std::vector<M3> homs_rti; // homographies
  std::vector<M3> homs_ti; // homographies

  /// the 2d point coordinates for each pixel
  //  Scalar *p2dSet;
  
  std::vector<Scalar> F00;
  std::vector<Scalar> F01;
  std::vector<Scalar> F10;
  std::vector<Scalar> F11;
  std::vector<Scalar> F0;
  std::vector<Scalar> F1;
  
  Scalar edgeScale;
  Scalar gamma;
  Scalar temporalWeight;
  double epsilon;

  Scalar tempJump;
  Scalar rotJump;
  Scalar depthJump;
  Scalar rotWeight;

  /// keeps track of the weights at the edges, test: use for non-submodular edge selection for occlusion mapping
  std::vector< std::map< int, Scalar> > weightMap;

  /// weights on the edges per pixel (4/8-neighbourhood)
  Scalar* halfEdgeX;
  Scalar* halfEdgeY;

  /// weight on the cross edges (8-neighbourhood)
  Scalar* halfEdgeXY;
  Scalar* halfEdgeiXY;

  /// used to adjust the threshold in case of 2d jumps of more than e.g. 1 pixel (given as parameter)
//  Scalar epix;

  /// cache for 2d-smoothness computation per pixel with 
  std::vector<P3> homoDisplace_l;
  std::vector<P3> homoDisplace_r;
  std::vector<P3> homoDisplace_rt;
};

#ifdef _use_safe_version_to_test_against_
#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANORPIXEL__cpp)
#include "EvalEnergyFull2FramePixel.cpp"
#endif
#else
#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANORPIXEL_UNSAFE__cpp)
#include "EvalEnergyFull2FramePixel_obsolete.cpp"
#endif
#endif
#endif