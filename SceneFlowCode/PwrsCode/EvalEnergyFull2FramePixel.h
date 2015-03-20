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

///////////////////////////////////////////////////////
////////// Smoothness term per pixel        ///////////
///////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANORPIXEL__
#define __ENERGY_ROTTRANORPIXEL__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "DataDefinitions.h"
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

  typedef Math::Mat3x3T<Scalar>    M3; 

  EvalEnergyFullFramePixel(int _width, int _height)
    : gamma(1.0), epsilon(0.00), normals(NULL), tra(), rot(), rotJump(15.),
    segImg(NULL), depthJump(sqrt(1300.)), nSegments(0), weightMap(), 
    width(_width), height(_height), centers(NULL), iK(1,0,0,0,1,0,0,0,1), 
    halfEdgeX(NULL), halfEdgeY(NULL), halfEdgeXY(NULL), halfEdgeiXY(NULL)
  {};

  ~EvalEnergyFullFramePixel(){};

  void setHom( Scalar* homs_, int nSegments )
  {
    homs_r.clear();  homs_r.resize(nSegments);
    homs_l.clear();  homs_l.resize(nSegments);
    homs_rt.clear(); homs_rt.resize(nSegments);

    for(int i =0; i< nSegments; i++)
    {
      homs_r  [i] = M3(&homs_[i*9              ]) * iK;
      homs_l  [i] = M3(&homs_[i*9+nSegments*9  ]) * iK;
      homs_rt [i] = M3(&homs_[i*9+nSegments*9*2]) * iK;
    }
  }

  void setSegImg    ( int _w, int _h, int* _segImg, int _nSegments )
  { 
    height=_h;
    width =_w;
    nSegments = _nSegments;
    segImg    = _segImg;
  };

  void set2dMotionMatrix ( Scalar* Kl, Scalar* Kr, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(Kl);
    Pr = M3(Kr) * M3(Rot);
    pr = M3(Kr) * P3(mc);
  };

  void setRotWeight (Scalar rotWeight_)         {rotWeight = rotWeight_;};
  void setGamma (Scalar gamma_)                 {gamma = gamma_;}

  void setRotJump   (Scalar rotJump_)           {rotJump = rotJump_;}
  void setDepthJump (Scalar depthJump_)         {depthJump = depthJump_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  void setK (Scalar* K_)
  {
    assert (homs_r.size() < 1); // this is to be called before setHoms
    iK=M3(K_);
    iK.invert();
  }

  void setCenters (Scalar* centers_) { centers = centers_;}

    /// weights on the edges, based on the input image
  void setHalfEdges( Scalar* _halfEdgeX, Scalar* _halfEdgeY) {halfEdgeX = _halfEdgeX; halfEdgeY = _halfEdgeY;}
  /// weights on the cross egges in a 8 neighbourhood
  void setCrossHalfEdges( Scalar* _halfEdgeXY, Scalar* _halfEdgeiXY) {halfEdgeXY = _halfEdgeXY; halfEdgeiXY = _halfEdgeiXY;}

  const std::vector<int>& getIdl()  { return Idl; }
  const std::vector<int>& getIdk()  { return Idk; }

  const std::vector<Scalar>& getF11()  { return F11; }
  const std::vector<Scalar>& getF00()  { return F00; }
  const std::vector<Scalar>& getF10()  { return F10; }
  const std::vector<Scalar>& getF01()  { return F01; }

  const std::vector<Scalar>& getF0()  { return F0; }
  const std::vector<Scalar>& getF1()  { return F1; }

  void prepare( int Npixel, int nSegments_);

  size_t edges_num() {return Idk.size();};

/// PURE 3d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution
#ifdef __USE3D__
  Scalar compute_score_combiDepth( int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize );
#else
  Scalar compute_score_3D( int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize );
#endif

  //compute_score_2D_homs
  /// PURE 2d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution
#ifndef __USE3D__
  Scalar compute_score_combiDepth(const int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize );
#else
  Scalar compute_score_2D(const int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize );
#endif

private:
  
  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2dDiff( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const 
  {
    P3 pixP_t0_l = Pl * p_t0;pixP_t0_l /= pixP_t0_l[2];
    P3 pixQ_t0_l = Pl * q_t0;pixQ_t0_l /= pixQ_t0_l[2];

    P3 pixP_t1_l = Pl * p_t1;pixP_t1_l /= pixP_t1_l[2];
    P3 pixQ_t1_l = Pl * q_t1;pixQ_t1_l /= pixQ_t1_l[2];

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
  
  std::vector<M3> homs_l;  // homographies
  std::vector<M3> homs_r;  // homographies
  std::vector<M3> homs_rt; // homographies
  
  std::vector<Scalar> F00;
  std::vector<Scalar> F01;
  std::vector<Scalar> F10;
  std::vector<Scalar> F11;
  std::vector<Scalar> F0;
  std::vector<Scalar> F1;
  
  Scalar gamma;
  double epsilon;
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
};

#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANORPIXEL__CPP)
#include "EvalEnergyFull2FramePixel.cpp"
#endif

#endif