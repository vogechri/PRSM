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

#ifndef __OcclusionMappingBufferVC__h
#define __OcclusionMappingBufferVC__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <vector>
#include "mex.h" // include after VectorT and Mat3x3T

// use iccv class
#include "SegmentGrowing_HeapHelper.h"
#include "Templates/HeapT.h"

#include "genHom.h"

/// an occlusion is only given if the depth deviation is more or less than 1 pct of the depth at the pixel
//#define _OccPct_ 0.01

typedef std::pair <int,int>  PixSegPair;

/// the list of the occPixelList's has to be sorted by the ownIds to be able to merge lists
template<typename Scalar>
struct occPixelListVC
{
  typedef std::pair <int,int>  PixSegPair;
  occPixelListVC( ) : ownId(-1), segId(-1), occluders(), fullOccluder(-1)
  {};

  occPixelListVC( int _pid, int _segid, std::vector < PixSegPair > _occluders, int _fullOccluder )
    : ownId(_pid), segId(_segid), occluders(_occluders), fullOccluder(_fullOccluder), occludes(-1)
  {};

  ~occPixelListVC() {};
  ////////////////////////////////////////////

  /// stored are pairs: pixel id, assign a 0 or 1, interpretation 0: occludes if id=0, encode x_i in poly: if x_i=1 not occluded
  std::vector < PixSegPair > occluders;

  int ownId;        // id of the occluded pixel GLOBAL ? NO!! - LOCAL !

  int occludes;

  int segId;        // seg assignment to pixel

  /// the pixel is always occluded - either from outside or from a fre pixel - no matter what is assigned ot the pixel
  int fullOccluder; // needed to encode a general occlusion (pixel id not really needed)
  // replace data term by occ pen or add unary poly with (occ-data)
};

template<typename Scalar>
struct oMapBufferStorageVC
{
public:

  oMapBufferStorageVC() : pixPos(-1), depth(-1.0), segId(-1), px(-1), py(-1) {};

  // with segid interesting for the lookup of new id assignments (a 1 to a pixel)
  oMapBufferStorageVC( int _pixPos, int _px, int _py, Scalar _depth, int _segId )
    : pixPos(_pixPos), px(_px), py(_py), depth(_depth), segId(_segId)
  {};

  ~oMapBufferStorageVC() {};
  ///////////////////////////////////////
  int px;
  int py;
  int pixPos;
  Scalar depth;
  int segId;
};

template<typename Scalar>
struct oMapBufferLookup
{
public:

  oMapBufferLookup( ) : vecPos(-1), pixPos(-1) {};// oob
  oMapBufferLookup( int _pixPos, int _vecPos) :  pixPos(_pixPos), vecPos(_vecPos)
  {};

  ~oMapBufferLookup() {};
  ///////////////////////////////////////
  int pixPos;
  int vecPos;
};
//////////////////////////////////////////////////////////////////////////

/*! Projects one moving plane segmentation into a 2nd image.
Exploits the moving plane information to project the moving plane segments.
Fills potential holes via seed growing with a priority queue.
The class is used to deliver an initial solution for other views of the scene.
*/
template<typename Scalar> class OcclusionMapBufferVC
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
  typedef std::pair <int,int>       Seg_ij;

  OcclusionMapBufferVC(int _maxPatchSize, bool doOcclusionStuff = true)
    : segImg(NULL), w(0), h(0), occPenalty(0), p2d(NULL), newIteration(0), 
    startX(0), endX(0), startY(0), endY(0), maxPatchSize(_maxPatchSize), locIds(NULL), genOccData(doOcclusionStuff),
    gHom2nd(NULL), gHom3rd(NULL)
  {}

  OcclusionMapBufferVC(int _maxPatchSize, int w_, int h_, bool doOcclusionStuff = true)
    : segImg(NULL), w(w_), h(h_), occPenalty(0), p2d(NULL), newIteration(0), 
    startX(0), endX(0), startY(0), endY(0), maxPatchSize(_maxPatchSize), locIds(NULL), genOccData(doOcclusionStuff),
    gHom2nd(NULL), gHom3rd(NULL)
  { }

  ~OcclusionMapBufferVC(){};

  void setP2d (Scalar* p2d_) { p2d = p2d_;}

  void setSegImg( int _w, int _h, int* _segImg, int _nSegments )
  { 
    h=_h;
    w=_w;
    nSegments = _nSegments;
    segImg    = _segImg;

    store0.clear();
    store0.resize(w*h);
    store1.clear();
    store1.resize(w*h);
    for (int i=0;i<w*h;i++)
    {
      store0[i].reserve(10);
      store1[i].reserve(10);
    }

    lookup0.clear();lookup0.resize(w*h);
    lookup1.clear();lookup1.resize(w*h);
  };

  void SetOccPenalty (Scalar _occPenalty) {occPenalty = _occPenalty;};

  void setSecondTestView( genHomoG<Scalar>* gHom, const std::vector<int>& solutionB )
  {
    gHom2nd = gHom;
    solution2nd.resize( solutionB.size() );
    for (int i=0;i< solution2nd.size();i++)
      solution2nd[i] = solutionB[i];
  }


  void setThirdTestView( genHomoG<Scalar>* gHom, const std::vector<int>& solutionB )
  {
    gHom3rd = gHom;
    solution3rd.resize( solutionB.size() );
    for (int i=0;i< solution3rd.size();i++)
      solution3rd[i] = solutionB[i];
  }

  /// the problem is that i have no data here; also occ and mvp are dependent on data, also true opt is my problem
  void extendSegmentation( genHomoG<Scalar>* gHom, const std::vector<P3>* normals, const mxArray* edges2_, 
                           Scalar *centers_, int nSegments, Scalar *centers2_, const std::vector<int>& solution, 
                           int* segImg, int* segImg2, std::vector<int>& solution2, int nSegments2, Scalar* Kl_, Scalar* Kl2_  );

private:

  Scalar getPenalty ( int ownSegment, int proposal, const M3& iHom_0, const P3& centre, 
                      const std::vector<int>& solution, const std::vector<P3>* normals, 
                      Scalar *centers_, int& seg1);

  ///////////////////////////////////////////////

  /// from the pixel ina the reference image to the pixel it is projected onto
  std::vector<oMapBufferLookup<Scalar> > lookup0;
  std::vector<oMapBufferLookup<Scalar> > lookup1;

  /// information about the pixel in the data image: depth for occlusions and id of the original pixel
  std::vector< std::vector< oMapBufferStorageVC<Scalar> > > store0;
  /// same as above but instead of the id also a iteration number 
  std::vector< std::vector< oMapBufferStorageVC<Scalar> > > store1;

  int startX; 
  int endX; 
  int startY;
  int endY;

  /// helper variable to identify correct occlusions in the current setting
  int newIteration;

  /// width
  int w;
  /// height
  int h;

  /// the amount of segmetns in the images
  int nSegments;

  /// a vector mapping the ids in the patch to local ids for the energy 
  std::vector<int>* locIds;

  Scalar *p2d;

  /// the segment ids for later test w.r.t. occlusions
  int* segImg;

  Scalar occPenalty;

  int maxPatchSize;

  bool genOccData;

  // used in projecting solutions from different views:
  genHomoG<Scalar>* gHom2nd;
  std::vector<int> solution2nd;
  genHomoG<Scalar>* gHom3rd;
  std::vector<int> solution3rd;
};

#if defined(INCLUDE_TEMPLATES) && !defined(__OMapBufVC__cpp)
#include "OcclusionMappingBufferVC.cpp"
#endif

//#undef _OccPct_
#endif // __OcclusionMapping__h
