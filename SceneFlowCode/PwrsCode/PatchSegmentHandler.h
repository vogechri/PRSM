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

//////////////////////////////////////////////////////////////////////////////
// Manage the partial instantion of the pixel graph in the pixel refinement //
//////////////////////////////////////////////////////////////////////////////

#ifndef __PatchSegmentHandler__h
#define __PatchSegmentHandler__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <algorithm>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Math;

#define __USE_FIXED_BORDERS__

/* 
1st set some centers, defining the positions of the patches.
needs random order of segments/centers
1st triggers :generateEdgeImage, init stuff
 ------------------------
2nd assign each pixel in the patch a number/id and 
mappings glob2loc, loc2glob, vector of global ids for the data generation
function : initMappings( int patch )
other functions:
---------
getglobIds() returns a vector of global ids in the patch
int getLocId( int globalID ) returns the local id for a global one
int getCurrentPatch(4x int ): returns the current segid and the 4 corners of the patch

not used: use precomputed input instead for smoothness -> more flexible
Scalar getEdgeWeight (int locId1,int locId2) returns the edge weight of these local ids
------------------------
*/ 
///////////////////////////////////////


template<typename Scalar> class PatchIdGenerator
{
public:

  const int static dx[4]; //= {-1,0,1,0};
  const int static dy[4];// = {0,1,0,-1};

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
  typedef Math::VectorT<int, 2>     P2i;

  typedef std::pair <int,int>       Seg_ij;

  /// segImg must be flipped (so c++ coordinates, not matlab ones)
  PatchIdGenerator(int _w, int _h, int _patchSize, Scalar* _centers, int* _initSegImg,
                    Scalar* _Kl, int _nInitSegments,  int* _segImg, Scalar* _refimg, bool doCensus_ = false, int cSize_ = 1 )
    : patchSize(_patchSize), Kl(_Kl), nInitSegments(_nInitSegments), currentPatch(-1), mrfIds(0), cRad(cSize_),
    currentSeg(-1), segImg(_segImg), w(_w), h(_h), startX(0), endX(0), startY(0), endY(0), doCensus(doCensus_)
  {
    setCensus(cSize_);

    setCenters(_centers);
    loc2glob.reserve( (2*patchSize+1)*(2*patchSize+1) );
    glob2loc.resize(w*h);
    edgeMap.reserve (w*h*dxyC.size());
    locIds.reserve ( (2*patchSize+1)*(2*patchSize+1) );
    globIds.reserve( (2*patchSize+1)*(2*patchSize+1) );
    setInitSegImg ( _initSegImg );
    setSegImg ( _segImg );
  }

  ~PatchIdGenerator(){};

  /// returns id of the center pixel in the patch, and the coordinates of the patch
  int getCurrentPatch( int& _startX, int& _endX, int& _startY, int& _endY )
  {
    _startX = startX;_endX = endX;
    _startY = startY;_endY = endY;
    return currentSeg;
  }

  /// defines together with the centers the proposal patches
  void setInitSegImg( int* _segImg )
  { 
    startSegs.clear();
    startSegs.resize(nInitSegments);

    for(int patch =0;patch<nInitSegments;patch++)
      startSegs[patch] = _segImg[ (centers[patch])[0]+(centers[patch])[1]*w ];
  }

  /// put expansion centers could on a grid - could also be appended to original list
  void gridCenters( int gridSize=0 )
  {
    /// everything else would just take too much time
    if ( gridSize < 5 ) gridSize = max( patchSize/2, 3);

    centers.clear();
    
    for (int dy=gridSize/2; dy<h ;dy += gridSize)
      for (int dx=gridSize/2; dx<w ;dx += gridSize)
        centers.push_back( P2i( dx, dy ) );

    nInitSegments = centers.size();

    printf("Currently : %d centers\n", nInitSegments);

    if (segImg)
      setInitSegImg( segImg  );
  }

  /// return the local id map, needed for occlusion handling and data computation (occMapBuffer computes/stores the 2d positions of the warped pixel)
  const std::vector<int>& getLocIds() {return locIds;}

  /// map a global id to a local one
  inline int getLocId( int globalID ) const 
  {
    int qx = globalID % w;
    int qy = globalID / w;

    if ( qx<startX || qy<startY || qx>=endX || qy>=endY )
      return -1;
    else
      return qx-startX+(qy-startY)*(endX-startX);
//      qx-startX+qy*(endX-startX);
  }

  inline int getLocIdFromGlobal( int globalID ) const
  {
    return glob2loc[globalID];
  }

  /// this is the mrf id!
  inline int getMrfId( int globalID ) const
  {
    int qx = globalID % w; 
    int qy = globalID / w;

    if ( qx<startX || qy<startY || qx>=endX || qy>=endY ) 
      return -1;
    else
      return locIds[qx-startX+(qy-startY)*(endX-startX)];
  }

  /// used for data computation globIds stores the global pixel positions of the pixel considered for data term evaluation
  const std::vector<int>& getGlobIds() {return globIds;}

  const std::vector<int>& getAllGlobIds() {return allGlobIds;}

  const std::vector<bool>& getFixedIds() {return freeFixedIds;}

  /// ids of the neighbours in the image space
  const std::vector<int>& getEdgeIds() {return edgeIds;}

  /// ids of the lower right neighbours in the image space
  const std::vector<int>& getOtherEdgeIds() {return otherEdgeIds;}

  int getNMrfVars(){return mrfIds;};

  int getNVars(){return globIds.size();};

  int getNProposals() {return nInitSegments;};

  /// prepare energy computation for a patch, new globIds for all! pixel in region, locids only for inner and non segment pixel used in smoothing
  void initMappingsNew( int patch, bool doEdgeMapping );
  
  /// return the number of new proposals (old ones and new ones) potentially increase expansion area instead? - also keep original positions(segment centers)
  int refreshMappings( bool onlyNew=true );

private:

  void setCensus(int cSize);

  /// edgeIds are only used for data term, so all edges in window are concerned here so getLocId not getMrfId
  void initEdgeMappingNew( );

  /// defines the current situation from which the solution is expanded - do not confuse with setInitSegImg
  void setSegImg( int* _segImg )
  { 
    segImg    = _segImg;
  }

  /// compute the 2d centers of the patches
  void setCenters (Scalar* centers_) 
  { 
    centers.clear();
    centers.resize(nInitSegments);
    for (int i=0; i<nInitSegments ;i++)

    {
      P3 cc = Kl * P3(&centers_[3*i]); cc = cc / cc[2];
      centers[i] = P2i( floor(cc[0]-0.5),  floor(cc[1]-0.5) );
    }
  }

  ///////////////////////////////////////////////

  /// pointer to the 2d centers
  std::vector<P2i> centers;

  M3 Kl;

  /// the radius of the census transform
  int cRad;

  /// the size of the pathces, that amount of pixel can get the id of its assigned segment
  int patchSize;

  /// the amount of free variables
  int mrfIds;

  /// the weights of the edges based on intensity difference
  std::vector<Scalar> edgeMap;
  
  /// map from local to global system
  std::vector<int> loc2glob;
  /// map from global to local system, could also be computed completely: qx=(px-startx), .. id=qx+qy*(endX-startX)
  std::vector<int> glob2loc;

  /// local ids of the pixel in the patch
  std::vector<int> locIds;

  /// global ids of the pixel in the patch, which are free
  std::vector<int> globIds;

  /// stores all the global ids of the local  region for data term only, later data term has to be trimmed to actual mrf ids
  std::vector<int> allGlobIds;

  bool doCensus;

  /// edge map: shows which edges are adjacent to other pixel: -1 if none pattern is upper left to upper right, left, so 4 edges per pixel
  std::vector<int> edgeIds;

  /// edge map: same as above but with the other 4 neighbours
  std::vector<int> otherEdgeIds;

  /// states whether a pixel included in the patch (==has a local id) should still be fixed to the current segment by adding a huge penalty if it changes
  std::vector<bool> fixedIds;

  /// states whether a pixel included in the optimization (==has a local id) should still be fixed to the current segment by adding a huge penalty if it changes
  std::vector<bool> freeFixedIds;

  int startX; 
  int endX;
  int startY;
  int endY;

  /// the initial assignments to the centers -> use these for expansion (or more)
  std::vector<int> startSegs;

  /// the current segment picked for expansion
  int currentSeg;

  /// the current patch picked for expansion
  int currentPatch;

  /// the image storing the segmentation
  int* segImg;

  /// width
  int w;
  /// height
  int h;

  /// the amount of segmetns in the images
  int nInitSegments;

  /// the neighbourhood in the censuss term
  std::vector<int> dxC;
  std::vector<int> dyC;
  std::vector<int> dxyC;
};


template<typename Scalar>
struct binaryDataScoresBuffer
{
  typedef std::pair <int,int>  PixSegPair;
  binaryDataScoresBuffer( int _w, int _h, std::vector< char >& _initDataScore, int cSize_ = 1 ) : w(_w), h(_h)
  {
    initDataScore( _initDataScore );
  };

  ~binaryDataScoresBuffer() {};

  /// initialize the data score 
  void initDataScore( std::vector<char>& _dataScore )
  {
    assert (_dataScore.size() == 2*dxyC.size()*w*h);
    dataScore.resize( _dataScore.size()/2 );
    for (int i=0;i<dataScore.size();i +=2)
      dataScore[i] = _dataScore[2*i];
  }

  ////////////////////////////////////////////
  ///width
  int w;
  ///height
  int h;
  /// the data score of the current solution, without considering occlusions
  std::vector< char > dataScore;
};


template<typename Scalar>
struct unaryDataScoresBuffer
{
  typedef std::pair <int,int>  PixSegPair;
  unaryDataScoresBuffer( int _w, int _h, std::vector< Scalar >& _initDataScore, int cSize_ = 1 ) : w(_w), h(_h)
  {
    initDataScore( _initDataScore );
  };

  ~unaryDataScoresBuffer() {};

  void initDataScore( std::vector<Scalar>& _dataScore ) // penalty score attached
  {
    assert (_dataScore.size() == 2*w*h);
    dataScore.resize( _dataScore.size() );
    for (int i=0;i<dataScore.size();i +=1)
    {
      dataScore[i]   = _dataScore[  i];
    }
  }

  /// writes the penalty scores stored in dataScore (in the buffer) to the vector to be as data score for a window to be optimized
  /// combines external score from afd with internally buffered score
  void constructDataScore( std::vector<Scalar>& newScores, const std::vector<Scalar>& trialScores, std::vector<int>& loc2Glob, PatchIdGenerator<Scalar>& genPatchIds );

  //////////////////////////////////////
  /// update the data score, based on a score map and the changes: 0 nothing changed 1 changed segment map
  void updateDataScore( const std::vector<Scalar>& trialScores, const std::vector<int>& loc2Glob, 
                        const PatchIdGenerator<Scalar>& genPatchIds, const std::vector<int>& changes );

 ////////////////////////////////////////////

  ///width
  int w;
  ///height
  int h;
  /// the data score of the current solution, without considering occlusions, penalties are attched to the end
  std::vector< Scalar > dataScore;
};

#if defined(INCLUDE_TEMPLATES) && !defined(__PatchSegmentHandler__cpp)
#include "PatchSegmentHandler.cpp"
#endif

#endif // __PatchSegmentHandler__h
