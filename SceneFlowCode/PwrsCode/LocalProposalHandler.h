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

////////////////////////////////////////////////////////////////
// Group segments to cover the whole image//
////////////////////////////////////////////////////////////////

#ifndef __LocalProposalHandler__h
#define __LocalProposalHandler__h

#include <limits>
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "Templates/HeapT.h"
#include "Segment_HeapHelper.h"

#include <list>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

using namespace std;
using namespace Math;

#define __USE_FIXED_BORDERS__

/*! 
  Generate local proposal sets build from region growing proposals spaced
  a certain amount away. Used is a grid with a given center distance.
  Segments lying on the grid centers are used for expansion.
  Growing is based on distance of the patch centres in 2d.

  As a results we get a massive reduction in full image proposals for the optimization
  leading to less data term and qpbo evaluations needed.

  Goal is to allow every segment to be expanded at least once.
  Therefore we keep track about how often a segment got already selected.

  We center the grid at a segment not yet selected. 
  We expand the segments using a priority queue with the minimal 
  distance of the patch centre to the grid as criterium.

  Potentially one could select the closesest segments to the grid centers, 
  or only close ones not yet selected as expanded proposal.
  Leading to fewer full proposal sets.
*/
///////////////////////////////////////

template<typename Scalar> class LocalProposalHandler
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
  typedef Math::VectorT<int, 2>     P2i;
  typedef Math::VectorT<int, 4>     P4i;

  typedef HeapStoreEntryT<double>  HeapStoreEntry;
  typedef vector< HeapStoreEntry > HeapEntryStore;

  /// segImg must be flipped (so c++ coordinates, not matlab ones)
  LocalProposalHandler(int _w, int _h, int _gridSize, Scalar* _centers,
                       Scalar* _Kl, int _nSegments,  int* _segImg)
    : gridSize(_gridSize), Kl(_Kl), nSegments(_nSegments), segImg(_segImg), w(_w), h(_h), 
      _heapEntryStore(), _heapInterface(&_heapEntryStore), _heap(_heapInterface), nGrids(0)
  {
    setCenters(_centers);
    setSegImg ( _segImg );
    pickedSegments.clear();
    pickedSegments.resize(nSegments,0);
    generateSmallBBox();
  }

  ~LocalProposalHandler(){};

  void setGridSize(int gSize)
  {
    gridSize = gSize;
  }

  /// return the compressed proposals
  std::vector<int> getProposalMap() {return proposalMaps;};

  /// input defines neighborhood relationship
  void generateProposals(const mxArray* edges_)
  {
    nGrids = 0;
    pickedSegments.clear();
    pickedSegments.resize(nSegments, 0);

    Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, 1) );
    int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, 1) )/5;
    for (int j =0; j < nIds ;j++)
    {
      int nId = (int) (edge [5*j ])-1;
      int test = nId;
    }

    proposalMaps.clear();
    int elements = nSegments * nSegments  / (w*h) * (gridSize*gridSize);
    proposalMaps.reserve( elements );
    for (int i = 0;i<nSegments; i++)
    {
      if (pickedSegments[i] > 0) continue; // segment already considered

      nGrids++;

      int startPos = proposalMaps.size();
      proposalMaps.resize( nSegments + proposalMaps.size() );

      _heap.clear();
      _heapEntryStore.clear();
      // cost w*h (larger than maximal) , startSeg : -1 == not touched yet, -2: fixed already 
      _heapEntryStore.resize( nSegments, HeapStoreEntry(w*h, -1) ); 

      // init the grid at segments i
      std::vector<int> startSegs;
      initGridJitter(i, startSegs, edges_, pickedSegments);

      // fix segments first
      for (int k=0; k< startSegs.size(); k++)
      {
        int cSeg = startSegs[k];
        proposalMaps[startPos + cSeg ] = cSeg;
        _heapEntryStore[cSeg].setStartSeg(-2); 
        pickedSegments[cSeg]++;
      }

      // put all neighs of the grid segments into the heap:
      for (int k=0; k< startSegs.size(); k++)
      {
        int sSeg = startSegs[k];
        // put neighs of cSeg into heap:

        Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, sSeg) );
        int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, sSeg) )/5;
        for (int j =0; j < nIds ;j++)
        {
          int nId = (int) (edge [5*j ])-1;

          HeapStoreEntry& hse = _heapEntryStore[nId];
          int bestSeg = hse.startSeg();

          // else already fixed
          if( bestSeg > -2) // replace if distance is smaller
          {
            Scalar cost = (centers[sSeg]-centers[nId]).sqrnorm();

            int heapid = hse.heapId();
            if (heapid<0 && bestSeg>0)
              printf("Weird\n");

            // if so enter or update heap:
            if (cost < hse.cost())
            {
              hse.setCost(cost);
              hse.setStartSeg(sSeg);
              if ( bestSeg == -1)
                _heap.insert( HeapEntryD (nId) );
              else // already exists 
                _heap.update( HeapEntryD( nId) );
            }
          }

        } // over neighs
      } // over grid Segments
      int numRed(0);
      while ( !_heap.empty() )
      {
        numRed++;

        HeapEntryD id = _heap.front();
        _heap.pop_front();
        int cSeg = id.pId();
        HeapStoreEntry& heapE = _heapEntryStore[ cSeg ];

        // fix segment , add neighbors if smaller
        int sSeg = heapE.startSeg();
        proposalMaps[startPos + cSeg ] = sSeg;

        // invalidate segments once and for all
        heapE.setStartSeg ( -2 );

        // loop neighs:
        Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, cSeg) );
        int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, cSeg) )/5;
        for (int j =0; j < nIds ;j++)
        {
          int nId = (int) (edge [5*j ])-1;

          HeapStoreEntry& hse = _heapEntryStore[nId];
          int bestSeg = hse.startSeg();

          // else already fixed
          if( bestSeg > -2) // replace if distance is smaller
          {
            Scalar cost = (centers[sSeg]-centers[nId]).sqrnorm();

            int heapid = hse.heapId();

            if (heapid<0 && bestSeg>0)
              printf("Weird\n");

            // if so enter or update heap:
            if (cost < hse.cost())
            {
              hse.setCost(cost);
              hse.setStartSeg(sSeg);
              if ( bestSeg == -1 )
                _heap.insert( HeapEntryD (nId) );
              else // already exists 
                _heap.update( HeapEntryD( nId ) );
//                _heap.update( HeapEntryD( hse.heapId() ) );
            }
          }
          // doing some stupid stuff
          int abs = 0;

        } // over neighs
      } // while heap not empty
    } // all segments are in a proposal map

    generateLargeBBox();
  }

  int getNGrids() {return nGrids;}

  /// returns how much on average a segment was covered
  Scalar efficiency()
  {
    int sumAll(0); 
    for (int i=0; i<pickedSegments.size(); i++)
      sumAll += pickedSegments[i];
    return Scalar(sumAll)/Scalar(pickedSegments.size());
  }

  /// input defines neighborhood relationship
  void generateProposalsSample(const mxArray* edges_, unsigned int trials = 20)
  {
    nGrids = 0;
    pickedSegments.clear();
    pickedSegments.resize(nSegments, 0);

    proposalMaps.clear();
    int elements = nSegments * nSegments  / (w*h) * (gridSize*gridSize);
    proposalMaps.reserve( elements );

    /// here the first n contain the free segments
    std::vector<int> freeSegments(nSegments, 0);
    std::vector<int> swappedSegs (nSegments, 0); // stores the position of segment i in freeSegments
    std::vector<int> sampleSegments(nSegments, 0);
    for (int i=0;i<nSegments;i++)
    {
      swappedSegs[i] = i;
      freeSegments[i] = i;
    }
    int nFreeSegments = nSegments-1;

    while( freeSegments.size() > 0)
    {
      nGrids++;
      std::copy( freeSegments.begin(), freeSegments.end(), sampleSegments.begin() );
      std::random_shuffle( sampleSegments.begin(), sampleSegments.end() );

      int maxHits = countGridJitter(sampleSegments[0], edges_, pickedSegments );
      
      int bestHit = sampleSegments[0];
      for (int trial=1;trial<std::min(trials, (unsigned int) (sampleSegments.size()));trial++)
      {
        int hits = countGridJitter(sampleSegments[trial], edges_, pickedSegments );
        if (hits > maxHits)
        {
          maxHits = hits;
          bestHit = sampleSegments[trial];
        }
      }

      int startPos = proposalMaps.size();
      proposalMaps.resize( nSegments + proposalMaps.size() );

      _heap.clear();
      _heapEntryStore.clear();
      // cost w*h (larger than maximal) , startSeg : -1 == not touched yet, -2: fixed already 
      _heapEntryStore.resize( nSegments, HeapStoreEntry(w*h, -1) ); 

      // init the grid at segments i
      std::vector<int> startSegs;
      initGridJitter(bestHit, startSegs, edges_, pickedSegments);

      // fix segments first
      for (int k=0; k< startSegs.size(); k++)
      {
        int cSeg = startSegs[k];
        proposalMaps[startPos + cSeg ] = cSeg;
        _heapEntryStore[cSeg].setStartSeg(-2); 

        // tausch: free position in freeSegments now filled with last free segment
        // 
        if (pickedSegments[cSeg] == 0)
        {
          int lastFree = freeSegments[ nFreeSegments-- ];
          freeSegments[ swappedSegs[cSeg] ] = lastFree;
          swappedSegs[ lastFree ] = swappedSegs[cSeg];
        }
        pickedSegments[cSeg]++;
      }
      freeSegments.resize( nFreeSegments+1 );
      sampleSegments.resize( nFreeSegments+1 );

      // put all neighs of the grid segments into the heap:
      for (int k=0; k< startSegs.size(); k++)
      {
        int sSeg = startSegs[k];
        // put neighs of cSeg into heap:

        Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, sSeg) );
        int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, sSeg) )/5;
        for (int j =0; j < nIds ;j++)
        {
          int nId = (int) (edge [5*j ])-1;

          HeapStoreEntry& hse = _heapEntryStore[nId];
          int bestSeg = hse.startSeg();

          // else already fixed
          if( bestSeg > -2) // replace if distance is smaller
          {
            Scalar cost = (centers[sSeg]-centers[nId]).sqrnorm();

            int heapid = hse.heapId();
            if (heapid<0 && bestSeg>0)
              printf("Weird\n");

            // if so enter or update heap:
            if (cost < hse.cost())
            {
              hse.setCost(cost);
              hse.setStartSeg(sSeg);
              if ( bestSeg == -1)
                _heap.insert( HeapEntryD (nId) );
              else // already exists 
                _heap.update( HeapEntryD( nId) );
            }
          }

        } // over neighs
      } // over grid Segments
      int numRed(0);
      while ( !_heap.empty() )
      {
        numRed++;

        HeapEntryD id = _heap.front();
        _heap.pop_front();
        int cSeg = id.pId();
        HeapStoreEntry& heapE = _heapEntryStore[ cSeg ];

        // fix segment , add neighbors if smaller
        int sSeg = heapE.startSeg();
        proposalMaps[startPos + cSeg ] = sSeg;

        // invalidate segments once and for all
        heapE.setStartSeg ( -2 );

        // loop neighs:
        Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, cSeg) );
        int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, cSeg) )/5;
        for (int j =0; j < nIds ;j++)
        {
          int nId = (int) (edge [5*j ])-1;

          HeapStoreEntry& hse = _heapEntryStore[nId];
          int bestSeg = hse.startSeg();

          // else already fixed
          if( bestSeg > -2) // replace if distance is smaller
          {
            Scalar cost = (centers[sSeg]-centers[nId]).sqrnorm();

            int heapid = hse.heapId();

            if (heapid<0 && bestSeg>0)
              printf("Weird\n");

            // if so enter or update heap:
            if (cost < hse.cost())
            {
              hse.setCost(cost);
              hse.setStartSeg(sSeg);
              if ( bestSeg == -1 )
                _heap.insert( HeapEntryD (nId) );
              else // already exists 
                _heap.update( HeapEntryD( nId ) );
            }
          }
          // doing some stupid stuff
          int abs = 0;

        } // over neighs
      } // while heap not empty
    } // all segments are in a proposal map

    generateLargeBBox();
  }


  int getNProposals() {return nSegments;};

  std::vector<P4i>& getLargeBBox() {  return largeBboxSeg; };

  std::vector< std::list< std::pair<int,int> > >& getPageSegPair() { return pageSegPair; };


private:


  /// defines the current situation from which the solution is expanded - do not confuse with setInitSegImg
  void setSegImg( int* _segImg )
  { 
    segImg    = _segImg;
  }
 
  /// compute the 2d centers of the patches
  void setCenters (Scalar* centers_) 
  { 
    centers.clear();
    centers.resize(nSegments);
    for (int i=0; i<nSegments ;i++)

    {
      P3 cc = Kl * P3(&centers_[3*i]); cc = cc / cc[2];
      centers[i] = P2( cc[0], cc[1] );
    }
  }

    /// select the pixel of the grid and the corresponding segments/proposals, setup the queue for expansion, returns a list of segments to start with
  void initGrid(int centralSegment, std::vector<int>& startSegs)
  {
    int startSeg = centralSegment;

    if (pickedSegments[startSeg] > 0)
      printf("ERROR LOCALPROPOSALHANDLER picked segment too often\n");

    int px = min(max(Scalar(0), floor(centers[centralSegment][0] - 0.5 )), Scalar(w));
    int py = min(max(Scalar(0), floor(centers[centralSegment][1] - 0.5 )), Scalar(h));

    P2i gridCenter( px, py );

    //startSegs = segImg[ px+py*w ];

    int sx = px-gridSize * ( px/gridSize );
    int sy = py-gridSize * ( py/gridSize );

    startSegs.clear();
    for (int ix = sx;ix < w;ix+=gridSize)
      for (int iy = sy;iy < h;iy+=gridSize)
      {
        int segId = segImg[ix*h + iy];
        startSegs.push_back(segId);
      }
  }

  void initGridJitter(int centralSegment, std::vector<int>& startSegs, const mxArray* edges_, const std::vector<int>& coveredSegs)
  {
    int startSeg = centralSegment;

    if (pickedSegments[startSeg] > 0)
      printf("ERROR LOCALPROPOSALHANDLER picked segment too often\n");

    int px = min(max(Scalar(0), floor(centers[centralSegment][0] - 0.5 )), Scalar(w));
    int py = min(max(Scalar(0), floor(centers[centralSegment][1] - 0.5 )), Scalar(h));

    P2i gridCenter( px, py );

    //startSegs = segImg[ px+py*w ];

    int sx = px-gridSize * ( px/gridSize );
    int sy = py-gridSize * ( py/gridSize );

    startSegs.clear();
    for (int ix = sx;ix < w;ix+=gridSize)
      for (int iy = sy;iy < h;iy+=gridSize)
      {
        int segId = segImg[ix*h + iy];
        
        // else look if neigh has better options
        if (coveredSegs[segId] > 0) 
        {
          Scalar bestDist(std::numeric_limits<Scalar>::max() );
          Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, segId) );
          int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, segId) )/5;
          for (int j =0; j < nIds ;j++)
          {
            int nId = (int) (edge [5*j ])-1;
            if (coveredSegs[nId] == 0) 
            {
              P2 gc (ix,iy);
              Scalar dist = (gc-centers[nId]).sqrnorm();
              if (bestDist > dist)
              {
                segId = nId;
                bestDist = dist;
              }
            }
          }
        }
        startSegs.push_back(segId);
      }
  }


    /// select the pixel of the grid and the corresponding segments/proposals, setup the queue for expansion, returns a list of segments to start with
  int countGrid(int centralSegment, const std::vector<int>& coveredSegs)
  {
    int freeSegs = 0;

    int startSeg = centralSegment;

    if (pickedSegments[startSeg] > 0)
      printf("ERROR LOCALPROPOSALHANDLER picked segment too often\n");

    int px = min(max(Scalar(0), floor(centers[centralSegment][0] - 0.5 )), Scalar(w));
    int py = min(max(Scalar(0), floor(centers[centralSegment][1] - 0.5 )), Scalar(h));

    P2i gridCenter( px, py );

    int sx = px-gridSize * ( px/gridSize );
    int sy = py-gridSize * ( py/gridSize );

    for (int ix = sx;ix < w;ix+=gridSize)
      for (int iy = sy;iy < h;iy+=gridSize)
      {
        int segId = segImg[ix*h + iy];
        if (coveredSegs[segId] == 0) 
          freeSegs++;
      }

    return freeSegs;
  }

  /// select the pixel of the grid and the corresponding segments/proposals, setup the queue for expansion, returns a list of segments to start with
  int countGridJitter(int centralSegment,  const mxArray* edges_, const std::vector<int>& coveredSegs)
  {
    int freeSegs = 0;

    int startSeg = centralSegment;

    if (pickedSegments[startSeg] > 0)
      printf("ERROR LOCALPROPOSALHANDLER picked segment too often\n");

    int px = min(max(Scalar(0), floor(centers[centralSegment][0] - 0.5 )), Scalar(w));
    int py = min(max(Scalar(0), floor(centers[centralSegment][1] - 0.5 )), Scalar(h));

    P2i gridCenter( px, py );

    int sx = px-gridSize * ( px/gridSize );
    int sy = py-gridSize * ( py/gridSize );

    for (int ix = sx;ix < w;ix+=gridSize)
      for (int iy = sy;iy < h;iy+=gridSize)
      {
        int segId = segImg[ix*h + iy];
        if (coveredSegs[segId] == 0) 
          freeSegs++;
        else
        {
          Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, segId) );
          int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, segId) )/5;
          for (int j =0; j < nIds ;j++)
          {
            int nId = (int) (edge [5*j ])-1;
            if (coveredSegs[nId] == 0) 
            {
              freeSegs++;
              break;
            }
          }
        }
      }

    return freeSegs;
  }

  void generateSmallBBox()
  {
    bboxSeg.clear();
    bboxSeg.resize( nSegments, P4i(w+1,h+1,-1,-1) );

    for (int i = 0;i<w*h; i++)
    {
      int px = i/h;
      int py = i%h;

      int segment = segImg[i];

      P4i& box = bboxSeg[segment];

      box[0] = min( box[0], px );
      box[1] = min( box[1], py );
      box[2] = max( box[2], px+1 );//+1 for a=0;a<end so set end to pix+1
      box[3] = max( box[3], py+1 );
    }
  }

  void generateLargeBBox()
  {

    pageSegPair.clear();
    pageSegPair.resize( nSegments );

    largeBboxSeg.clear();
    largeBboxSeg.resize( nSegments, P4i(w+1,h+1,-1,-1) );

    for (int i=0; i< proposalMaps.size(); i++)
    {
      int startSeg = proposalMaps[i];
      int segment = i % nSegments;

      pageSegPair[startSeg].push_back ( std::pair<int,int> (i/nSegments, segment) );

      P4i& box = largeBboxSeg[startSeg];

      box[0] = min( box[0], (bboxSeg[segment])[0] );
      box[1] = min( box[1], (bboxSeg[segment])[1] );
      box[2] = max( box[2], (bboxSeg[segment])[2] );
      box[3] = max( box[3], (bboxSeg[segment])[3] );
    }
  }

  ///////////////////////////////////////////////

  /// counts how manygrids were necessary: the fewer the better
  int nGrids;

  /// pointer to the 2d centers
  std::vector<P2> centers;

  M3 Kl;

  /// the size of the pathces, that amount of pixel can get the id of its assigned segment
  int gridSize;

  /// keeps track of the amount a segment occurs in the maps
  std::vector<int>  pickedSegments;

  /// all the proposal maps joint together
  std::vector<int>  proposalMaps;

  /// the image storing the segmentation
  int* segImg;

  /// width
  int w;
  /// height
  int h;

  /// the bounding box for a single segment
  std::vector<P4i> bboxSeg;

  /// the full bounding box for a single segment
  std::vector<P4i> largeBboxSeg;

  /// stores where the segment i got assigned to page n and segment m
  std::vector< std::list< std::pair<int,int> > > pageSegPair;

  /// the amount of segmetns in the images
  int nSegments;

    /// stores relevant information about HE collapse
  HeapEntryStore _heapEntryStore;

  /// the heap-Interface structure
  HeapInterfaceD < HeapEntryD, HeapEntryStore > _heapInterface;

  /// the heap for selecting the element with the least cost
  Utils::HeapT< HeapEntryD, HeapInterfaceD < HeapEntryD, HeapEntryStore > > _heap;
};

#endif // __PatchSegmentHandler__h
