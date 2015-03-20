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
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#ifndef __OcclusionMappingBufferSeg__h
#define __OcclusionMappingBufferSeg__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <set>
#include <limits>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Math;

/// an occlusion is only given if the depth deviation is more or less than 1 pct of the depth at the pixel
#define _OccSegPct_ 0.05

typedef std::pair <int,int>  PixSegPair;

/// the list of the occPixelList's has to be sorted by the ownIds to be able to merge lists
template<typename Scalar>
struct occSegList
{
  typedef std::pair <int,int>  SegLabelPair;
  occSegList( ) : mpId(-1), segId(-1), label(-1), fullOcclusion(0), score(0)
  {};

  occSegList( int _mpid, int _segid, int _label, std::vector < SegLabelPair > _occluders, int _fullOcclusion )
    : mpId(_mpid), segId(_segid), label(_label), occluders(_occluders), fullOcclusion(_fullOcclusion), score(0)
  {};

  ~occSegList() {};
  ////////////////////////////////////////////

  bool operator<(const occSegList<Scalar>& _rhs) const 
  {
    return (segId < _rhs.segId) || ( (segId == _rhs.segId) && (label < _rhs.label));
  }

  bool operator>(const occSegList<Scalar>& _rhs) const 
  {
    return (segId > _rhs.segId) || ( (segId == _rhs.segId) && (label > _rhs.label));
  }

  bool operator==(const occSegList<Scalar>& _rhs) const 
  {
    return ( (segId == _rhs.segId) && (label == _rhs.label));
  }

  bool isLabel0OccludedbyLabel0()
  {
    if (label) return false;
    if (fullOcclusion) return true;

    for(int i=0;i<occluders.size();i++)
      if (occluders[i].second == 0)
        return true;

    return false;
  }

  /// stored are pairs: pixel id, assign a 0 or 1, interpretation 0: occludes if id=0, encode x_i in poly: if x_i=1 not occluded
  std::vector < SegLabelPair > occluders;

  int mpId;        // id of the moving plane
  int segId;       // segment id
  int label;       // label
  int fullOcclusion;

  Scalar score;
};

template<typename Scalar>
struct patchEntry
{
  public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;

  patchEntry() : label(-1), pPix(0., 0., 1.), depth(0), segId(-1), hom() , iHom(), 
    iHom_otherLabel(), mp_id(-1), mpi_id(-1), iDepth(0), 
    normal(0,0,1), normalI(0,0,1) {};

  // with segid interesting for the lookup of new id assignments (a 1 to a pixel)
  patchEntry( int _label, int _segId, M3 _Hom , M3 _iHom, Scalar _px, Scalar _py, Scalar _depth, 
              M3 _iHom_otherLabel, int _pid, int _ipid, Scalar _iDepth, P3 _normal, P3 _normalI)
    : label(_label), pPix(_px, _py, 1.), depth(_depth), segId(_segId), hom(_Hom) , iHom(_iHom ), 
    iHom_otherLabel(_iHom_otherLabel), mp_id(_pid), mpi_id(_ipid), iDepth(_iDepth), 
    normal(_normal), normalI(_normalI)
  {};

  ~patchEntry() {};
  ///////////////////////////////////////
  P3 pPix;
  int label;
  Scalar depth;
  int segId;
  M3 hom;
  M3 iHom;
  M3 iHom_otherLabel; // idea: fast check if center is occluded by segment with label 0 and 1 -> always occluded
  int mp_id;//id of the moving plane
  int mpi_id;//id of the moving plane when assigned other laberl
  Scalar iDepth; // depth of the other label assigne,tn
  /// normal of the trafo patch
  P3 normal;
  /// normal of the other label patch
  P3 normalI;
};
//////////////////////////////////////////////////////////////////////////

template<typename Scalar>
bool sortpatchEntry (patchEntry<Scalar> i, patchEntry<Scalar> j) { return (i.pPix[0] < j.pPix[0]); };

template<typename Scalar>
bool sortOccluders ( std::pair <int,int> i, std::pair <int,int> j) { return (i.first < j.first) || ((i.first == j.first) && (i.second < j.second)); };

// idea:
// forward map all! segment centers, sort these by x value
// check for occlusions by testing all projected cneters in a 50 pixel neighborhood
// by using the backward homography = invert it
//
// if the 3d center of plane A (segment sA) projects back to 
// segment sB using plane B there can be an occlusion (if depth pB < depth plane A)
//
// the list discriminates submodular edges -> at once conversion
// and other edges: accumulated and converted by new means
///////////////////////////////////////
template<typename Scalar> class OcclusionMapBufferSeg
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;

  typedef std::pair <int,int>       Seg_ij;
  typedef std::pair <int,int>       SegLbl;

  OcclusionMapBufferSeg( int w_, int h_, genHomography<Scalar>* _gHom, int* _segImg, Scalar* _centers, Scalar* _Kl)
    : w(w_), h(h_), occPenalty(0), p2d(NULL), segImg(_segImg), centers(_centers), gHom (_gHom), Kl(_Kl)
  {
    iKl  = Kl;iKl.invert();
    iKlt = iKl.transpose();
  }

  ~OcclusionMapBufferSeg(){};

  std::vector<SegLbl>& getFullOcc() {return fullOcclusions;}
  /// list of all occlusions observed for the view
  std::vector< occSegList<Scalar> >& getOccList() 
  {
    return  occlusionList;
  };

  /// first project all centers and sort them
  void buildList ( std::vector<int>& currentSolution, std::vector<int>& trialVec )
  {
    projPatchCenters.clear();
    projPatchCenters.resize( 2*currentSolution.size() );

    // first project current solutions:
    for (int i=0;i<currentSolution.size();i++)
    {
      int pid   = currentSolution[i];
      int pid_t = trialVec[i]; 
      // from pixel to pixel:
      M3 Hom_0  = gHom->getHom_n( pid, i );
      M3 Hom_1  = gHom->getHom_n( pid_t, i );
      M3 iHom_0 = Hom_0;iHom_0.invert();
      M3 iHom_1 = Hom_1;iHom_1.invert();
      P3 p0     = Hom_0 * Kl * P3( &centers[3*i] );// center is in camera coords
      P3 p1     = Hom_1 * Kl * P3( &centers[3*i] );

      // should work on pixel: 1. image to camera coords 2. scalar product 3.invert
      P3 vn_0_  = gHom->getViewNormal( pid, i );
      P3 vn_0  = iKlt * gHom->getViewNormal( pid, i ); // this is vn_o' * p == normal' * (IKlt*p), since IKlt is transposed
      P3 vn_1  = iKlt * gHom->getViewNormal( pid_t, i );      

      p0 /= p0[2];
      p1 /= p1[2];

      Scalar depth  = 1/(vn_0|p0);
      Scalar idepth = 1/(vn_1|p1);

      projPatchCenters[2*i]   = patchEntry<Scalar>(0, i, Hom_0, iHom_0, p0[0], p0[1], depth, iHom_1, pid, pid_t, idepth, vn_0, vn_1);
      projPatchCenters[2*i+1] = patchEntry<Scalar>(1, i, Hom_1, iHom_1, p1[0], p1[1], idepth, iHom_0, pid_t, pid, depth, vn_1, vn_0);
    }

    std::sort( projPatchCenters.begin(), projPatchCenters.end(), sortpatchEntry<Scalar> );
  }

  void generateOcclusionLists()
  {
    const Scalar maxDiff   (50.*2);
    const Scalar maxDiff_2 (maxDiff*maxDiff/2);

    fullOcclusions.clear();
    occlusionList.clear();

    for (int i=0;i<projPatchCenters.size();i++)
    {
      patchEntry<Scalar>& pe = projPatchCenters[i];
      P3 &start = pe.pPix;
      int j=i-1;int k=i+1;

      // build the oclcusionList inserting stuff:
      occSegList< Scalar > pixList;
      pixList.mpId  = pe.mp_id;
      pixList.segId = pe.segId;
      pixList.label = pe.label;
      // has to contain its own id for some reason
      // pixList.occluders.push_back(PixSegPair( pe.segId, pe.label ) );

      bool fullOcc = false;
      // backwards:
      while( (j>=0) && (fabs(start[0]- projPatchCenters[j].pPix[0]) < maxDiff) && !fullOcc )
      {
        patchEntry<Scalar> & po = projPatchCenters[j];
        P3 &other = po.pPix;

        if (pe.segId == po.segId) // not any sense in this, one label per segi
        {
          j--;
          continue;
        }
        // normal oclcusion
        int occlusionCheck = occlusionCheckHom(start, other, pe, po, maxDiff, maxDiff_2 );
        if (occlusionCheck == 1) // normal oclcusion case : insert into list
           pixList.occluders.push_back(PixSegPair( po.segId, po.label ) );

        // full occlusion: segment always occluded
        if (occlusionCheck == 2) 
        {
          pixList.fullOcclusion = 1;
//          pixList.occluders.clear(); // temporarily
          fullOcclusions.push_back( SegLbl( pe.segId, pe.label ) );fullOcc =true;
        };

        j--;
      }

      // backwards:
      while( (k<projPatchCenters.size()) && (fabs(start[0]- projPatchCenters[k].pPix[0]) < maxDiff) && !fullOcc )
      {
        patchEntry<Scalar> & po = projPatchCenters[k];
        P3 &other = po.pPix;

        if (pe.segId == po.segId) // not any sense in this only one active 
        {
          k++;
          continue;
        }
          
        // normal oclcusion
        int occlusionCheck = occlusionCheckHom(start, other, pe, po, maxDiff, maxDiff_2 );
        if (occlusionCheck == 1) // normal oclcusion case : insert into list
           pixList.occluders.push_back(PixSegPair( po.segId, po.label ) );

        // full occlusion: segment always occluded
        if (occlusionCheck == 2)
        {
          pixList.fullOcclusion = 1;
          fullOcclusions.push_back( SegLbl( pe.segId, pe.label ) );fullOcc =true;
        };

        k++;
      }
      ///

      if ( (!fullOcc && pixList.occluders.size() > 0 ) || pixList.fullOcclusion )
      {
        std::sort( pixList.occluders.begin(), pixList.occluders.end(), sortOccluders<Scalar> );
        occlusionList.push_back(pixList);
      }
    }

    //  must be sorted otherwise merging process fails outside (one could sort there but then the list is sorted twice)
    std::sort( occlusionList.begin(), occlusionList.end() );
  }


  /// actually the depth check must be done differnetly - storing the trafo normal and scalar product
  int occlusionCheckHom (P3& start, P3& other, patchEntry<Scalar>& pe, patchEntry<Scalar>& po, Scalar maxDiff, Scalar maxDiff_2 )
  {
    if ( (start[0]- other[0])*(start[0]- other[0]) + (start[1]- other[1])*(start[1]- other[1]) < maxDiff_2)
    {
      // check back projection: in segment ?
      P3 backProj = po.iHom * start;backProj[0] /= backProj[2];backProj[1] /= backProj[2];

      // round:
      int ppx = floor( backProj[0]-0.5 );
      int ppy = floor( backProj[1]-0.5 );

      Scalar otherDepth = Scalar(1.)/(po.normal|start);

      // TODO
      // projects on original segment: segImg[ ppy*w + ppx ] == po.segId 
      // from a differnet moving plane: pe.mp_id != po.mp_id -> problems if restartet !!, covered by _OccSegPct_
      // depth occludes one another and distance larger than threshold

      if ( ppx>=0 && ppx<=w-1 && ppy>=0 && ppy<=h-1 && segImg[ ppy*w + ppx ] == po.segId &&
           //pe.mp_id != po.mp_id && 
           (pe.depth - _OccSegPct_ * pe.depth > Scalar(1.)/(po.normal|start)) ) // flipped normals in smoothness eval
      {
        P3 backProj_L = po.iHom_otherLabel * start;backProj_L[0] /= backProj_L[2];backProj_L[1] /= backProj_L[2];
        int ppx_L = floor( backProj_L[0]-0.5 );
        int ppy_L = floor( backProj_L[1]-0.5 );
        if ( ppx_L>=0 && ppx_L<=w-1 && ppy_L>=0 && ppy_L<=h-1 && segImg[ ppy_L*w + ppx_L ] == po.segId && 
          //pe.mp_id != po.mpi_id && 
          (pe.depth - _OccSegPct_ * pe.depth > Scalar(1.)/(po.normalI|start))) // flipped normals
        {
          return 2;// fullOcclusion: segment occluded anyway
        }
        return 1; // simple occlusion
      }
    }
    return 0; // no occlusion
  }

  void SetOccPenalty (Scalar _occPenalty) {occPenalty = _occPenalty;};

  /// first project all centers and sort them
  void buildSingleList ( std::vector<int>& currentSolution )
  {
    projPatchCenters.clear();
    projPatchCenters.resize( currentSolution.size() );

    // first project current solutions:
    for (int i=0;i<currentSolution.size();i++)
    {
      int pid   = currentSolution[i];
      // from pixel to pixel:
      M3 Hom_0  = gHom->getHom_n( pid, i );
      M3 iHom_0 = Hom_0;iHom_0.invert();
      P3 p0     = Hom_0 * Kl * P3( &centers[3*i] );// center is in camera coords

      // should work on pixel: 1. image to camra coords 2. scalar product 3.invert
      P3 vn_0_  = gHom->getViewNormal( pid, i );
      P3 vn_0  = iKlt * gHom->getViewNormal( pid, i );

      p0 /= p0[2];

      Scalar depth  = 1/(vn_0|p0);

      P3 p0_cam = iKl * p0;
      Scalar  wft_depth = 1/(vn_0_|p0_cam);

      projPatchCenters[i]   = patchEntry<Scalar>(0, i, Hom_0, iHom_0, p0[0], p0[1], depth, iHom_0, pid, pid, depth, vn_0, vn_0);
    }

    std::sort( projPatchCenters.begin(), projPatchCenters.end(), sortpatchEntry<Scalar> );
  }



private:

  ///////////////////////////////////////////////

  /// temporal list for fast z buffering
  std::vector<patchEntry<Scalar> > projPatchCenters;

  genHomography<Scalar>* gHom;

  /// width
  int w;
  /// height
  int h;

  M3 Kl;
  M3 iKl;
  M3 iKlt;

  Scalar *p2d;

  /// the segment ids for later test w.r.t. occlusions
  int* segImg;

  // life would be simple - here recomputed inside
//  std::vector<M3> homs; // homographies
//  std::vector<P3> vNoms;// view normals, to compute the depth at that pixel quickly

  Scalar occPenalty;

  /// a segment is occluded if the center is:
  Scalar* centers;

  /// vector of all segments which are anyway occluded
  std::vector<SegLbl> fullOcclusions;
  /// list of all occlusions observed for the view
  std::vector< occSegList<Scalar> > occlusionList;

  /// combined list (occlusionList and fullOcclusions) for conveniance
  std::vector< occSegList<Scalar> > fullOcclusionList;
};

#undef _OccSegPct_
#endif // __OcclusionMapping__h
