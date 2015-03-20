////////////////////////////////////////////////////////////////
// handle occlusions for an image pair: z-Buffer with memory //
// generate lists of (potential) occluders                   //
////////////////////////////////////////////////////////////////

#ifndef __OcclusionMappingBuffer__cpp
#define __OcclusionMappingBuffer__cpp

#include "OcclusionMappingBuffer.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <set>
#include <limits>
#include <math.h>

#include "QPBO_Converter.h" // convert occusion configuration to binary terms

using namespace std;
using namespace Math;

/// an occlusion is only given if the depth deviation is more than 1 pct of the depth at the pixel
#define _OccPixPct_ 0.015

typedef std::pair <int,int>  PixSegPair;


/// the lists of the occPixelList's has to be sorted by the ownIds to be able to merge lists
template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
generateOcclusionLists( int segId, int aoix, int aoiX, int aoiy, int aoiY )
  {
    freeOcclusionsSeg0.clear();
    freeOcclusionsSeg1.clear();

    // extra part:
    outsideOccluded.clear();//list of indices of pixel outside occluded by free pixel 
    outerOcclusions.clear();

    for (int j=aoiy;j<aoiY;j++)
    {
      int pos = j*w;int pId = pos+startX;
      for (int ii=startX;ii<endX;ii++, pId++)
      {
        int locId = (*locIds)[(endX-startX)*(j-aoiy)+(ii-startX)]; 
        if (locId < 0) continue;
        occPixelList< Scalar > pixList0;
        occPixelList< Scalar > pixList1;
        pixList0.ownId = locId;
        pixList1.ownId = locId;
        // has to contain its own id for some reason
        pixList0.occluders.push_back(PixSegPair( locId, 0) );
        pixList1.occluders.push_back(PixSegPair( locId, 1) );
        int retVal = getOcclusionLists( pixList0, pixList1, pId, segId, aoix, aoiX, aoiy, aoiY, locId );
        if (retVal & 1)
          freeOcclusionsSeg0.push_back( pixList0 );
        if (retVal & 2)
          freeOcclusionsSeg1.push_back( pixList1 );
      }
    }

    std::set<int>::iterator s_it(outsideOccluded.begin()), s_end(outsideOccluded.end());
    for (; s_it != s_end;s_it++ )
    {
      int pId = *s_it;
      occPixelList< Scalar > outerPixelList;
      outerPixelList.ownId = pId;
      if( generateExtraOcclusions ( pId, outerPixelList, segId, aoix, aoiX, aoiy, aoiY ) )
        outerOcclusions.push_back( outerPixelList ); // copy
    }
  }


template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
init()
  {
    idX.resize(w*h);
    idY.resize(w*h);

    idXCurSeg.resize(w*h);
    idYCurSeg.resize(w*h);

    // the id map: idea: for 2x use in warping we need that for the original/reference image
    if (p2d == NULL || segImg==NULL)
    {
      for (int i=0;i<w*h;i++)
      {
        idX[i]       = 1+i%w;idY[i]       = 1+i/w;
        idXCurSeg[i] = 1+i%w;idYCurSeg[i] = 1+i/w;
      }
      return;
    }

    /// reserve some memory
    for (int i=0;i<w*h;i++)
    {
      store0[i].reserve(15);
      store1[i].reserve(15);
    }

    // precompute the current situation
    int lastSeg = -1;
    M3 H;P3 N;
    if (segImg != NULL && homs.size() == nSegments )
    {
      for (int i=0;i<w*h;i++)
      {
        int segI = segImg[i];
        if (lastSeg != segI)
        {
          lastSeg = segI;
          H = homs [segI];
          N = vNoms[segI];
        }
        /// get the new position
        P3 pixpos = H * P3(&p2d[3*((i%w)*h+int(i/w))]);pixpos /= pixpos[2];

        idX[i] = pixpos[0];
        idY[i] = pixpos[1];

        int px ( floor(pixpos[0]-0.5) );
        int py ( floor(pixpos[1]-0.5) );

        if (px>=0 && px < w && py>=0 && py < h)
        {
          int pixId = px+py*w;
          Scalar depth = (Scalar(-1.)/(N|P3(&p2d[3*((pixId%w)*h+int(pixId/w))])));

          store0[pixId].push_back( oMapBufferStorage<Scalar>( i, px,py, depth, segI ) );
          lookup0[i] = oMapBufferLookup<Scalar>( pixId, store0[pixId].size()-1 );
        }
        else
          lookup0[i] = oMapBufferLookup<Scalar>( -1, -1 );
      }
    }
    else
      printf("OcclusionMapBuffer not properly initialized\n");
  }

/// convert to directly usable stuff
template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
setHomvNom( Scalar* homs_, Scalar* vNoms_ )
  {
    assert(segImg != NULL);
    assert(nSegments > 0);

    homs.clear(); homs.resize(nSegments);
    vNoms.clear(); vNoms.resize(nSegments);

    for(int i =0; i< nSegments; i++)
    {
      homs [i] = M3(&homs_[i*9]);
      vNoms[i] = P3(&vNoms_[i*3]);
    }
  }

template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
setSegImg( int _w, int _h, int* _segImg, int _nSegments )
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


  ////////////////////////////////////////////////////////////////////////////
  // functions which are called several times
  ////////////////////////////////////////////////////////////////////////////

  /// clears the part out of the 
template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
initNewSeg( int segid1, int _startX, int _endX, int _startY, int _endY )
  {
    outsideOccluded.clear();

    startX = _startX;
    endX   = _endX;
    startY = _startY;
    endY   = _endY;

    if (p2d == NULL || segImg==NULL)
      return;

    newIteration++;
    M3 H = homs [segid1];
    P3 N = vNoms[segid1];

    if (genOccData)
    {
      // delete old entries first
      for (int j=startY;j<endY;j++)
      {
        int pos = j*w;int i = pos+startX;
        for (int ii=startX;ii<endX;ii++,i++)
        {
          if ( lookup1[i].pixPos >=0 )
            store1[lookup1[i].pixPos].clear();

          if ( lookup0[i].pixPos >=0 )
            store1[lookup0[i].pixPos].clear();

          lookup1[i] = oMapBufferLookup<Scalar>(-1,-1);
        }
      }
    }

    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++,i++)
      {
        int segI = segImg[i];
        /// get the new position
        P3 pixpos = H * P3( &p2d[3*((i%w)*h+int(i/w))] );pixpos /= pixpos[2]; 

        idXCurSeg[i] = pixpos[0];
        idYCurSeg[i] = pixpos[1];

        if (genOccData)
        {
          /// nothing new
          if (segid1 == segI) 
            continue;

          int px = floor(pixpos[0]-0.5);
          int py = floor(pixpos[1]-0.5);

          if (px>=0 && px < w && py>=0 && py < h)
          {
            int pixId = px+py*w;
            Scalar depth = (Scalar(-1.)/(N|P3(&p2d[3*((pixId%w)*h+int(pixId/w))])));

            // I only do assignments here, so there is NO WAY 
            // that there are pixel from a different segment in the map
            store1[pixId].push_back( oMapBufferStorage<Scalar>( i, px,py, depth, segid1 ) );
            lookup1[i] = oMapBufferLookup<Scalar>( pixId, store1[pixId].size()-1 );
          }
          //        else // oob
        }
      }
    }
  }

  // update : the other pixel have an id in the mrf stored here HERE updates must be 0: no update 1:changed
  /// clears some stuff and updates according to segment map, assume updates encodes a 0 for no change and 1 else
//  void finishNewSeg( int segid1, int startX, int endX, int startY, int endY, std::vector<int>& changes )
template<class Scalar>
void 
OcclusionMapBuffer<Scalar>::
finishNewSeg( int segid1, std::vector<int>& loc2glob, std::vector<int>& changes )
  {
    if (p2d == NULL || segImg==NULL)
      return;

    M3 H = homs [segid1];
    P3 N = vNoms[segid1];

    for (int ii=0;ii<loc2glob.size();ii++)
    {
      int i = loc2glob[ii];
      if (changes[ii] <=0 || lookup0[i].pixPos < 0) continue;

      if (genOccData)
      {
        std::vector< oMapBufferStorage<Scalar> >& tmp = store0[lookup0[i].pixPos];
        if (tmp.size() > 1)
        {
          int vp = lookup0[i].vecPos;
          tmp[vp] = tmp.back();
          lookup0[i].vecPos = -1;
          lookup0[i].pixPos = -1;
          lookup0[tmp[vp].pixPos].vecPos = vp;
          tmp.pop_back();
        }
        else
        {
          tmp.clear();
        }
      }
      /// get the new position
      P3 pixpos = H * P3( &p2d[3*((i%w)*h+int(i/w))] );pixpos /= pixpos[2];

      idX[i] = pixpos[0];
      idY[i] = pixpos[1];

      if (genOccData)
      {
        int px = floor(pixpos[0]-0.5);
        int py = floor(pixpos[1]-0.5);

        if (px>=0 && px < w && py>=0 && py < h)
        {
          int pixId = px+py*w;
          Scalar depth = (Scalar(-1.)/(N|P3(&p2d[3*((pixId%w)*h+int(pixId/w))])));

          store0[pixId].push_back( oMapBufferStorage<Scalar>( i, px,py, depth, segid1 ) );
          lookup0[i] = oMapBufferLookup<Scalar>( pixId, store0[pixId].size()-1 );
        }
        else
          lookup0[i] =oMapBufferLookup<Scalar>(); // redundant
      }
    }
  }

  // projPix: the pixel is projected onto that pixel, info the list of pixel assigned a 0, falling on the same pixel 
template<class Scalar>
bool
OcclusionMapBuffer<Scalar>::
constructList0( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY)//, int projPix)
  {
    // first the 0 part - original segments
    int pp = lookup0[pId].pixPos;
    int vp = lookup0[pId].vecPos;
    const std::vector <oMapBufferStorage<Scalar> >& info = store0[pp];

    ///--------------------------------------------------
    // to find out whether an outside pixel is occluded 
    // encodes whether there is a different occluder as a free variable, if so do not store the largest occluded outside pixel
    int maxFreeOccluder =  1; 
    int outerOccluded   = -1;
    Scalar outerOccludedDepth = Scalar( std::numeric_limits<Scalar>::max() );
    ///--------------------------------------------------

    const oMapBufferStorage<Scalar>& pix = info[vp];

    typename std::vector<oMapBufferStorage<Scalar> >::const_iterator it(info.begin()), it_end(info.end());

    for (; it != it_end; it++)
    {
      if ( pId == it->pixPos ) continue; // no self occlusions of segments

      int locPixId = isFree( *it, aoix, aoiX, aoiy, aoiY );
      if ( locPixId != -1)
      {
        if (it->segId == pix.segId) continue; // no self occlusions of segments

        bool occ0(false);
        // occlusion:
        if (pix.depth - _OccPixPct_ * pix.depth > it->depth && pix.segId != it->segId)
        {
          occ0 = true;
//          pixList.occluders.push_back( occPixelList<Scalar>::PixSegPair(it->pixPos ,0) );
          pixList.occluders.push_back( typename occPixelList<Scalar>::PixSegPair(locPixId,0) );
          maxFreeOccluder = 0;
        }

        // does the same pixel provoke also an occlusion when assigned to 1
        if (occ0 && lookup1[it->pixPos].pixPos == pp)
        {
          const oMapBufferStorage<Scalar> entry = (store1[ lookup1[it->pixPos].pixPos]) [lookup1[it->pixPos].vecPos];

          if ( entry.depth - _OccPixPct_ * pix.depth > entry.depth )
          {
            pixList.occluders.clear();
            pixList.fullOccluder = it->pixPos;
            return true;
          }
        }
      }
      else // potentially always occluded from outside, missing from inside
      {
        // always occluded from outside: but that is not a full occlusion - yes it is
        if (pix.depth - _OccPixPct_ * pix.depth > it->depth)// && (pix.segId != it->segId)) // done below already in the else case
        {
          pixList.occluders.clear();
          pixList.fullOccluder = it->pixPos;
          return true; // why not so far ? the pix is always occluded no matter what
          //maxFreeOccluder = 0; // not considered anyway, since we return
        }
        else  // always occludes a pixel outside, make sure that only the one with the largest depth is used here !
        {
          //
          // I NEED TO ADD A VARIABLE FOR THE LARGEST OUTSIDE PIXEL WHICH 
          // IS OCCLUDED AND ADD EDGES TO ENCODE A POTENTIAL (NON-)OCCLUSION
          // only store the id if our guy here is on top of all other free pixel in the list
          //
          // cases to check: another inside pixel occludes the same one
          // another outside pixel occludes the same one 
          // ensure that the largest one is considered here only (slanted surface: multiple largest possible)
          if ((pix.depth + _OccPixPct_ * pix.depth < it->depth) && (pix.segId != it->segId))
          {
            if ( it->depth < outerOccludedDepth )
            {
              outerOccluded      = it->pixPos;// this guy is potentially a pixel outside, which is occluded
              outerOccludedDepth = it->depth;
            }
          }
          else // same depth pixel exists outside -> never an occlusion (we can assume that fullOccluder == -1)
          {
            maxFreeOccluder    =  0;
            outerOccluded      = -1;
            outerOccludedDepth = it->depth;
          }
        }

        // Done:
        ////////////////// !!!!!!!!!!!! MISSIGN CASE TODO todo ///////////////////////
        // occlusion of a pixel outsied of the aoi, can be handled within the data term, (NO ACTUALLY NOT: add extra variables)
        // if a pixel gets assigned a new label an occlusion might be gone 
        // such that the data term at that position is back in place
        // encode by altering the data term of the occlusing pixel accordingly
        ////////////////// !!!!!!!!!!!! MISSIGN CASE TODO todo ///////////////////////
      }
    }

    // check out the assigned ones for occlusions - this store only contains the pixel assigned a 1 - new segment
    const std::vector <oMapBufferStorage<Scalar> >& info1 = store1[pp];
    it = info1.begin(); it_end = info1.end();

    for (; it != it_end; it++)
    {
      if (it->pixPos == pId) continue; // no pixel can occlude itself      
      if (it->segId == pix.segId) continue; // no pixel can occlude itself
      assert (it->segId != pix.segId); // no sense can only contain pixel from seg expSeg

      int locPixId = isFree( *it, aoix, aoiX, aoiy, aoiY );

      // occlusion:
      if (pix.depth - _OccPixPct_ * pix.depth > it->depth && (locPixId>=0))
        pixList.occluders.push_back(typename occPixelList<Scalar>::PixSegPair(locPixId ,1) );
    }

    // outside pixel which is occluded by free variables
    if ( maxFreeOccluder && outerOccluded >= 0 )
      outsideOccluded.insert( outerOccluded );

    if (outerOccluded >= 0)
      pixList.occludes = outerOccluded;

    return (pixList.occluders.size() > 1);
  }

  // projPix: the pixel is projected onto that pixel, info the list of pixel assigned a 0, falling on the same pixel 
template<class Scalar>
bool
OcclusionMapBuffer<Scalar>::
constructList1( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY)//, int projPix)
  {
    // pixel occludes some other pixel which would else be visible:
    //int potOcclusion = -1;
    Scalar potDepth  = 100000000;

    // ony the 0 part - original segments
    int pp = lookup1[pId].pixPos;
    int vp = lookup1[pId].vecPos;

    M3 H = homs [expSeg];
    P3 pixpos = H * P3( &p2d[3*((pId%w)*h+int(pId/w))] );pixpos /= pixpos[2];
    int px = floor(pixpos[0]-0.5);
    int py = floor(pixpos[1]-0.5);

    if (pp<0) return false;

    const std::vector <oMapBufferStorage<Scalar> >& info = store0[pp];

    const oMapBufferStorage<Scalar>& pix = (store1[pp])[vp];
    typename std::vector<oMapBufferStorage<Scalar> >::const_iterator it(info.begin()), it_end(info.end());

    ///--------------------------------------------------
    // to find out whether an outside pixel is occluded 
    // encodes whether there is a different occluder as a free variable, if so do not store the closest occluded outside pixel
    // occluded by the closest free pixel
    int maxFreeOccluder =  1; // is there an inner/free pixel which is closer to the camera ?
    int outerOccluded   = -1;
    Scalar outerOccludedDepth = Scalar( std::numeric_limits<Scalar>::max() ); // the depth of the occluded pixel, only consider the closest one
    ///--------------------------------------------------

    assert(expSeg == pix.segId);

    for (; it != it_end; it++)
    {
      if (it->pixPos == pId) continue; // no pixel can occlude itself, no self occlusions of segments
       int locPixId = isFree( *it, aoix, aoiX, aoiy, aoiY );
      if (locPixId != -1)
      {
        if (it->segId == pix.segId) continue; // no self occlusions of segments

        assert ((expSeg != it->segId)); // since any free pixel can not be from expSeg

        // occlusion:
        if ((pix.depth - _OccPixPct_ * pix.depth > it->depth) && (expSeg != it->segId))
        {
          pixList.occluders.push_back(typename occPixelList<Scalar>::PixSegPair(locPixId,0) );
          maxFreeOccluder = 0;
        }
      }
      else // potentially always occluded from outside, missing from inside
      {
        // always occluded: there is an occlusion from outside
        if ( pix.depth - _OccPixPct_ * pix.depth > it->depth && (expSeg != it->segId) )
        {
          pixList.occluders.clear();
          pixList.fullOccluder = it->pixPos;
          return true;
        }
        else // occludes a pixel outside
        {
          if ((pix.depth + _OccPixPct_ * pix.depth < it->depth) && (pix.segId != it->segId))
          {
            if ( it->depth < outerOccludedDepth )
            {
              outerOccluded      = it->pixPos;// this guy is potentially a pixel outside, which is occluded
              outerOccludedDepth = it->depth;
            }
          }
          else // same depth pixel exists -> never an occlusion
          {
            maxFreeOccluder    =  0;
            outerOccluded      = -1;
            outerOccludedDepth = it->depth;
          }
        }
      }
    }

    int testBreak = 0;
    if (outerOccluded == 148111)
      testBreak = 999;
    
    if ( maxFreeOccluder && outerOccluded>=0 )
      outsideOccluded.insert( outerOccluded );

    if (outerOccluded >= 0)
      pixList.occludes = outerOccluded;

    return (pixList.occluders.size() > 1); // 1 since we add the 'identity' at the beginning 
  }

/// assumes that pId is already the id of a pixel closest to the camera of all fixed pixel and occluded by one (or more) free pixel ! 
template<class Scalar>
bool
OcclusionMapBuffer<Scalar>::
generateExtraOcclusions ( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY )
  {
    // first the 0 part - original segments
    const int pp = lookup0[pId].pixPos;
    const int vp = lookup0[pId].vecPos;
    const std::vector <oMapBufferStorage<Scalar> >& info = store0[pp];

    pixList.ownId = pId;

    const oMapBufferStorage<Scalar>& pix = info[vp];
    typename std::vector<oMapBufferStorage<Scalar> >::const_iterator it(info.begin()), it_end(info.end());

    // only consider free pixel which potentially occluded this pixel: case 
    for (; it != it_end; it++)
    {
      if ( it->segId == pix.segId || it->pixPos == pId ) continue; // no self occlusions of segments
      int locPixId = isFree( *it, aoix, aoiX, aoiy, aoiY );

      if ( locPixId == -1 ) 
      {
        assert( pix.depth <= it->depth );
        continue;
      }
      // occlusion by a free pixel:
      if (pix.depth - _OccPixPct_ * pix.depth > it->depth )
      {
        pixList.occluders.push_back( typename occPixelList<Scalar>::PixSegPair(locPixId,0) );

        // does the same pixel provoke also an occlusion when assigned to 1
        if (lookup1[it->pixPos].pixPos == pp)
        {
          const oMapBufferStorage<Scalar> entry = (store1[ lookup1[it->pixPos].pixPos]) [lookup1[it->pixPos].vecPos];

          // this is a indifferent case since the outside pixel is occluded no matter what th eassignments are 
          // therefore we can omit any action here
          if ( pix.depth - _OccPixPct_ * pix.depth > entry.depth && entry.segId != pix.segId )
          {
            pixList.occluders.clear();
            pixList.fullOccluder = it->pixPos;
            return false;
          }
        }
      }
    }

    // no selfocclusion
    if (expSeg == pix.segId)
      return (pixList.occluders.size() >= 1);

    // segments which are assigned a one:
    // check out the assigned ones for occlusions - this store only contains the pixel assigned a 1 - new segment
    const std::vector <oMapBufferStorage<Scalar> >& info1 = store1[pp];
    it = info1.begin(); it_end = info1.end();

    for (; it != it_end; it++)
    {
      int locPixId = isFree( *it, aoix, aoiX, aoiy, aoiY );
      // occlusion: by free pixel
      if ((pix.depth - _OccPixPct_ * pix.depth > it->depth) && (locPixId>=0))
        pixList.occluders.push_back( typename occPixelList<Scalar>::PixSegPair(locPixId ,1) );
    }
    return (pixList.occluders.size() >= 1);
  }

  /// returns a list for each pixel asked, expSeg is to check if a pixel within the update area aoi is fixed
template<class Scalar>
int
OcclusionMapBuffer<Scalar>::
getOcclusionLists( occPixelList< Scalar >& pixList0, occPixelList< Scalar >& pixList1, int pId, int expSeg, int aoix, int aoiX, int aoiy, int aoiY, int locId )
  {
    if (p2d == NULL || segImg==NULL)
      return 0;

    int retVal(0);

    ///////////// old segment assignement:
    int pp = lookup0[pId].pixPos;
    int vp = lookup0[pId].vecPos;

    /// oob projection
    if (pp < 0 || vp < 0) return 0;

    const std::vector <oMapBufferStorage<Scalar> >& info = store0[pp];
    assert( info[vp].pixPos == pId );

    int segId = info[vp].segId;

    // now run along info vector and accumulate occlusions - or check for 
    if (constructList0( pId, pixList0, expSeg, aoix, aoiX, aoiy, aoiY))
      retVal +=1;
    if ( (segId != expSeg) && constructList1( pId, pixList1, expSeg, aoix, aoiX, aoiy, aoiY))
      retVal +=2;

    return retVal;
    /////////////////////////////////////////////////
  }

template<class Scalar>
void
OcclusionMapBuffer<Scalar>::
getCurrentOcclusions( Scalar* oMap )
  {
    memset( oMap, 0, sizeof(Scalar) * w * h );

    for (int j=0;j<h;j++)
    {
      int pos = j*w;int i = pos;
      for (int ii=0;ii<w;ii++,i++)
      {
        if ( lookup0[i].pixPos >=0 )
        {
          std::vector< oMapBufferStorage<Scalar> > info = store0[lookup0[i].pixPos];
          oMapBufferStorage<Scalar> pix = store0[lookup0[i].pixPos][lookup0[i].vecPos];

          typename std::vector<oMapBufferStorage<Scalar> >::const_iterator it(info.begin()), it_end(info.end());
          for (; it != it_end; it++)
            if (pix.depth - _OccPixPct_ * pix.depth > it->depth && pix.segId != it->segId)
            {
              oMap[i] = 1;//occluded
              break;
            }

        }
        else
          oMap[i] = Scalar(-1); // oob
      }
    }
  }

/*
// debugging:
template<class Scalar>
void
OcclusionMapBuffer<Scalar>::
getDepthMap( Scalar* outDepth )
{
    for (int i=0 ; i<w*h; i++)
      if ( lookup0[i].pixPos >= 0 )
        outDepth[i] = 100000;
    
    for (int i=0 ;i<w*h; i++)
      if ( lookup0[i].pixPos >= 0 )
      {
        outDepth[i] = min( outDepth[i], (store0[lookup0[i].pixPos])[lookup0[i].vecPos].depth );
        if ( lookup0[i].vecPos >= store0[lookup0[i].pixPos].size() )
          printf("Impossible \n");
      }
  }
*/

#undef _OccPixPct_
#endif // __OcclusionMapping__h
