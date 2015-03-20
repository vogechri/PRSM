/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

#ifndef __OMapBufVC__cpp
#define __OMapBufVC__cpp

#include "OcclusionMappingBufferVC.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include "SegmentGrowing_HeapHelper.h"
#include "Templates/HeapT.h"

#include "genHom.h"
#include "mex.h" // include after VectorT and Mat3x3T

using namespace std;
using namespace Math;
///////////////////////////////////////

// rotation centered at segment center not in origin
#define _centeredViews_
#define  __impPenalty__ 200.

/// does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled: p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
template<class Scalar>
void OcclusionMapBufferVC<Scalar>::
extendSegmentation( genHomoG<Scalar>* gHom, const std::vector<Math::VectorT<Scalar, 3> >* normals, const mxArray* edges2_, Scalar *centers_, int nSegments, Scalar *centers2_, const std::vector<int>& solution, int* segImg, int* segImg2, std::vector<int>& solution2, int nSegments2, Scalar* Kl_, Scalar* Kl2_  )
  {
    solution2.resize(nSegments2, -1);
    std::vector<Scalar> depth ( nSegments2, 1000000000 );
    std::vector<Scalar> scores( nSegments2, 1000000000 );
    M3 Kl(Kl_);
    M3 Kl2(Kl2_);
    M3 iKl  = Kl;iKl.invert();

    M3 iKlt = iKl;iKlt=iKlt.transpose();

    // if second view solution check is on:
    // hom maps from to pixel but nom acts on image world coords
    std::vector<P3> normals2nd;
    if (gHom2nd != NULL )
    {
    normals2nd.resize(normals->size());
    for (int i=0;i< normals->size();i++)
    {
#ifdef _centeredViews_
      P3 vn_0_  = gHom2nd->getViewNormalC0( i );
//      P3 vn_0  = iKlt * gHom2nd->getViewNormalC0( i );
#else
      P3 vn_0_  = gHom2nd->getViewNormal( i, i );
//      P3 vn_0  = iKlt * gHom2nd->getViewNormal( i, i );
#endif
      normals2nd[i] = vn_0_;
    }
    }
    std::vector<P3> normals3rd;
    if (gHom3rd != NULL )
    {
    normals3rd.resize( normals->size());
    for (int i=0;i< normals3rd.size();i++)
    {
#ifdef _centeredViews_
      P3 vn_0_  = gHom3rd->getViewNormalC0( i );
//      P3 vn_0  = iKlt * gHom2nd->getViewNormalC0( i );
#else
      P3 vn_0_  = gHom2nd->getViewNormal( i, i );
//      P3 vn_0  = iKlt * gHom2nd->getViewNormal( i, i );
#endif
      normals3rd[i] = vn_0_;
    }
    }

    ////// step1: project forward. and keep the one with lowest depth !? later push relentlessly
    for ( int i = 0; i < nSegments ;i++ )
    {
      int proposal = solution[i];
      P3 centerA = P3( &centers_[ 3*i ] );
#ifdef _centeredViews_
      M3 Hom_0  = gHom->getHomC0( proposal );// getHom is also ok if normals are not negative
#else
      M3 Hom_0  = gHom->getHom( proposal, i );// getHom is also ok if normals are not negative
#endif
      P3 pp = Kl * centerA;pp /= pp[2];

      M3 iHom_0 = Hom_0;iHom_0.invert();
      P3 p0     = Hom_0 * Kl * P3( &centers_[3*i] );// center is in camera coords

      // should work on pixel: 1. image to camra coords 2. scalar product 3.invert
#ifdef _centeredViews_
      P3 vn_0_  = gHom->getViewNormalC0( proposal );
      P3 vn_0  = iKlt * gHom->getViewNormalC0( proposal );
#else
      P3 vn_0_  = gHom->getViewNormal( proposal, i );
      P3 vn_0  = iKlt * gHom->getViewNormal( proposal, i );
#endif

      p0 /= p0[2];

      Scalar depth_p   = 1/(vn_0|p0);
#ifdef _DEBUG
      Scalar depth_p_  = 1/(vn_0_|p0);

      P3 p0_cam = iKl * p0;
      Scalar  wft_depth = 1./(vn_0_|p0_cam);

      // or control: also check: 
      Scalar depth_p_orig   = 1/( (*normals)[proposal] | (P3( &centers_[3*i] )) );
      int wrong = 0 ;
      if ( fabs(depth_p_orig - depth_p) > 0.1)
        wrong =1;
#endif

//      int pixId = int( p0[0]-0.5 ) +int ( p0[1]-0.5 )*w;
      int pixId = int( p0[0]-0.5 )*h +int ( p0[1]-0.5 );

      int inBounds = (int( p0[0]-0.5 ) >= 0) && (int( p0[0]-0.5 ) < w) && (int( p0[1]-0.5 ) >= 0) && (int( p0[1]-0.5 ) < h);
      if (!inBounds) continue; 
      int seg2 = segImg2 [pixId];
      // lookup depth:
      if (depth_p < depth[seg2]) // new foreground segment
      {
        depth[seg2]     = depth_p;
        solution2[seg2] = proposal;

        Scalar penalty = 0;

        if (gHom2nd != NULL)
        {
           int dummy=0;
           M3 Hom_0_2nd  = gHom2nd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 2nd
           // i need from solution filled view to 2nd:
           // to first then to 2nd
           Hom_0_2nd = Hom_0_2nd * iHom_0;
           penalty = getPenalty ( seg2, proposal, Hom_0_2nd, Kl2 * P3( &centers2_[ 3*seg2 ] ), solution2nd, &normals2nd, centers_, dummy);
        };
        if (gHom3rd != NULL)
        {
          int dummy=0;
           M3 Hom_0_3rd  = gHom3rd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 3rd
           // i need from solution filled view to 3rd:
           // to first then to 3rd
           Hom_0_3rd = Hom_0_3rd * iHom_0;
           penalty += getPenalty ( seg2, proposal, Hom_0_3rd, Kl2 * P3( &centers2_[ 3*seg2 ] ), solution3rd, &normals3rd, centers_, dummy);
        };

        scores[seg2]    = penalty;
      }
      // push proposals anyway into Q - due to occlusion? neigh could be of different nature?
      // but what with oob == occ - > need oracle, like 
    }
    //////////////////////
//    return;
    typedef HeapStoreEntryTg<Scalar>       HeapStoreEntry;
    typedef std::vector< HeapStoreEntry >  HeapEntryStore;
    typedef HeapInterfaceDg < HeapEntryG,  HeapEntryStore > HeapInterface;

    /// stores relevant information about HE collapse
    HeapEntryStore _heapEntryStore;

    /// the heap-Interface structure
    HeapInterface _heapInterface(&_heapEntryStore);

    /// the heap for selecting the element with the least cost
    Utils::HeapT< HeapEntryG, HeapInterface > _heap(_heapInterface);

    _heapEntryStore.reserve( 10*nSegments2);//, HeapStoreEntry(0, -1, -1) );

    // push into Q those which are not marked
    for ( int i = 0; i < nSegments ;i++ )
    {
      int proposal = solution[i];

      P3 centerA = P3( &centers_[ 3*i ] );
#ifdef _centeredViews_
      M3 Hom_0  = gHom->getHomC0( proposal );// getHom is also ok if normals are not negative
#else
      M3 Hom_0  = gHom->getHom( proposal, i );// getHom is also ok if normals are not negative
#endif
      M3 iHom_0 = Hom_0;iHom_0.invert();
      P3 p0     = Hom_0 * Kl * P3( &centers_[3*i] );// center is in camera coords

      // should work on pixel: 1. image to camra coords 2. scalar product 3.invert
#ifdef _centeredViews_
      P3 vn_0_  = gHom->getViewNormalC0( proposal );
      P3 vn_0  = iKlt * gHom->getViewNormalC0( proposal );
#else
      P3 vn_0_  = gHom->getViewNormal( proposal, i );
      P3 vn_0  = iKlt * gHom->getViewNormal( proposal, i );
#endif

      p0 /= p0[2];

//      int pixId2 = int( p0[0]-0.5 ) +int ( p0[1]-0.5 )*w;
      int pixId2 = int( p0[0]-0.5 )*h +int ( p0[1]-0.5 );

      // ignore oob problems
      int inBounds = (int( p0[0]-0.5 ) >= 0) && (int( p0[0]-0.5 ) < w) && (int( p0[1]-0.5 ) >= 0) && (int( p0[1]-0.5 ) < h);
      if (! inBounds ) continue;

      int seg2 = segImg2 [pixId2];// segment in the second image, which is already set!

      // later:
      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges2_, seg2) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges2_, seg2) );

      // neighs in 2nd view:
      for (int j=0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        // project back
        P3 centerB = P3( &centers2_[ 3*id ] );

        // now get a score when projecting back from 2nd view:
        int dummy=1;
        Scalar penalty = getPenalty ( id, proposal, iHom_0, Kl2 * centerB, solution, normals, centers_, dummy);

        if (gHom2nd != NULL)
        {
           M3 Hom_0_2nd  = gHom2nd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 2nd
           // i need from solution filled view to 2nd:
           // to first then to 2nd
           Hom_0_2nd = Hom_0_2nd * iHom_0;
           penalty += getPenalty ( id, proposal, Hom_0_2nd, Kl2 * centerB, solution2nd, &normals2nd, centers_, dummy);
        };
        if (gHom3rd != NULL)
        {
           M3 Hom_0_3rd  = gHom3rd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 3rd
           // i need from solution filled view to 3rd:
           // to first then to 3rd
           Hom_0_3rd = Hom_0_3rd * iHom_0;
           penalty += getPenalty ( id, proposal, Hom_0_3rd, Kl2 * centerB, solution3rd, &normals3rd, centers_, dummy);
        };

        _heapEntryStore.push_back ( HeapStoreEntry( penalty, id, proposal, i) );// why not a second criterium here; like distance or so
        _heap.insert( HeapEntryG (_heapEntryStore.size()-1) );
      }
    }

    int dead = 0;// debug mode:
    while ( !_heap.empty() ) 
    {
      HeapEntryG hid = _heap.front();
      _heap.pop_front();
      int sId = hid.pId(); // store id
      HeapStoreEntry& heapE = _heapEntryStore[ sId ];

//      if (heapE.cost() >= 9999.) break;// just randomness

      int proposal = heapE.pixNr();
      int     seg2 = heapE.segNr(); 
      int     seg1 = heapE.segNr_left(); 

      if ( scores[seg2] > heapE.cost())
//      if ( (solution2[seg2] == -1) || ( scores[seg2] > heapE.cost() && ( heapE.cost() == 0 || scores[seg2] >= __impPenalty__ ) ) )
      {
        scores[seg2] = heapE.cost();
        solution2[seg2] = proposal;
      }

      if (solution2[seg2] == -1)
      {
        scores[seg2] = heapE.cost();
        solution2[seg2] = proposal;
      }

#ifdef _centeredViews_
      M3 Hom_0  = gHom->getHomC0( proposal );// getHom is also ok if normals are not negative
#else
      M3 Hom_0  = gHom->getHom( proposal, seg1 );
#endif
      M3 iHom_0 = Hom_0;iHom_0.invert();

      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges2_, seg2) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges2_, seg2) );

      for (int j=0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        // not ideal to stop it here - but no other choice, would run inf else
//        if (solution2[id] != -1) continue; // this would be problematic to leave out
        // else push into Q
        // NEW PUSH ONLY IF LOWER SCORE

        if (dead)
          return;

        // project back
        P3 centerB = P3( &centers2_[ 3*id ] );

        // now get a score when projecting back from 2nd view:
        Scalar penalty = getPenalty ( id, proposal, iHom_0, Kl2 * centerB, solution, normals, centers_, seg1);

        if (gHom2nd != NULL)
        {
          int dummy=0;
           M3 Hom_0_2nd  = gHom2nd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 2nd
           // i need from solution filled view to 2nd:
           // to first then to 2nd
           Hom_0_2nd = Hom_0_2nd * iHom_0;
           penalty += getPenalty ( id, proposal, Hom_0_2nd, Kl2 * centerB, solution2nd, &normals2nd, centers_, dummy);
        };
        if (gHom3rd != NULL)
        {
          int dummy=0;
           M3 Hom_0_3rd  = gHom3rd->getHomC0( proposal );// getHom is also ok if normals are not negative
           // from canonical to 3rd
           // i need from solution filled view to 3rd:
           // to first then to 3rd
           Hom_0_3rd = Hom_0_3rd * iHom_0;
           penalty += getPenalty ( id, proposal, Hom_0_3rd, Kl2 * centerB, solution3rd, &normals3rd, centers_, dummy);
        };

         if (scores[id] > penalty ) // must be identical to condition above
//         if ( scores[id] > penalty && ( penalty == 0 || scores[id] >= __impPenalty__ ) ) 
         {
          _heapEntryStore.push_back ( HeapStoreEntry( penalty, id, proposal, seg1) );// why not a second criterium here; like distance or so
          _heap.insert( HeapEntryG (_heapEntryStore.size()-1) );
         }
      }
    }
    // finally just grow. that oob area should just be grown. ohterwise 9999 delivers random stuff or worse.
  }

template<class Scalar>
Scalar OcclusionMapBufferVC<Scalar>::
getPenalty ( int ownSegment, int proposal, const M3& iHom_0, const P3& centre, const std::vector<int>& solution, const std::vector<P3>* normals, Scalar *centers_, int& seg1)
  {
    Scalar impPenalty = 250.;
    Scalar oobScore   = 100.;// == occpen 
    Scalar missSeg    = 100.; // < occPen

    int area =1; // full segment
    Scalar freePix = 1;// no clue - all

    P3 p0     = iHom_0 * centre;// center is in camera coords
    p0 /= p0[2];// in pixel first image.

    int inBounds = (int( p0[0]-0.5 ) >= 0) && (int( p0[0]-0.5 ) < w) && (int( p0[1]-0.5 ) >= 0) && (int( p0[1]-0.5 ) < h);
    // 1. oob check - if oob oob penalty here
    if( ! inBounds )
      return oobScore;

    // proceed : find case occ oob, wrong segi
    int pixId = int( p0[0]-0.5 )*h +int ( p0[1]-0.5 );
    seg1   = segImg [pixId];

    if( solution[ seg1 ] == proposal )
      return 0.;

    //
    P3 normal2Seg1   = (*normals)[proposal];// normal by projected plane
    Scalar depth0    = 1/(normal2Seg1 | P3( &centers_[ 3*seg1 ] ));//?
    P3 normal1Seg1   = (*normals)[ solution[ seg1 ] ];//stored in orig view
    Scalar depth1    = 1/(normal1Seg1 | P3( &centers_[ 3*seg1 ]) );

    Scalar __OCC_THRESH__ = 0.1;
    Scalar maxDepthThresh = 0.1;
    Scalar _minDepth_     = 1.0;// visible - always true if canonical view is involved

    // here check depth occ imp or ok'ish 
    Scalar depthThresh = max ( maxDepthThresh, fabs(depth0  * __OCC_THRESH__) );
    // let depths by positive -> depthDiff>0 implies an occlusion case, sign change: < 0 occlusion case
    Scalar depthDiff   = depth0 - depth1;
    if ( depth0 > Scalar(-_minDepth_) )
      return Scalar(10)*impPenalty;// not possible 

    // occ:
    if ( depthDiff < -depthThresh )
      return oobScore;

    if ( depthDiff > depthThresh )
      return impPenalty;

    return missSeg;
  }


#endif