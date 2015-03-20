/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
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

//////////////////////////////////////////////////////////////////////////////
// Manage the partial instantion of the pixel graph in the pixel refinement //
//////////////////////////////////////////////////////////////////////////////

#ifndef __PatchSegmentHandler__cpp
#define __PatchSegmentHandler__cpp

#include "PatchSegmentHandler.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <algorithm>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

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

  /// prepare energy computation for a patch, new globIds for all! pixel in region, locids only for inner and non segment pixel used in smoothing
template<class Scalar>
void
PatchIdGenerator<Scalar>::
initMappingsNew( int patch, bool doEdgeMapping )
  {
    int cSize = cRad-1;

    if (patch >= 0)
    {
      startX = max(0,   (centers[patch])[0]-patchSize);
      endX   = min(w-1, (centers[patch])[0]+patchSize)+1;
      startY = max(0,   (centers[patch])[1]-patchSize);
      endY   = min(h-1, (centers[patch])[1]+patchSize)+1;
      currentSeg   = startSegs[patch];
      currentPatch = patch;
    }
    else // full image:
    {
      startX = 0;
      endX   = w;
      startY = 0;
      endY   = h;
      currentSeg   = startSegs[0];
      currentPatch = 0;
    }


    loc2glob.clear();
    loc2glob.resize( (endX-startX) * (endY-startY) );
    locIds.resize(   (endX-startX) * (endY-startY), -1 );
    globIds.clear();  // only globs of mrf ids
    allGlobIds.clear();//all the glob ids in the window

    freeFixedIds.clear();
    fixedIds.clear();
    if (doCensus)
       fixedIds.resize( (endX-startX) * (endY-startY), false );

    int locId = 0;
    int locPos = 0;
    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++, i++, locPos++)
      {
        glob2loc[i]      = locPos;// do not need that
        loc2glob[locPos] = i;
        allGlobIds.push_back( i ); // all for data term - could be reduced but this is not nice

        int theLocalId = getLocId( i ); // should be equal to locPos

        // do not fix pixel at image borders no matter what!
        if (segImg[i] == currentSeg)
        {
            locIds[locPos] = -1;
//          loc2glob[locPos] = -1;
        }
        else
        {
          if ( (ii<=startX+cSize && startX >cSize) || (ii>=endX-(cSize+1) && endX<(w-cSize)) || 
               (pos <= (startY+cSize)*w && startY >cSize) || (pos >= (endY-(cSize+1))*w && endY < (h-cSize)) )
          {
            fixedIds[locPos] = true;
            locIds[locPos] = -1;
          }
          else
          {
            globIds.push_back( i ); // global ids for all mrf pixel this is getting ugly
            locIds[locPos] = locId;
            locId++;
          }
          // if it is on the boundary fix it, if the boundary is not at image borders
        }
      }
    }

    mrfIds = locId;

    if (doCensus && doEdgeMapping)
    {
#ifdef __USE_FIXED_BORDERS__
      freeFixedIds.clear();
      freeFixedIds.reserve(globIds.size());
      for(int i=0;i< locIds.size();i++)
        if( locIds[i] >=0 )
          freeFixedIds.push_back( fixedIds[i] );
#else
      freeFixedIds.clear();
      freeFixedIds.reserve(globIds.size(), false);
#endif      
      initEdgeMappingNew();
    }
  }
  
  /// return the number of new proposals (old ones and new ones) potentially increase expansion area instead? - also keep original positions(segment centers)
template<class Scalar>
int
PatchIdGenerator<Scalar>::
refreshMappings( bool onlyNew )
{
    // need to update segmentation as well - if more than 2 steps are done
    // 0. which segments are still present ? mark the conquered ones and remove 
    // 1. new centers: center of expansion proposal, startSegs: segment of expansion proposal
    // 2. startSegs
    
    int border = 2;

    printf("refreshMappings, started with %d segments\n", startSegs.size());

    int maxSeg(0);
    for(int i=0;i<startSegs.size();i++)
      maxSeg = std::max(startSegs[i], maxSeg);

    std::vector<bool> presentS(2*maxSeg, false);

    for(int i=0;i<w*h;i++)
      presentS[segImg[i]] = true;

    // now expressed in startSeg positions:
    std::vector<bool> present(startSegs.size(), false);

    for(int i=0;i<startSegs.size();i++)
      if (presentS[ startSegs[i]] )
      present[i] = true;


    int num(0);
    for(int i=0;i<present.size();i++)
          if (present[i]) num++;
     printf("refreshMappings, presents %d \n", num);

    int last( present.size()-1 ), first( 0 ), lastinValid( present.size()-1 );
    while( last > first )
    {
      while( first <present.size()-1  && present[first] )
        first++;
      while( last > 0 && !present[last] )
        last--;

      // swap
      if( !present[first] && present[last] && last > first )
      {
        P2i tempP      = centers[first];
        centers[first] = centers[last];
        centers[last]  = tempP;

        int tempS        = startSegs[first];
        startSegs[first] = startSegs[last];
        startSegs[last]  = tempS;
        present[first]   = true;
        present[last ]   = false;
        lastinValid      = last;
      }
    }
    startSegs.resize(lastinValid);
    centers.resize(lastinValid);

    printf("refreshMappings, after removal %d segments\n", startSegs.size());

    // insert new centers: how to find these positions
    //centers[patch];
    //startSegs[patch];

    // how can i get the fully expanded positions? grow from the start position find pixel which are far away (Manhattan distance)
    // these can be several for a single segment
    // now pick pixel and set a startSegment,center
    // remove all pixel covered by these new pixel/expansion areas
    // repeat until all is covered

    // simple approach: check all pix in area - also could use Q

    std::vector<bool> coveredpixel(w*h, false);
    for(int patch=0; patch< startSegs.size(); patch++)
    {
      int seg = startSegs[patch];

      int startX = max(0,   (centers[patch])[0]-patchSize+border);
      int endX   = min(w-1, (centers[patch])[0]+patchSize-border);
      int startY = max(0,   (centers[patch])[1]-patchSize+border);
      int endY   = min(h-1, (centers[patch])[1]+patchSize-border);

      for (int j=startY;j<=endY;j++)
      {
        int pos = j*w;int i = pos+startX;
        for (int ii=startX;ii<=endX;ii++, i++)
        {
          if ( seg == segImg[i])
            coveredpixel[i] = true;
        }
      }
    }
    /////
    // gather all uncovered pixel
    std::vector<int> uncovered;
    uncovered.reserve(10000);
    for(int i=0; i< h*w; i++)
      if (!coveredpixel[i])
        uncovered.push_back(i);

    random_shuffle( uncovered.begin(), uncovered.end() );

    printf("# uncovered pixel : %d\n", uncovered.size());

    int newPatches = 0;
    for (int uc = 0; uc<uncovered.size();uc++)
    {
      int expCenter = uncovered[uc];
      if (coveredpixel[ expCenter ])
         continue;
      // else : insert as a expansion proposal and cover the other pixel
      int seg = segImg[expCenter];
      startSegs.push_back( seg );

      int qx = expCenter % w;
      int qy = expCenter / w;
      centers.push_back( P2i( qx, qy ) );
      newPatches++;
//      printf("pos uncovered pixel : %d, %d\n", qx,qy);

      int startX = max(0,   qx-patchSize+border);
      int endX   = min(w-1, qx+patchSize-border);
      int startY = max(0,   qy-patchSize+border);
      int endY   = min(h-1, qy+patchSize-border);

      for (int j=startY;j<=endY;j++)
      {
        int pos = j*w;int i = pos+startX;
        for (int ii=startX;ii<=endX;ii++, i++)
        {
          if ( seg == segImg[i])
            coveredpixel[i] = true;
        }
      }
    }

    // idea is to apply expansion only on patches of interest by putting these in front of the vector
    if (onlyNew)
    {
      for(int i=0;i<newPatches;i++)
      {
        P2i tempP               = centers[i];
        centers[i]              = centers[i+lastinValid];
        centers[i+lastinValid]  = tempP;

        int tempS                 = startSegs[i];
        startSegs[i]              = startSegs[i+lastinValid];
        startSegs[i+lastinValid]  = tempS;
      }
      nInitSegments = newPatches;
      //nInitSegments = startSegs.size()-lastinValid;
    }
    else
    {
      nInitSegments = startSegs.size();
    }

    printf("refreshMappings, after insertion %d segments\n", startSegs.size());
    return startSegs.size();
  }

template<class Scalar>
void
PatchIdGenerator<Scalar>::
setCensus(int cSize)
  {
    int censusSize = ((2*cSize+1)*(2*cSize+1)-1)/2;
    cRad = cSize;
    dxC.clear();
    dyC.clear();
    dxyC.clear();
    dxC.resize(  censusSize );
    dyC.resize(  censusSize );
    dxyC.resize( censusSize  );

    int sumAll=0;
    for (int xx = -cSize; xx <= cSize; xx++ )
      for (int yy = -cSize; yy <= cSize; yy++ )
        if (sumAll < censusSize)
        {
          dxC[sumAll] = xx;
          dyC[sumAll] = yy;
          dxyC[sumAll] = w*xx+yy;
          sumAll++;
        }
  }

/// edgeIds are only used for data term, so all edges in window are concerned here so getLocId not getMrfId
template<class Scalar>
void
PatchIdGenerator<Scalar>::
initEdgeMappingNew( )
  {
    int cSize = dxyC.size();

    edgeIds.clear();
    edgeIds.resize(cSize * allGlobIds.size(), -1);

    otherEdgeIds.clear();
    otherEdgeIds.resize(cSize * allGlobIds.size(), -1);


    for (int i=0; i< allGlobIds.size(); i++ )
    {
      int gId = allGlobIds[i];

      for(int j = 0; j< cSize;j++)
      {
        // left up -1
        int nid  = getLocId(gId + dxyC[j]);
        edgeIds[cSize*i+j] = nid;

        nid  = getLocId(gId - dxyC[j]);
        otherEdgeIds[cSize*i+j] = nid;
      }
    }
  }


  /// writes the penalty scores stored in dataScore (in the buffer) to the vector to be as data score for a window to be optimized
  /// combines external score from afd with internally buffered score
template<class Scalar>
void
unaryDataScoresBuffer<Scalar>::
constructDataScore( std::vector<Scalar>& newScores, const std::vector<Scalar>& trialScores, std::vector<int>& loc2Glob, PatchIdGenerator<Scalar>& genPatchIds )
  {
    int nPixel = loc2Glob.size();

    int nbufferScores = dataScore.size()/2;
    int nDataScores   = trialScores.size()/2;
    // assert that newScores is of the right size, we fill every second entry here with the stored values
    // can fill with -1(no edge), 0 no data penalty or 1 data penalty

    int nVars = genPatchIds.getNMrfVars();
    // order is data penalty @ mrf ==0 and data,penalty @ mrf ==1
    newScores.resize( 4*nVars );

    int validPixNr(0);
    for (int i=0;i<loc2Glob.size();i++)
    {
      int gId = loc2Glob[i];
      int nid  = genPatchIds.getMrfId(gId);

      // left up -1
      if( nid>=0 )// must be
      { 
        validPixNr++;
//        newScores[validPixN] = penaltyScore[gId]; // new scores for 0 proposal

        newScores[nid]         = dataScore[gId   ]; // new scores for 0 proposal
        newScores[nid+nVars]   = dataScore[gId + nbufferScores]; // the data penalties attached after the data scores

        int pid = genPatchIds.getLocId( gId );
        int global2Local = genPatchIds.getLocIdFromGlobal(gId);

        // need a map from the global or local mrf coords to the window coords
        newScores[nid+2*nVars] = trialScores[pid   ]; // new scores for 0 proposal
        newScores[nid+3*nVars] = trialScores[pid + nDataScores]; // the data penalties attached after the data scores

        int test1 = trialScores[pid   ];
        int test2 = trialScores[i     ];
      }
        /// map local to global ids: ok, but find also neighbor value: look at 01, 10, 11 - all cases are possible
    }
  }

  //////////////////////////////////////
  /// update the data score, based on a score map and the changes: 0 nothing changed 1 changed segment map
template<class Scalar>
void
unaryDataScoresBuffer<Scalar>::
updateDataScore( const std::vector<Scalar>& trialScores, const std::vector<int>& loc2Glob, 
                 const PatchIdGenerator<Scalar>& genPatchIds, const std::vector<int>& changes )
  {
    int nbufferScores = dataScore.size()/2;
    int nDataScores   = trialScores.size()/2;

    int nPixel = loc2Glob.size();
    for (int i=0;i<loc2Glob.size();i++)
    {
      if ( changes[i]<=0 ) continue;
//      if ( changes[nid]<=0 ) continue;

      int gId = loc2Glob[i];
      int nid = genPatchIds.getMrfId( gId );

      int pid = genPatchIds.getLocId( gId );

      int global2Local = genPatchIds.getLocIdFromGlobal(gId);

      dataScore[gId]               = trialScores[ pid ];
      dataScore[gId+nbufferScores] = trialScores[ pid + nDataScores];
      /// map local to global ids: ok, but find also neighbor value: look at 01, 10, 11 - all cases are possible
    }
  }


#endif // __PatchSegmentHandler__h
