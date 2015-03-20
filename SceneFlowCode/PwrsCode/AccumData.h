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

#ifndef _ACCUM_DATA_H
#define _ACCUM_DATA_H

#include <algorithm>
#include "genWarp.h"
#include "AccumDataWarp.h"
#include "OobInfo.h"

/// accumulate data scores for all image pairs
template<typename Scalar> class accumulate2FrameData
{
public:

  typedef Math::Mat3x3T< Scalar>     M3;
  typedef Math::VectorT< Scalar, 3>  P3;
  typedef Math::VectorT<int, 4>      P4i;


  accumulate2FrameData( int N_, int M_ )
    : nSegments(0), nNormals(0), N(N_), M(M_), allScores(), allPartialScores(), freeAllPartialScores(), freeVariables(), 
      doStoreLocaldataScores(false), localdataScores()
  {};

  ~accumulate2FrameData(){};

  void addWarp( genWarp<Scalar> *gw_ )                  {warpVec.push_back( gw_ );};
  void addAccumulateWarp( accumulateWarp<Scalar> *aw_ ) {aWarpVec.push_back( aw_ );};

  /// special class which stores the warped coordinates in a vector idx, idy, one for each image
  void addOcclusionMappingBuffer( OcclusionMapBuffer<Scalar>* oMapBuf_) {oMapBufVec.push_back( oMapBuf_ );}

  std::vector<Scalar>& getAllScores() {return allScores;}; 
  std::vector<Scalar>& getPartialScores(int i=0) {return allPartialScores[i];}; 

  /// special case if boundary pixel are nailed to a fixed value
  std::vector< std::vector<Scalar> >& getFreePartialScores() {return freeAllPartialScores;};
  /// return the amount of free variables currently available
  std::vector< std::vector<int> >&    getNFreeVariables() {return freeVariables;};

  // for initialisation
  Scalar getBestConfiguration(int nNormals, std::vector<int>& bestDataConfig )
  {
    // init all 0
    bestDataConfig.resize( nSegments, 0);
    if (allScores.size() < 1 && allPartialScores.size() < 1)
      return 0;

    Scalar sum (0);
    
    if (allScores.size() > 1)
    {
      for (int i = 0; i < nSegments ;i++)
      {
        Scalar bestScore( allScores[i] );
        int bestId = 0;
        for (int n = 1; n < nNormals; n++)
          if ( bestScore > allScores[n*nSegments+i] )
          {
            bestId = n;
            bestScore=allScores[n*nSegments+i];
          }

          bestDataConfig[i] = bestId;
          sum += bestScore;
      }
    }
    else
    {
      int nPartialScores =  (int)(allPartialScores.size());

      for (int i = 0; i < nSegments ;i++)
      {
        Scalar locsum(0);
        for (int np =0;np < nPartialScores; np++)
          locsum += (allPartialScores[np])[i];
        
        Scalar bestScore( locsum );
        int bestId = 0;
        for (int n = 1; n < nNormals; n++)
        {
          Scalar locsum(0);
          for (int np =0;np < nPartialScores; np++)
            locsum += (allPartialScores[np])[n*nSegments+i];
          if ( bestScore > locsum )
          {
            bestId = n;
            bestScore = locsum;
          }
        }
        bestDataConfig[i] = bestId;
        sum += bestScore;
      }
    }
    return sum;
  }

    void computeVecScores( const mxArray*& Segments, int nLabels, int* segImg, 
      const std::vector<int>& SegIds,
      const std::vector<bool>& oobInfoDepth, const std::vector<bool>& oobInfoFlow, 
      const std::vector<bool>& oobInfoFlowR, const std::vector<bool>& oobInfoDepthT, 
      const std::vector<bool>& oobInfoCrossLRT )
  {
    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
    allScores.clear();
    allScores.resize( nSegments, 0 );
#ifdef _occHandlingOn_
    freeAllPartialScores.clear();
    freeAllPartialScores.resize( aWarpVec.size() );
    freeVariables.clear();
    freeVariables.resize( aWarpVec.size() );
#endif

    if ( doStoreLocaldataScores )
    {
      localdataScores.clear();
      localdataScores.resize( aWarpVec.size() );
    }

    assert(SegIds.size() >= nSegments);

#pragma omp parallel for schedule (static)
      for (int i=0; i < warpVec.size(); i++) // warpVec.size()-1
        (warpVec[i])->warp_noOmp_patchBased( Segments, SegIds );

#pragma omp parallel for schedule (static)
      for (int i=0; i < aWarpVec.size(); i++)
      {
        switch (i)
        {
        case 0:
          aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, 0,0,N,M, oobInfoDepth ); // l r
          break;
        case 1:
          aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, 0,0,N,M, oobInfoFlow ); // l lt
          break;
        case 2:
          aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, 0,0,N,M, oobInfoFlowR ); // rt + rt1
          break;
        case 3:
          aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, 0,0,N,M, oobInfoDepthT ); // lt1 rt1
          break;
        case 4:
          aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, 0,0,N,M, oobInfoCrossLRT ); // lt1 rt1
          break;
        default:
          printf("Something does not add up - timeout\n");
          break;
        }
      }

      for (unsigned int i=0; i < aWarpVec.size(); i++)
      {
        std::transform((aWarpVec[i]->getScores()).begin(), (aWarpVec[i]->getScores()).end(), allScores.begin(), allScores.begin(), plus<Scalar>());
#ifdef _occHandlingOn_
        freeAllPartialScores[i].resize( aWarpVec[i]->getFreeScores().size() );
        freeVariables[i].resize( aWarpVec[i]->getFreeScores().size() );
        std::copy((aWarpVec[i]->getFreeScores()).begin(), (aWarpVec[i]->getFreeScores()).end(), freeAllPartialScores[i].begin() );
        std::copy((aWarpVec[i]->getFreeVariables()).begin(), (aWarpVec[i]->getFreeVariables()).end(), freeVariables[i].begin() );
#endif
        // sometime used for jittering ??
        if ( doStoreLocaldataScores )
        {
          localdataScores[i].resize( aWarpVec[i]->getFreeScores().size() );
          std::copy((aWarpVec[i]->getFreeScores()).begin(), (aWarpVec[i]->getFreeScores()).end(), localdataScores[i].begin() );
        }
      }
    };

    std::vector<Scalar>& getDataScores(int i=0) {return localdataScores[i];}; 

    void storeLocaldataScores( bool onOff=true )
    {
      doStoreLocaldataScores = onOff;
      if ( !onOff )
      {
        localdataScores.clear();
      }
    }

    /// compute per bounding boxes
    void computeAllVecScores_boxed( const mxArray*& Segments, int nLabels, int* segImg, 
      std::vector< std::list< std::pair<int,int> > >&   pageSegPairs, 
      std::vector< P4i >& bboxes,
      int propList, // proposals to evaluate just a number of how many
      const std::vector<int>& AllSegIds,
      const std::vector<bool>& oobInfoDepth, const std::vector<bool>& oobInfoFlow, 
      const std::vector<bool>& oobInfoFlowR, const std::vector<bool>& oobInfoDepthT, 
      const std::vector<bool>& oobInfoCrossLRT )
    {
      nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
      allScores.clear();
      allScores.resize( AllSegIds.size(), 0 );
#ifdef _occHandlingOn_
      freeAllPartialScores.clear();
      freeAllPartialScores.resize( aWarpVec.size() );
      freeVariables.clear();
      freeVariables.resize( aWarpVec.size() );

      for (int i =0; i<aWarpVec.size(); i++ )
      {
        freeAllPartialScores[i].resize( AllSegIds.size(), 0 );
        freeVariables[i].resize( AllSegIds.size(), 0 );
      }
#endif
      assert(AllSegIds.size() >= nSegments);

      int iterations = propList;

      int nCovers = AllSegIds.size() / propList;
      printf("n covers: %d\n", nCovers);

      std::vector<int> SegIds (AllSegIds.begin(), AllSegIds.begin()+nSegments );
      // init once if used more often:
      (warpVec[warpVec.size()-1])->warp_noOmp_patchBased( Segments, SegIds );

      for (int its=0;its < iterations; its++)
      {
        int proposal = its;
        P4i bbox = bboxes[proposal % nSegments];

#pragma omp parallel for schedule (static)
        for (long int i=0; i < warpVec.size()-1; i++) // warpVec.size()-1
          (warpVec[i])->warp_noOmp_patchBased( Segments, proposal, bbox[0], bbox[1], bbox[2], bbox[3], segImg );

#pragma omp parallel for schedule (static)
        for (int i=0; i < aWarpVec.size(); i++)
        {
          switch (i)
          {
          case 0:
            aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, bbox[0], bbox[1], bbox[2], bbox[3], oobInfoDepth ); // l r
            break;
          case 1:
            aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, bbox[0], bbox[1], bbox[2], bbox[3], oobInfoFlow ); // l lt
            break;
          case 2:
            aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, bbox[0], bbox[1], bbox[2], bbox[3], oobInfoFlowR ); // rt + rt1
            break;
          case 3:
            aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, bbox[0], bbox[1], bbox[2], bbox[3], oobInfoDepthT ); // lt1 rt1
            break;
          case 4:
            aWarpVec[i]->computeFullScoresCensus3OMP_Box( Segments, segImg, bbox[0], bbox[1], bbox[2], bbox[3], oobInfoCrossLRT ); // l rt1
            break;
          default:
            printf("Something does not add up - timeout\n");
            break;
          }
        }


        std::list< std::pair<int,int> > plist = pageSegPairs[proposal % nSegments];
        std::list< std::pair<int,int> >::iterator p_it = plist.begin(), p_end = plist.end();

        int position = int(proposal / nSegments) * nCovers * nSegments;//pageSegPairs.size();

        for( ; p_it != p_end; p_it++)
          for (unsigned int i=0; i < aWarpVec.size(); i++)
          {
            allScores[position+p_it->first*nSegments + p_it->second] += (aWarpVec[i]->getScores())[p_it->second];
#ifdef _occHandlingOn_
            freeAllPartialScores[i][position+p_it->first*nSegments + p_it->second] = (aWarpVec[i]->getFreeScores())[p_it->second];
            freeVariables[i][position+p_it->first*nSegments + p_it->second] = (aWarpVec[i]->getFreeVariables())[p_it->second];
#endif
          }
      }
    };

  // this is new: take care of still need the box here for 8 neigh check, box ids (or M) define the memory addresses
  /// ids are the GLOBAL ids in order to grab the correct idx warped indices (getIdx, getSegIdx), edgeids could be the local (-1:no vertex, x>0 mrf id =x)
  void computeScorePairWindowCensusNew_short( 
    std::vector< int >* ids, std::vector< int >* localIds, int Mstep,//std::vector< int >&otherEdgeIds
    const std::vector<bool>& oobInfoDepth, const std::vector<bool>& oobInfoFlow, 
    const std::vector<bool>& oobInfoFlowR, const std::vector<bool>& oobInfoDepthT, bool doOldIds = false )
  {
    int nPixel(0);
    allScores.clear( );

    std::vector< int > trueIds;
    std::vector< int > truelocalIds;
    bool trialSegment = true;
    if ( localIds == NULL )
    {
      trueIds.resize(N*M);
      truelocalIds.resize(N*M);
      for( int i=0;i<N*M;i++ )
        trueIds[i] = i; // or other way round - no idea
      ids = &trueIds;
      std::copy(trueIds.begin(), trueIds.end(), truelocalIds.begin() );
      localIds = &truelocalIds;
      Mstep = M;
      trialSegment = false;
    }

    nPixel = ids->size();

    allPartialScores.resize( aWarpVec.size() ); // for each image pair
    for( int i=0;i<aWarpVec.size();i++ )
      (allPartialScores[i]).resize(2*nPixel, 0);

#pragma omp parallel for schedule (static)
      for (int i=0; i < warpVec.size(); i++)
      {
        std::vector< Scalar > idx0;
        std::vector< Scalar > idy0;

        // ids used to return only those idx which are free variables - those listed in ids
        // which are those in the globids container!
        if (!trialSegment)
        {
          oMapBufVec[i]->getFullIdx( idx0 );
          oMapBufVec[i]->getFullIdy( idy0 );
        }
        else
        {
          if (!doOldIds)
          {
          oMapBufVec[i]->getSegIdx( idx0, *ids );
          oMapBufVec[i]->getSegIdy( idy0, *ids );
          }
          else
          {
          oMapBufVec[i]->getIdx( idx0, *ids );
          oMapBufVec[i]->getIdy( idy0, *ids );          
          }
        }
      
        /// produces two images instead of one here OMG and TWO ids and ... !!!
        (warpVec[i])->warp_noOmp_patchBased( idx0, idy0);//, idx1, idy1 );
      }

#pragma omp parallel for schedule (static)
      for (int i=0; i < aWarpVec.size(); i++)
      {
        switch (i)
        {
        case 0:
          aWarpVec[i]->computeFullScoresMappedCensusNew_short( *ids, *localIds, Mstep, oobInfoDepth); // l r , &otherEdgeIds
          break;
        case 1:
          aWarpVec[i]->computeFullScoresMappedCensusNew_short( *ids, *localIds, Mstep, oobInfoFlow); // l lt , &otherEdgeIds
          break;
        case 2:
          aWarpVec[i]->computeFullScoresMappedCensusNew_short( *ids, *localIds, Mstep, oobInfoFlowR); // rt + rt1 , &otherEdgeIds
          break;
        case 3:
          aWarpVec[i]->computeFullScoresMappedCensusNew_short( *ids, *localIds, Mstep, oobInfoDepthT); // lt1 rt1  , &otherEdgeIds
          break;
        default:
          printf("Something does not add up - timeout\n");
          break;
        }

        std::vector<Scalar> &tmp = allPartialScores[i];
        std::copy( (aWarpVec[i]->getScores()).begin(), (aWarpVec[i]->getScores()).end(), tmp.begin() + 0*nPixel ); // first scores data
        std::copy( (aWarpVec[i]->getScores2()).begin(), (aWarpVec[i]->getScores2()).end(), tmp.begin() + 1*nPixel );// 2nd scores penalties
      }
  };

private:

  int nSegments;
  int nNormals;
  int N;
  int M;

  std::vector<Scalar> allScores;
  std::vector< std::vector<Scalar> > allPartialScores;

  /// unused (currently)
  std::vector< std::vector<Scalar> > localdataScores;
  /// store data stores per image pair
  bool doStoreLocaldataScores;

  /// some pixel are fixed and contribute by penalties if oob or not 
  std::vector< std::vector<Scalar> > freeAllPartialScores;
  /// the variables which are given a data term, so these are to be replcaed if there is an occlusion
  std::vector< std::vector<int> > freeVariables;

  std::vector< genWarp<Scalar>* >        warpVec;
  std::vector< accumulateWarp<Scalar>* > aWarpVec;

  std::vector< OcclusionMapBuffer<Scalar>* >  oMapBufVec;
};
#endif