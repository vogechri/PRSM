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

/// idea is that this function just need a sceond proposal for each segment 
/// it applies thes only in this exact spot
/// allowing gradient descent like algorithms
#ifndef __QPBO_OCC_2FRAME__
#define __QPBO_OCC_2FRAME__
////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////
#define _occHandlingOn_

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "DataDefinitions.h"
#include "AccumData.h" //#include "GenerateWarp_logic.h"
#include "EvalEnergyFull2Frame.h"

#include "EvalEnergyStore.h"

#include "OcclusionMappingBufferSegments.h"
#include "BinaryConversion.h"

#include "LocalProposalHandler.h"

#include <algorithm>
#include <list>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T
#include "../QPBO-v1.3.src/QPBO.h"//#include "QPBO.h"

#include <ctime>

using namespace std;
using namespace Math;

//#define __Small_Improvement__ 0001
#define __Small_Improvement__ 0.501

// also use cross view
//#define _use_5th_view

//#define __DEBUG_ENERGIES__

  //////////////////////////////////////////////////////////////////////////////////////

template<typename Scalar>
void getEnergyVectors ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda, std::vector<Scalar>& solDataScores, std::vector<Scalar>& solSmoothScores)
{
  solSmoothScores.clear();
  solSmoothScores.resize(nSegments, 0);
  solDataScores.clear();
  solDataScores.resize(nSegments, 0);
  /////
  for ( int j = 0; j< nSegments; j++ )
    solDataScores[j] = dataScores[j];

  evalSmooth.compute_score_Fuse( current, current);
//  evalSmooth.compute_score_combiDepth( 0, current);//, oobT, oobS );
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  const std::vector<int>& idk=evalSmooth.getIdk();
  const std::vector<int>& idl=evalSmooth.getIdl();

  for ( int i = 0; i < f00.size(); i++)
  {
    solSmoothScores[idk[i]] += f00[i] * lambda;
    solSmoothScores[idl[i]] += f00[i] * lambda;
  }
}


template<typename Scalar>
Scalar getEnergy ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda)//, std::vector<bool> oobS,  std::vector<bool> oobT)
{
  Scalar score(0.);
  for ( int j = 0; j< nSegments; j++ )
    score += dataScores[j];

  evalSmooth.compute_score_Fuse( current, current);
//  evalSmooth.compute_score_combiDepth( 0, current);//, oobT, oobS );
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  for ( int i = 0; i < f00.size(); i++)
    score+= f00[i] * lambda;

  return score;
}

template<typename Scalar>
Scalar getDataEnergy ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores )
{
  Scalar score(0.);
  for ( int j = 0; j< nSegments; j++ )
    score += dataScores[j];

  return score;
}

template<typename Scalar>
Scalar getSmoothEnergy ( std::vector<int>& current, int nSegments, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda)//, std::vector<bool> oobS, std::vector<bool> oobT)
{
  Scalar score(0.);

  evalSmooth.compute_score_Fuse( current, current);
//  evalSmooth.compute_score_combiDepth( 0, current);//, oobT, oobS );
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  for ( int i = 0; i < f00.size(); i++)
    score+= f00[i] * lambda;

  return score;
}

template<typename Scalar>
Scalar getOcclusionEnergy ( std::vector<int>& currentSolution, Scalar occThresh, 
                            OcclusionMapBufferSeg<Scalar>& ombs1, OcclusionMapBufferSeg<Scalar>& ombs2, OcclusionMapBufferSeg<Scalar>& ombs3,
                            std::vector< std::vector<Scalar> >& freeScoresCurrent, std::vector< std::vector<int> >& freeCurrent)
{
  Scalar occPenalty = occThresh;//0.75; // best thing is to do nothing !
        /////////////////////////////////////
        // now add the occlusion edges 
        ombs1.buildList(currentSolution, currentSolution );
        ombs1.generateOcclusionLists();
        ombs2.buildList(currentSolution, currentSolution );
        ombs2.generateOcclusionLists();
        ombs3.buildList(currentSolution, currentSolution );
        ombs3.generateOcclusionLists();

        std::vector< occSegList<Scalar> >& _occList  = ombs1.getOccList();
        std::vector< occSegList<Scalar> >& _occList2 = ombs2.getOccList();
        std::vector< occSegList<Scalar> >& _occList3 = ombs3.getOccList();

        BinaryConverter<Scalar> bc_0( occPenalty, _occList , freeCurrent[0], freeCurrent[0] );
        BinaryConverter<Scalar> bc_1( occPenalty, _occList2, freeCurrent[1], freeCurrent[1] );
        BinaryConverter<Scalar> bc_2( occPenalty, _occList, _occList3,  freeCurrent[2], freeCurrent[2] );
        BinaryConverter<Scalar> bc_3( occPenalty, _occList2, _occList3, freeCurrent[3], freeCurrent[3] );

        bc_0.setPenaltiesPerSegment( freeScoresCurrent[0], freeScoresCurrent[0] );
        bc_1.setPenaltiesPerSegment( freeScoresCurrent[1], freeScoresCurrent[1] );
        bc_2.setPenaltiesPerSegment( freeScoresCurrent[2], freeScoresCurrent[2] );
        bc_3.setPenaltiesPerSegment( freeScoresCurrent[3], freeScoresCurrent[3] );

        Scalar penalty(0);
        ////////////////////////
        penalty += bc_0.generateUnaries0L0();
        penalty += bc_1.generateUnaries0L0();
        penalty += bc_2.generateUnaries0L0();
        penalty += bc_3.generateUnaries0L0();

        return penalty;
}

// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void evalQPBO2FramesNew( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

  size_t N,M;
  size_t elements;
  char buffer[256];

  Scalar thresh = 0.1;
  Scalar lambda = 0.02;
  Scalar theta  = 100.0;
  Scalar dispMax = 32.;
  Scalar rotJump(15.);//   setRotJump (rotJump)
  Scalar depthJump(sqrt(1300.));//   setRotJump (rotJump)
  Scalar occThresh(0.05/0.75);
  Scalar motMax = 220; // the maximal motion in the scene

  Scalar* halfEdgeX(NULL);
  Scalar* halfEdgeY(NULL);
  Scalar* halfEdgeXY(NULL);
  Scalar* halfEdgeiXY(NULL);
  int* origSegImg(NULL);

  Scalar* Kr_(NULL);
  int* segImg(NULL);

  Scalar* img1    = (Scalar*)  mxGetPr(prhs[0]);
  Scalar* img2    = (Scalar*)  mxGetPr(prhs[1]);
  Scalar* img3    = (Scalar*)  mxGetPr(prhs[2]);
  Scalar* img4    = (Scalar*)  mxGetPr(prhs[3]);

  Scalar* p2d_    = (Scalar*)  mxGetPr(prhs[4]); // K_l^-1
  Scalar* K_      = (Scalar*)  mxGetPr(prhs[5]); // K_r
  const mxArray* Segments  = (prhs[6]);

  Scalar* normals = (Scalar*) mxGetPr(prhs[ 7]);
  thresh          = (Scalar) *mxGetPr(prhs[ 8]);
  lambda          = (Scalar) *mxGetPr(prhs[ 9]);
  theta           = (Scalar) *mxGetPr(prhs[10]);

  const mxArray* edges    = (prhs[11]);
  const mxArray* weights  = (prhs[12]);

  // optional: always 0,0,0 if not needed
  Scalar* mc_       =  (Scalar*) mxGetPr(prhs[13]);// right cam

  int* oobSegsR  = NULL;
  int* oobSegsLT = NULL;
  int* oobSegsRT = NULL;

  // optional:
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[14]);
  Scalar* centers   =  (Scalar*) mxGetPr (prhs[15]);
  Scalar* MC_       =  (Scalar*) mxGetPr (prhs[16]);//right cam

  if (nrhs > 17)
    rotJump       =  (Scalar) *mxGetPr(prhs[17]);

  if (nrhs > 18)
    depthJump     =  (Scalar) *mxGetPr(prhs[18]);

  if (nrhs > 19)
    dispMax     =  (Scalar) *mxGetPr(prhs[19]);

  //////////////////
  if (nrhs > 20)
    occThresh   =  (Scalar) *mxGetPr(prhs[20]);

  if (nrhs > 21)
    segImg      =  (int*) mxGetPr(prhs[21]); // for census - not flipped

  if (nrhs > 22)
    oobSegsR   =  (int*) mxGetPr(prhs[22]);

  if (nrhs > 23)
    oobSegsLT   =  (int*) mxGetPr(prhs[23]);

  if (nrhs > 24)
    oobSegsRT  =  (int*) mxGetPr(prhs[24]);

  if (nrhs > 25)
    motMax     =  (Scalar) *mxGetPr(prhs[25]);

  if (nrhs > 26)
    halfEdgeX       =  (Scalar*) mxGetPr(prhs[26]);

  if (nrhs > 27)
    halfEdgeY       =  (Scalar*) mxGetPr(prhs[27]);

  if (nrhs > 28)
    halfEdgeXY      =  (Scalar*) mxGetPr(prhs[28]);

  if (nrhs > 29)
    halfEdgeiXY     =  (Scalar*) mxGetPr(prhs[29]);

  if (nrhs > 30)
    origSegImg      =  (int*) mxGetPr(prhs[30]);

  if (nrhs > 31)
    Kr_             = (Scalar*)  mxGetPr(prhs[31]); // K_r
  /// if Kr == Kl
  if (Kr_==NULL) Kr_ = K_;

  //////////////////////////
  Scalar outOfBoundsThresh = occThresh;

  printf("rotJump %f \n", rotJump);
  printf("depthJump %f \n", depthJump);
  printf("dispMax %f \n", dispMax);
  printf("motMax %f \n", motMax);
  printf("occThresh %f \n", occThresh);
  printf("outOfBoundsThresh %f \n", outOfBoundsThresh);
  printf("data thresh %f \n", thresh);
  printf("motion smooth theta %f \n", theta);
  printf("general smoothness lambda %f \n", lambda);
  printf("size int %d \n", sizeof(int) );

  Scalar* output, *output2 = NULL;
  Scalar* output3, *output4 = NULL;

  M           = (size_t) mxGetM(prhs[0]);
  N           = (size_t) mxGetN(prhs[0]);

  int nNormals = max( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );
  int nStep    = min( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );

  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  printf ("nSegments %d, nNormals %d \n", nSegments, nNormals);
  //  mexEvalString("drawnow");

  std::vector<bool> oobOutsideLR(M*N, false);
  if (oobSegsR != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsR[i] ) oobOutsideLR[i]=true;//disp

  std::vector<bool> oobOutsideLLT(M*N, false);
  if (oobSegsLT != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsLT[i] ) oobOutsideLLT[i]=true;//flow

  // maybe use later
  std::vector<bool> oobOutsideLRT(M*N, false);
  if (oobSegsRT != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsRT[i] ) oobOutsideLRT[i]=true;//left right t+1

  /// these could be merged with left, right t+1 
  std::vector<bool> oobOutsideRRT(M*N, false); //right right t+1
  if (oobSegsRT !=NULL)
  {
    for (int i= 0; i< M*N; i++)
      if ( oobSegsR[i] || oobSegsRT[i] ) oobOutsideRRT[i]=true;
  }
  else 
    if (oobSegsR != NULL)
    {
      for (int i= 0; i< M*N; i++)
        if ( oobSegsR[i] ) oobOutsideRRT[i]=true;
    }

  /// rrt
  std::vector<bool> oobOutsideLTRT(M*N, false);  //left t+1 - right t+1
  if (oobSegsRT != NULL  )
  {
    for (int i= 0; i< M*N; i++)
      if ( oobSegsLT[i] || oobSegsRT[i] ) oobOutsideLTRT[i]=true;
  }
  else // no third view
  {
    if (oobSegsRT != NULL)
      for (int i= 0; i< M*N; i++)
        if ( oobSegsLT[i] ) oobOutsideLTRT[i]=true;
  }

  mwSize dims[2];dims[0] = N;dims[1] = M;
  // Create the output array
  {
    // scores
    plhs[0]      =  mxCreateDoubleMatrix(nSegments, 1, mxREAL);
    output  = (Scalar*) mxGetPr(plhs[0]);
  }

  ///////////////////////////////////////////////////////

  EvalEnergyStore<Scalar> enStore;
  enStore.prepareNormals( normals, nNormals, nStep );
  enStore.setRotTra( rotations, nNormals );

  genHomography<Scalar> gHom1 (K_, MC_, mc_, Kr_);//t,t, left right
  genHomography<Scalar> gHom2 (K_); // t,t+1 left
  genHomography<Scalar> gHom3 (K_, MC_, mc_, Kr_);// t,t+1 left right

  gHom1.setNormals( enStore.getNormals() );

  gHom2.setNormals( enStore.getNormals() );
  gHom2.setRotTra(  enStore.getTra(), enStore.getRot() );
  
  gHom3.setNormals( enStore.getNormals() );
  gHom3.setRotTra(  enStore.getTra(), enStore.getRot() );

  gHom2.setCenters( centers );
  gHom3.setCenters( centers );

  genWarp<double> gWarp1(N,M, img2, p2d_);// do warps
  genWarp<double> gWarp2(N,M, img3, p2d_);// do warps
  genWarp<double> gWarp3(N,M, img4, p2d_);// do warps
  genWarp<double> gWarp0(N,M, img1);// no warps -> orig image

  gWarp1.setHom( &gHom1 );//rt
  gWarp2.setHom( &gHom2 );//lt1
  gWarp3.setHom( &gHom3 );//rt1

  accumulateWarp<Scalar> aW01(N,M,thresh,outOfBoundsThresh); // disparity
  accumulateWarp<Scalar> aW02(N,M,thresh,outOfBoundsThresh); // flow
  accumulateWarp<Scalar> aW13(N,M,thresh,outOfBoundsThresh); // flow right
  accumulateWarp<Scalar> aW23(N,M,thresh,outOfBoundsThresh); // disparity t+1

#ifdef _use_5th_view
  accumulateWarp<Scalar> aW03(N,M,thresh,outOfBoundsThresh); // cross view lt->rt1
  aW03.setWarp1(&gWarp3);  aW03.setWarp2(&gWarp0);// to enable testing of maxDisp
  aW03.setRefImage(img1);
#endif


  aW01.setMaxDisp(dispMax); // if set a test is done before adding up values
  aW01.setWarp1(&gWarp1);  aW01.setWarp2(&gWarp0);// to enable testing of maxDisp
  aW02.setWarp1(&gWarp2);  aW02.setWarp2(&gWarp0);// takes indices from the first warp
  aW02.setMaxMot(motMax);
  aW13.setWarp1(&gWarp1);  aW13.setWarp2(&gWarp3);
  aW23.setWarp1(&gWarp2);  aW23.setWarp2(&gWarp3);

  accumulate2FrameData<double> aFD(N,M);
  // 4(3) images to be warped
  aFD.addWarp(&gWarp1); 
  // 4 equations from images:
  aFD.addAccumulateWarp(&aW01);  

  aW01.setRefImage(img1);
  aW02.setRefImage(img1);
  aW13.setRefImage(img1);
  aW23.setRefImage(img1);

  aFD.addWarp(&gWarp2);  aFD.addWarp(&gWarp3);  aFD.addWarp(&gWarp0);
  aFD.addAccumulateWarp(&aW02);  aFD.addAccumulateWarp(&aW13);  aFD.addAccumulateWarp(&aW23);

#ifdef _use_5th_view
  aFD.addAccumulateWarp(&aW03);
#endif

  typedef genWarp<Scalar>::P3 P3;

  int runs = nNormals/nSegments;

  std::vector< std::vector<int> > secondProposals;
  secondProposals.resize(runs+1);
  std::vector< std::vector<Scalar> > dataScores;
  dataScores.resize(runs+1);

  std::vector< std::vector< std::vector<int> > >    freeTrial       ;
  std::vector< std::vector< std::vector<Scalar> > > freeScoresTrial ;
  freeTrial.resize(runs+1);
  freeScoresTrial.resize(runs+1);

  std::vector<int> currentSolution( nSegments, 0 );
  std::vector<int> oldSolution( nSegments, 0 );

  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      currentSolution[i]=i;

  // first data term, before normal flipping: all 0
  int kkk = 0;
  while(nNormals >= (kkk+1)*nSegments )
  {
    secondProposals[kkk].resize(nSegments,0);
    for ( int i = 0;i < nSegments;i++ )
      (secondProposals[kkk])[i]=i+kkk*nSegments;

#ifdef __QUIET__
	  std::clock_t start(std::clock());
#endif

    aFD.computeVecScores (Segments, nNormals, segImg, secondProposals[kkk], oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);

    dataScores[kkk]      = aFD.getAllScores();
    freeTrial[kkk]       =  aFD.getNFreeVariables();
    freeScoresTrial[kkk] =  aFD.getFreePartialScores();

    Scalar fullDataScore(0);

    for(int i = 0; i < dataScores[kkk].size(); i++)
      fullDataScore += (dataScores[kkk])[i];
#ifdef __QUIET__
    printf("Time Data scores are %d: %f \n", kkk, fullDataScore);
    std::clock_t stop(std::clock());
    printf("Data took %f\n", double(stop-start)/CLOCKS_PER_SEC );
#endif
    kkk++;
  }
  std::vector<Scalar> currentData = dataScores[0];
  std::vector<Scalar> saveStartData = dataScores[0];

  aFD.computeVecScores (Segments, nNormals, segImg, currentSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);
  std::vector<Scalar> dataScoresRepeat = aFD.getAllScores();

  std::vector< std::vector<int> >&    freeCurrent       =  aFD.getNFreeVariables();
  std::vector< std::vector<Scalar> >& freeScoresCurrent =  aFD.getFreePartialScores();
  //--
  OcclusionMapBufferSeg<Scalar> ombs1 (N,M, &gHom1, origSegImg, centers, K_ );// K_: left camera matrix
  OcclusionMapBufferSeg<Scalar> ombs2 (N,M, &gHom2, origSegImg, centers, K_ );
  OcclusionMapBufferSeg<Scalar> ombs3 (N,M, &gHom3, origSegImg, centers, K_ );
  HigherOrderConverter<Scalar> hoc(nSegments);
  //--

#ifdef __QUIET__
  for (int i=0;i<nSegments;i++)
  {
    if ( dataScoresRepeat [i] != currentData[i] )
    {
      printf("ALERT Start diff score: %d:%d - %f %f, start: %f\n", i, currentSolution[i], dataScoresRepeat [i], currentData[i], saveStartData[i]);
    }
  }
#endif

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  // now smoothness cost setup:

  EvalEnergyFullFrame<Scalar> evalSmooth;
  evalSmooth.setRotJump  (rotJump);
  evalSmooth.setDepthJump(depthJump);

  evalSmooth.setGamma ( 1.0 );
  evalSmooth.setRotWeight( theta );
  evalSmooth.set2dMotionMatrix ( K_, Kr_, MC_, mc_, 1.0); // jump at 0.5 pixel -> needs to set jump smaller

  enStore.flipNormals();
  evalSmooth.setNormals( enStore.getNormals() );
  evalSmooth.setRotTra( enStore.getRot(), enStore.getTra() );

  if ( segImg != NULL)
  {
    evalSmooth.setHalfEdges(halfEdgeX, halfEdgeY);
    if ( halfEdgeiXY != NULL)
      evalSmooth.setCrossHalfEdges(halfEdgeXY, halfEdgeiXY);

    evalSmooth.prepareOwnWeights( nSegments, edges, centers, origSegImg, N, M );
  }
  else
    evalSmooth.prepare( nSegments, edges, weights, centers );

  int numEdges = evalSmooth.edges_num( );

  Scalar OLD_ENERGY = std::numeric_limits<Scalar>::max();
  Scalar ENERGY(0), occEnergy(0);

  /// fixed anyway
  const std::vector<int>& idk=evalSmooth.getIdk();
  const std::vector<int>& idl=evalSmooth.getIdl();
  const std::vector<Scalar>& f11=evalSmooth.getF11();
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  const std::vector<Scalar>& f10=evalSmooth.getF10();
  const std::vector<Scalar>& f01=evalSmooth.getF01();     
  std::vector<int> noSolutions;noSolutions.reserve(nSegments);
  
  int maxValue = 1 << 29;
  int non_sub = 0;
  int non_sol = 0;

  QPBO<int> q(4*nSegments, 4*numEdges); // max number of nodes & edges

  ENERGY = getEnergy ( currentSolution, nSegments, dataScores[0], evalSmooth, theta, lambda);
#ifdef __QUIET__  
  printf("------------ Start Energy is %f --------------\n", ENERGY);
  printf("Start Energies (full=data+smooth) %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, dataScores[0], evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, dataScores[0] ), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda) );
#endif

  std::clock_t fullStart(std::clock());
  int fullIts(0);
  // fuse all proposals in a row
  while( fullIts < 3)
  {
    fullIts++;
    kkk=1;
  while(nNormals >= (kkk+1)*nSegments )
  {
      {
        Scalar occPenalty = occThresh; // best thing is to do nothing !
        /////////////////////////////////////

        // now add the occlusion edges 
        ombs1.buildList(currentSolution, secondProposals[kkk] );
        ombs1.generateOcclusionLists();
        ombs2.buildList(currentSolution, secondProposals[kkk] );
        ombs2.generateOcclusionLists();
        ombs3.buildList(currentSolution, secondProposals[kkk] );
        ombs3.generateOcclusionLists();
        std::vector< occSegList<Scalar> >& _occList  = ombs1.getOccList();
        std::vector< occSegList<Scalar> >& _occList2 = ombs2.getOccList();
        std::vector< occSegList<Scalar> >& _occList3 = ombs3.getOccList();

        BinaryConverter<Scalar> bc_0( occPenalty, _occList , freeCurrent[0], freeTrial[kkk][0] );
        BinaryConverter<Scalar> bc_1( occPenalty, _occList2, freeCurrent[1], freeTrial[kkk][1] );
        BinaryConverter<Scalar> bc_2( occPenalty, _occList, _occList3,  freeCurrent[2], freeTrial[kkk][2] );
        BinaryConverter<Scalar> bc_3( occPenalty, _occList2, _occList3, freeCurrent[3], freeTrial[kkk][3] );

        bc_0.setPenaltiesPerSegment( freeScoresCurrent[0], freeScoresTrial[kkk][0] );
        bc_1.setPenaltiesPerSegment( freeScoresCurrent[1], freeScoresTrial[kkk][1] );
        bc_2.setPenaltiesPerSegment( freeScoresCurrent[2], freeScoresTrial[kkk][2] );
        bc_3.setPenaltiesPerSegment( freeScoresCurrent[3], freeScoresTrial[kkk][3] );

        ////////////////////////
        bc_0.generateHigherOrderTerms();
        bc_1.generateHigherOrderTerms();
        bc_2.generateHigherOrderTerms();
        bc_3.generateHigherOrderTerms();

        hoc.clear();
        hoc.addNonSubs ( bc_0.get_nSubmods() );
        hoc.addSubs ( bc_0.get_Submods(), bc_0.getUnaries() );
        hoc.addNonSubs ( bc_1.get_nSubmods() );
        hoc.addSubs ( bc_1.get_Submods(), bc_1.getUnaries() );
        hoc.addNonSubs ( bc_2.get_nSubmods() );
        hoc.addSubs ( bc_2.get_Submods(), bc_2.getUnaries() );
        hoc.addNonSubs ( bc_3.get_nSubmods() );
        hoc.addSubs ( bc_3.get_Submods(), bc_3.getUnaries() );

        hoc.convert();
        q.AddNode( hoc.get_nVars() );
      }
      ////////////
  Scalar worstDataScore(0);
  Scalar bestDataScore(0);

  for (int i=0; i< nSegments;i++)
  {
    worstDataScore += max( currentData[i], dataScores[kkk][i] );
    bestDataScore  += min( currentData[i], dataScores[kkk][i] );
  }

  Scalar worstSmoothScore = evalSmooth.compute_score_Fuse( currentSolution, secondProposals[kkk] );

  Scalar scale = (Scalar)(maxValue) / (lambda * worstSmoothScore + worstDataScore);

  Scalar scaleSmo = scale * lambda;
  // add unary data terms:
  for ( int j = 0; j< nSegments; j++ )
  {
    q.AddUnaryTerm(j, scale*(currentData[j]     ), 
                      scale*(dataScores[kkk][j] ) );
  }      

  for( int j=0; j < idk.size(); j++ )
  {
    if (f00[j]+f11[j] > f10[j]+f01[j])
      non_sub++;

    q.AddPairwiseTerm(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);
  }

   // change QPBO for occlusion handling
   int nVars = hoc.get_nVars();
   if (nVars > nSegments || hoc.get_nEdges() > 0 || hoc.get_b0().size() + hoc.get_b1().size() > 0)
      {
        // probing takes ages - do improve instead !
        if ( hoc.get_nVars() > 6000 || hoc.get_nEdges()> 4*numEdges) 
        {
          kkk++;
          hoc.clear();
          q.Reset();
          continue;
        }
        std::vector< std::pair<int, Scalar> >& occ_b0 = hoc.get_b0();
        std::vector< std::pair<int, Scalar> >& occ_b1 = hoc.get_b1();

        std::vector< edgePenalty<Scalar> >& occ_b00 = hoc.get_b00();
        std::vector< edgePenalty<Scalar> >& occ_b01 = hoc.get_b01();
        std::vector< edgePenalty<Scalar> >& occ_b10 = hoc.get_b10();
        std::vector< edgePenalty<Scalar> >& occ_b11 = hoc.get_b11();

        for (int j=0; j< occ_b00.size();j++)
          q.AddPairwiseTerm(occ_b00[j].seg_i, occ_b00[j].seg_j, scale*occ_b00[j].f11,  0,0,0);
        for (int j=0; j< occ_b01.size();j++)
          q.AddPairwiseTerm(occ_b01[j].seg_i, occ_b01[j].seg_j, 0, scale*occ_b01[j].f11, 0,0);
        for (int j=0; j< occ_b10.size();j++)
          q.AddPairwiseTerm(occ_b10[j].seg_i, occ_b10[j].seg_j, 0,0, scale*occ_b10[j].f11, 0);
        for (int j=0; j< occ_b11.size();j++)
          q.AddPairwiseTerm(occ_b11[j].seg_i, occ_b11[j].seg_j, 0,0,0, scale*occ_b11[j].f11 );

        for (int j=0; j< occ_b0.size();j++)
          q.AddUnaryTerm(occ_b0[j].first, scale*(occ_b0[j].second), 0);
        for (int j=0; j< occ_b1.size();j++)
          q.AddUnaryTerm(occ_b1[j].first, 0, scale*(occ_b1[j].second));

        hoc.clear();
        // otherwise all edges should be unique/called only once!
        q.MergeParallelEdges();
      }
      /////////////////////////////

  q.Solve();
  q.ComputeWeakPersistencies();

#define useProbe
#ifndef useProbe
  for (int j = 0; j < nSegments;j++)
  {
    int x = q.GetLabel(j);
    if (x<0)
      non_sol++;

    if (x==1) {currentSolution[j] = secondProposal[j];};
  }

#else
  // try to improve if something is not solved completely:
  int doImprove = 0;
#ifdef _DEBUG
  std::vector<int> noSolution ( nSegments, 0 );
#endif
  noSolutions.clear();
#ifdef _DEBUG
  std::vector<int> new1old0(nSegments, 0);
  std::list<int>   newListed;
#endif
  int allNew (0), n_unsolved(0);
  for (int j = 0; j < nSegments;j++)
  {
    oldSolution[j] = currentSolution[j];
    int x = q.GetLabel(j);
    if (x<0)
    {
      doImprove = 1;
      non_sol++;
#ifdef _DEBUG
      noSolution[j]++; 
#endif
      noSolutions.push_back(j);
      n_unsolved++;
    }
    if (x==1) {
      currentSolution[j] = secondProposals[kkk][j];
      currentData[j]     = dataScores[kkk][j];
      allNew++;
#ifdef _DEBUG
      new1old0[j]=1;
      newListed.push_back(j);
#endif
    };
  }

  for (int j = nSegments; j < q.GetNodeNum();j++)
     if ( q.GetLabel(j) <0 )
       n_unsolved++;

  int doProbe = 1;int impIts =0;int allImpIts = 0;
  if (doProbe && doImprove && n_unsolved < limitUnsolved)
  {
#ifdef __QUIET__
    printf(" Start probing :%d vars, %d undef ", nVars, noSolutions.size() );mexEvalString("drawnow");
#endif
    for (int j = 0; j < nVars;j++)
    {
      int x = q.GetLabel(j);
      if (x>=0)
        q.SetLabel(j,x);
    }

    int* mapping(NULL);
    if (maxProbeRuns>0)
    {
        QPBO<int>::ProbeOptions options;
        options.C = 200000;
        options.dilation = 1;
        options.weak_persistencies = 1;
        mapping = (int*) malloc( sizeof(int) * nVars );

      q.Probe( mapping, options );
      q.ComputeWeakPersistencies();
#ifdef __QUIET__
      printf(" probe 0 " );      mexEvalString("drawnow");
#endif
    int* mapping2 = (int*) malloc( sizeof(int) * nVars );
    for (int i=1;i < maxProbeRuns; i++)
    {
      q.Probe( mapping2, options );
      q.ComputeWeakPersistencies();
#ifdef __QUIET__
      printf("/%d ", i );      mexEvalString("drawnow");
#endif
      q.MergeMappings( nVars, mapping, mapping2 );
    }
    free(mapping2);
    }

    const int maxInner = 10;//10;
    const int maxReps = 20; int imp = 0;
    while (imp < maxInner && impIts < maxReps)
    {
      impIts++;
      for(imp=0; imp<maxInner; imp++ )
        if ( q.Improve() )  {break;printf("Improved after probing\n");};

      allImpIts += imp;
      if (imp == maxInner) break;// out of while
    }

    for (int j=0;j < noSolutions.size(); j++)
    {
      int pos (noSolutions[j]), x(-1);

      if (maxProbeRuns>0)
         x = q.GetLabel( mapping[pos]/2 );
      else
         x = q.GetLabel(pos);

      if( (x>=0) && ( maxProbeRuns<=0 || (x+ mapping[pos]) % 2) )
      {
        currentSolution[pos] = secondProposals[kkk][pos];
        currentData[pos] = dataScores[kkk][pos];
        allImpIts++;
        impIts++;
        allNew++;
      }
    }
    if (maxProbeRuns > 0 )
      free( mapping );
  }
#endif

  if (allNew > 0)
  {
    enStore.flipNormals();
    aFD.computeVecScores (Segments, nNormals, segImg, currentSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);
    enStore.flipNormals();
    currentData = aFD.getAllScores();
    freeCurrent       =  aFD.getNFreeVariables();
    freeScoresCurrent =  aFD.getFreePartialScores();
  }

  OLD_ENERGY = ENERGY;
  occEnergy = getOcclusionEnergy ( currentSolution, occThresh, ombs1, ombs2, ombs3, freeScoresCurrent, freeCurrent );

  ENERGY = getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda) + occEnergy;

#ifdef __QUIET__
  printf("Energies (full=data+smooth+occ) comp %f = (%f,%f,%f)\n", ENERGY, getDataEnergy ( currentSolution, nSegments, currentData), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda), occEnergy );
#endif

  q.Reset();
  kkk++;
  }
}

#ifdef __QUIET__
  std::vector<int> realSolution( nSegments, 0 );
  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      realSolution[i]=i;
  // not working since flipped
  enStore.flipNormals();
  aFD.computeVecScores (Segments, nNormals, segImg, currentSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
  std::vector<Scalar> dataScores1 = aFD.getAllScores();

  std::vector<int> realSolution2( nSegments, 0 );
  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      realSolution2[i]=i+nSegments;

  aFD.computeVecScores (Segments, nNormals, segImg, realSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
  std::vector<Scalar> dataScores1Test0 = aFD.getAllScores();
  aFD.computeVecScores (Segments, nNormals, segImg, realSolution2, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
  std::vector<Scalar> dataScores1Test1 = aFD.getAllScores();

  for (int i=0;i<nSegments;i++)
  {
    if ( dataScores1 [i] != currentData[i] )
    {
      printf("diff score: seg:%d : sol:%d - eval:%f store%f, start: %f, stored by id %f\n", i, currentSolution[i], dataScores1 [i], currentData[i], saveStartData[i], dataScores[currentSolution[i]/nSegments][i]);
      printf("new score: %d:%d - %f %f\n", i, i+nSegments, dataScores1Test0[i], dataScores1Test1[i] );
    }
  }

  occEnergy = getOcclusionEnergy ( currentSolution, occThresh, ombs1, ombs2, ombs3, freeScoresCurrent, freeCurrent );

  Scalar d0 = getDataEnergy ( currentSolution, nSegments, currentData);
  Scalar s0 = getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda);
  Scalar d1 = getDataEnergy (  currentSolution, nSegments, dataScores1);
  Scalar s1 = getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda);
  printf("Energies (full=data+smooth+occ) comp %f = (%f,%f,%f) 0:N %f = (%f,%f,%f)\n", ENERGY, d0+s0+occEnergy, d0,s0,occEnergy, d1,s1,occEnergy );

  printf("QPBO took %f\n", double(std::clock()-fullStart)/CLOCKS_PER_SEC);
  printf("ENERGY: %f  NON submodularity detected %d, no solutions: %d\n", ENERGY, non_sub, non_sol);
#endif

  for (int i=0;i < currentSolution.size(); i++)
    output[i] =currentSolution[i];

    ombs1.buildSingleList( currentSolution );
    ombs1.generateOcclusionLists();
    ombs2.buildSingleList( currentSolution );
    ombs2.generateOcclusionLists();
    ombs3.buildSingleList( currentSolution );
    ombs3.generateOcclusionLists();

    std::vector<OcclusionMapBufferSeg<Scalar>::SegLbl> testFO1 = ombs1.getFullOcc();
    std::vector<OcclusionMapBufferSeg<Scalar>::SegLbl> testFO2 = ombs2.getFullOcc();
    std::vector<OcclusionMapBufferSeg<Scalar>::SegLbl> testFO3 = ombs3.getFullOcc();

    if (nlhs > 1)
    {
      plhs[1]           =  mxCreateDoubleMatrix(testFO1.size(), 1, mxREAL);
      Scalar *outputDI  = (Scalar*) mxGetPr(plhs[1]);
      for (int i=0;i < testFO1.size(); i++)
        outputDI[i] = testFO1[i].first;
    }
    if (nlhs > 2)
    {
      plhs[2]           =  mxCreateDoubleMatrix(testFO2.size(), 1, mxREAL);
      Scalar *outputDI  = (Scalar*) mxGetPr(plhs[2]);
      for (int i=0;i < testFO2.size(); i++)
        outputDI[i] = testFO2[i].first;
    }
    if (nlhs > 3)
    {
      plhs[3]           =  mxCreateDoubleMatrix(testFO3.size(), 1, mxREAL);
      Scalar *outputDI  = (Scalar*) mxGetPr(plhs[3]);
      for (int i=0;i < testFO3.size(); i++)
        outputDI[i] = testFO3[i].first;
    }

}
#undef useProbe
#undef __Small_Improvement__
#endif // __QPBO_OCC_2FRAME__
