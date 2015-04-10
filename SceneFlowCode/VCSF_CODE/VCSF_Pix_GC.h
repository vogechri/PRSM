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

#ifndef __VCSF_PIX_GC__
#define __VCSF_PIX_GC__
///////////////////////////////////////////////////////////////
/// Compute mapping S from pixels to segments given fixed P ///
///////////////////////////////////////////////////////////////
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "DataDefinitionsVC.h"

#ifdef _qpboVersion_
#include "QPBOCutting.h" // use QPBO to solve
#else
#include "GraphCutting.h"  // use LSA Aux to solve
#endif

#ifdef WIN32
#include <windows.h>
#endif
// TODO: centers should again be distinguished from the expansion centers!
#ifdef max 
#undef max
#endif
#ifdef min
#undef min
#endif

#include <numeric>

#include "EvalEnergyFull2FramePixel.h"
#include "EvalEnergyStore.h"
#ifdef _MotionChecking_
  #include "SensibleMotionCheck.h"
#endif
#include "DataSegments.h"

#include <vector>
#include <algorithm>

#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

#include <ctime>
#include <cstdlib>

using namespace std;
using namespace Math;

// fast for debugging:
//#define _shortCut_forTesting_
//#define __Energy_Check__Data__ 
//////////////////////////////////////////////////////////////////////////////////////

template<typename Scalar>
void addUnariesBinaries( GraphCutContainer<Scalar>& q, const typename std::vector< Unary<Scalar> >& unaries, 
  const typename std::vector< Binary<Scalar> >& binaries, 
  const Scalar scale, int& non_sub)//, std::vector<Scalar>& min0min1max0max1)
{
  //  int nVars = min0min1max0max1.size()/4;
  for (int k=0; k< unaries.size();k++)
  {
    q.AddUnaryTerm(unaries[k].var, scale*( unaries[k].node[0] ) ,scale*( unaries[k].node[1] ) );
  }
  for( int j=0; j < binaries.size(); j++ )
  {
    typename Datasegments<Scalar>::P4 ff (scale * binaries[j].edge);
    if (ff[0]+ff[3] > ff[1]+ff[2])
      non_sub++;

    if (ff[1] != ff[0] || ff[0] != ff[2]  || ff[0] != ff[3] )
    {
      q.AddPairwiseTermLSAAux(binaries[j].segI, binaries[j].segJ, ff);
    }
  }
};


template<typename Scalar>
void addSmoothnessConstraints ( GraphCutContainer<Scalar>& q, const std::vector<Scalar>& f0, const std::vector<Scalar>& f1, const std::vector<Scalar>& f00, 
  const std::vector<Scalar>& f01, const std::vector<Scalar>& f10, const std::vector<Scalar>& f11, 
  const std::vector<int>& idk, const std::vector<int>& idl, int start, int end, Scalar& scaleSmo, int& non_sub, const int nVars)//,
{
  for (int j=start; j<end; j++ )
  {
    q.AddUnaryTerm(j, scaleSmo*f0[j], scaleSmo*f1[j]);
  }

  for( int j=0; j < idk.size(); j++ )
  {

    if (f00[j]+f11[j] > f10[j]+f01[j])
      non_sub++;

    typename Datasegments<Scalar>::P4 edge(scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);

    if (f11[j] != f00[j] || f00[j] != f01[j] || f00[j] != f10[j] )
      q.AddPairwiseTermLSAAux(idk[j], idl[j], edge );
    //      q.AddPairwiseTermTrunc(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);
  }
}

// std::min(cam1,cam2) map to be unique - so ?? 
int getQuadId(int cam1, int time1, int cam2, int time2, int nCams, int nTimeSteps) 
{
  int c1,c2,t1,t2;
  if ( cam1 < cam2 || (cam1 == cam2 && time1< time2) )
  {c1 = cam1;c2=cam2;t1=time1;t2=time2;}
  else
  {c1 = cam2;c2=cam1;t1=time2;t2=time1;}

  return c1*nCams*nTimeSteps*nTimeSteps + c2*nTimeSteps*nTimeSteps + t1*nTimeSteps + t2;
}

int getQuadId(Datasegments<double>::P4i ctct, int nCams, int nTimeSteps) 
{
  return getQuadId(ctct[0], ctct[1], ctct[2], ctct[3], nCams, nTimeSteps);
}

Datasegments<double>::P4i getQuadIdInv(int id, int nCams, int nTimeSteps) 
{
  int c1,c2,t1,t2;

  t2  = id % nTimeSteps;
  id -= t2;id /= nTimeSteps;
  t1  = id % nTimeSteps;
  id -= t1;id /= nTimeSteps;
  c2  = id % nCams;
  id -= c2;id /= nCams;
  c1  = id % nCams;
  id -= c1;
  assert(id==0);
  return Datasegments<double>::P4i(c1,t1,c2,t2);
}

template<typename Scalar, typename camTimePair>
Scalar getSmoothnessEnergy( std::vector< EvalEnergyFullFramePixel<Scalar> >& evalSmooth, 
  const std::vector<camTimePair>& smoothCams, 
  const std::vector< int >& dsgId, const int nCams, const int nTimeSteps, 
  const std::vector < Datasegments<Scalar>* >& dsgX, const Scalar lambda )
{
  Scalar allSmoothScore =0;
  for(int sc=0;sc<evalSmooth.size();sc++)
  {

    int pid =0;
    camTimePair ct = smoothCams[sc];
    // a) find the boxes:
    //
    int qid_10 = getQuadId( ct[0], ct[1], ct[0], ct[3], nCams, nTimeSteps);
    int qid_01 = getQuadId( ct[0], ct[1], ct[2], ct[1], nCams, nTimeSteps);
    int qid_11 = getQuadId( ct[0], ct[1], ct[2], ct[3], nCams, nTimeSteps);

    int qid10 = dsgId[ qid_10 ];// saem cam other time
    int qid01 = dsgId[ qid_01 ];// same time other cam
    int qid11 = dsgId[ qid_11 ];
    // b) get whether fw or bw direction of dsg
    //

    typename Datasegments<Scalar>::P4i dummy1 = getQuadIdInv( qid_01, nCams, nTimeSteps );
    typename Datasegments<Scalar>::P4i dummy2 = getQuadIdInv( qid_10, nCams, nTimeSteps );
    typename Datasegments<Scalar>::P4i dummy3 = getQuadIdInv( qid_11, nCams, nTimeSteps );
    bool i01 = getQuadIdInv( qid_01, nCams, nTimeSteps ) == typename Datasegments<Scalar>::P4i( ct[0], ct[1], ct[2], ct[1]);
    bool i10 = getQuadIdInv( qid_10, nCams, nTimeSteps ) == typename Datasegments<Scalar>::P4i( ct[0], ct[1], ct[0], ct[3]);
    bool i11 = getQuadIdInv( qid_11, nCams, nTimeSteps ) == typename Datasegments<Scalar>::P4i( ct[0], ct[1], ct[2], ct[3]);

    dataElem<Scalar> de10 = dsgX[qid10]->getStoredElem(pid);
    dataElem<Scalar> de01 = dsgX[qid01]->getStoredElem(pid);
    dataElem<Scalar> de11 = dsgX[qid11]->getStoredElem(pid);

    const std::vector<int>& loc2glob = i01 ? de01.seg2varFW : de01.seg2varBW;

    Scalar SmoothScore = evalSmooth[sc].compute_score_combiDepth_consistentScore(  
      i01 ? dsgX[qid01]->getvHoms() : dsgX[qid01]->getviHoms(), 
      i10 ? dsgX[qid10]->getvHoms() : dsgX[qid10]->getviHoms(), 
      i11 ? dsgX[qid11]->getvHoms() : dsgX[qid11]->getviHoms() );
    allSmoothScore += SmoothScore;
#ifdef _writeEnergies_
    printf("cam/time (%d,%d) smooth Score: %.2f\n", ct[0], ct[1], SmoothScore*lambda);
#endif
  }
  return allSmoothScore;
}

//////////////////////////////////////////////////////////////////////////////
// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void run_VCSF( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

  size_t N,M;
  size_t elements;
  char buffer[256];

  if (nrhs < 6)
    mexErrMsgTxt("Must have at least 6 input arguments");

  Scalar thresh = 0.1;
  Scalar lambda = 0.02;
  Scalar theta  = 100.0;
  Scalar dispMax = 32.;
  Scalar pixJump = 1.0;
  Scalar rotJump(15.);          //   setRotJump (rotJump)
  Scalar depthJump(sqrt(1300.));//   setRotJump (rotJump)
  Scalar occThresh(0.8);
  Scalar outOfBoundsThresh(0.8);
  Scalar motMax = 220; // the maximal motion in the scene
  Scalar gamma = 1.0;
  Scalar gridSize = 16;
  Scalar patchSize = 18;
  Scalar ot = 0.015; // min set anyway 0.015 good enough - maybe reduce min as well
  Scalar tempPotts = 0.25;//0.1//0.25;// less here in the first place, later more
  Scalar segWeight = 0.1;
  Scalar doAuto =1.; // update autoScores based on data term

  Scalar* halfEdgeX(NULL);
  Scalar* halfEdgeY(NULL);
  Scalar* halfEdgeXY(NULL);
  Scalar* halfEdgeiXY(NULL);
  std::vector<Scalar*> pict;
  std::vector<int*> segImg_t; //
  int *proposals(NULL);       // defines the proposal to expand
  int nTimeSteps = 2;

  pict.push_back ( (Scalar*)  mxGetPr(prhs[0]) );
  pict.push_back ( (Scalar*)  mxGetPr(prhs[1]) );

  mwSize nDims = mxGetNumberOfDimensions(prhs[0]);
  assert(nDims == 3);
  const mwSize* imgDims = mxGetDimensions(prhs[0]);
  const int nCams   = imgDims[2];
  assert( mxGetNumberOfElements(prhs[0]) == mxGetNumberOfElements(prhs[1]) );
  M = imgDims[0];
  N = imgDims[1];

  Scalar* parameters = (Scalar*) mxGetPr(prhs[ 4]);// all parameters:
  int nParameters    = (size_t)  mxGetNumberOfElements(prhs[4]);

  if ( nParameters > 0) thresh            = parameters[0];
  if ( nParameters > 1) lambda            = parameters[1];
  if ( nParameters > 2) theta             = parameters[2];
  if ( nParameters > 3) rotJump           = parameters[3];
  if ( nParameters > 4) depthJump         = parameters[4];
  if ( nParameters > 5) dispMax           = parameters[5];
  if ( nParameters > 6) motMax            = parameters[6];
  if ( nParameters > 7) outOfBoundsThresh = parameters[7];
  if ( nParameters > 8) occThresh         = parameters[8];
  if ( nParameters > 9) gridSize          = parameters[9];
  if ( nParameters > 10) patchSize        = parameters[10];
  if ( nParameters > 11) segWeight        = parameters[11];
  if ( nParameters > 12) tempPotts        = parameters[12];
  if ( nParameters > 13) ot               = parameters[13];
  if ( nParameters > 14) doAuto           = parameters[14];

  segImg_t.push_back( (int*) mxGetPr(prhs[2])); 
  segImg_t.push_back( (int*) mxGetPr(prhs[3]));

  std::vector< Scalar* > K_( nCams );
  assert( nCams*9 == mxGetNumberOfElements(prhs[5]) );
  for(int i=0;i<nCams;i++)
    K_[i] = &(((Scalar*)  mxGetPr(prhs[5]))[9*i]);

  std::vector< Scalar* > MC_( nCams-1 );
  std::vector< Scalar* > mc_( nCams-1 );
  assert( (nCams-1)*9 == mxGetNumberOfElements(prhs[6]) );
  assert( (nCams-1)*3 == mxGetNumberOfElements(prhs[7]) );
  for(int i=0;i<nCams-1;i++)
  {
    MC_[i] = &(((Scalar*)  mxGetPr(prhs[6]))[9*i]);
    mc_[i] = &(((Scalar*)  mxGetPr(prhs[7]))[3*i]);
  }

  Scalar* normals   =  (Scalar*) mxGetPr (prhs[8]);
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[9]);
  // these centers define for canonical view where to expand - can also be some area in some other view
  Scalar* centers   =  (Scalar*) mxGetPr (prhs[10]);

  int* oobSegsR  = NULL;
  int* oobSegsLT = NULL;
  int* oobSegsRT = NULL;

  proposals          =  (int*) mxGetPr(prhs[11]);// seg2prob mapping, it is actually center to mvp proposal mapping

  // here all autoscores are put in one, also for ? t0, nSegments t1, t-1 t2 ..
  Scalar* autoScore = NULL;
  autoScore = (Scalar*) mxGetPr(prhs[12]); // turn off if not needed - or set to all 1's or constant stuff

  /////////////
  if (nrhs > 13)
    halfEdgeX       =  (Scalar*) mxGetPr(prhs[13]);

  if (nrhs > 14)
    halfEdgeY       =  (Scalar*) mxGetPr(prhs[14]);

  if (nrhs > 15)
    halfEdgeXY      =  (Scalar*) mxGetPr(prhs[15]);

  if (nrhs > 16)
    halfEdgeiXY     =  (Scalar*) mxGetPr(prhs[16]);

  const int pastTimeNr = 17;
  // past time step:
  Scalar *rotations_glob(NULL);
  if (nrhs > pastTimeNr)
    pict.push_back ( (Scalar*)  mxGetPr(prhs[pastTimeNr]) );
  if (nrhs > pastTimeNr+1)
    segImg_t.push_back( (int*) mxGetPr(prhs[pastTimeNr+1]) ); 
  if (nrhs > pastTimeNr+2)
    rotations_glob      =  (Scalar*) mxGetPr(prhs[pastTimeNr+2]);

  if (nrhs > pastTimeNr+2)
    nTimeSteps++;

  Scalar* rotations_glob_t1t2(NULL);
  if (nrhs > pastTimeNr+3)
    pict.push_back ( (Scalar*)  mxGetPr(prhs[pastTimeNr+3]) );
  if (nrhs > pastTimeNr+4)
    segImg_t.push_back( (int*)  mxGetPr(prhs[pastTimeNr+4]) ); 
  if (nrhs > pastTimeNr+5)
    rotations_glob_t1t2 = (Scalar*) mxGetPr(prhs[pastTimeNr+5]);

  if (nrhs > pastTimeNr+5)
    nTimeSteps++;


  // currently i expand at every center, given by segment id
  // a proposal:
  // expansion point = center[propId%nSegments], proposal = seg2Prop[propId]
  //
  // faster: parallel expansion and graph optimization (expand 4 patches far apart) is a must here !!!
  //  std::vector< vector<int> > nSegments( nTimeSteps );// per timestep, nCam # of Segments
  std::vector< vector<Scalar*> >  imgs( nTimeSteps );// per timestep nCam image
  std::vector< vector<int*> >   segImg( nTimeSteps );// per timestep nCam image
  int maxSegs=0;
  for( int j=0;j<nTimeSteps;j++ )
    for( int i=0;i<nCams;i++ )
    {
      //      nSegments[j].push_back( getNSegments( &segImg_t[j][N*M*i], M,N ) );
      segImg[j].push_back ( (int*) malloc(N*M*sizeof(int)));
      memcpy( segImg[j][i], &segImg_t[j][N*M*i], N*M*sizeof(int));
      imgs[j].push_back   ( &pict[j][N*M*i] );
      //      maxSegs = std::max  ( nSegments[j][i], maxSegs );
    }

    std::vector< vector<Scalar*> > autoScores( nTimeSteps );// start 2 timesteps - can be more
    int temp=0;
    for( int j=0;j<nTimeSteps;j++ )
    {
      autoScores[j].resize(nCams, NULL);
      if (autoScore != NULL)
        for(int i=0;i<nCams;i++)
        {
          if ( N*M*(j*nCams+i+1) <= mxGetNumberOfElements(prhs[12]) )
            (autoScores[j])[i] = &( autoScore[N*M*(j*nCams+i)] );
        }
    }

    std::vector< vector<Scalar*> > heX( nTimeSteps );// start 2 timesteps - can be more
    std::vector< vector<Scalar*> > heY( nTimeSteps );// start 2 timesteps - can be more
    std::vector< vector<Scalar*> > heXY( nTimeSteps );// start 2 timesteps - can be more
    std::vector< vector<Scalar*> > heYX( nTimeSteps );// start 2 timesteps - can be more

    for (int i=0;i< nTimeSteps; i++)
    {
      heX[i].resize(nCams, NULL);
      heY[i].resize(nCams, NULL);
      heXY[i].resize(nCams, NULL);
      heYX[i].resize(nCams, NULL);
    }

    assert ( (nTimeSteps*nCams)*(N-1)*M     <= mxGetNumberOfElements(prhs[13]) );
    assert ( (nTimeSteps*nCams)*N*(M-1)     <= mxGetNumberOfElements(prhs[14]) );
    assert ( (nTimeSteps*nCams)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[15]) );
    assert ( (nTimeSteps*nCams)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[16]) );

    for (int j=0;j< nTimeSteps; j++)
    {
      for(int i=0;i<nCams;i++)
      {
        if ( (i+1)*(N-1)*M <= mxGetNumberOfElements(prhs[13]) )
          (heX[j])[i] = &( halfEdgeX[(i+j*nCams)*(N-1)*M] );
        if ( (i+1)*(M-1)*N <= mxGetNumberOfElements(prhs[14]) )
          (heY[j])[i] = &( halfEdgeY[(i+j*nCams)*(M-1)*N] );
        if ( (i+1)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[15]) )
          (heXY[j])[i] = &( halfEdgeXY[(i+j*nCams)*(N-1)*(M-1)] );
        if ( (i+1)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[16]) )
          (heYX[j])[i] = &( halfEdgeiXY[(i+j*nCams)*(N-1)*(M-1)] );
      }
    }

    // a grid if nothing else specified:
    // input: layout for one! timestep; copy this layout for all time steps
    // and add connections from t to t+1
    typedef Datasegments<Scalar>::P4i camTimePair;
    std::vector<camTimePair> dataCams;
    for (int j=1;j<nCams;j++)
      dataCams.push_back( camTimePair(0,0, j,0) );// cam 00 to other at t0
    for (int j=0;j<nCams;j++)
      dataCams.push_back( camTimePair(j,0, j,1) );// cam j to j at time t0 to t1
    if (nTimeSteps > 2)
      for (int j=0;j<nCams;j++)
        dataCams.push_back( camTimePair(j,0, j,2) );// cam j to j at time t0 to t -1
    if (nTimeSteps > 3)
      for (int j=0;j<nCams;j++)
        dataCams.push_back( camTimePair(j,1, j,3) );// cam j to j at time t1 to t2
    for(int i=1;i<nTimeSteps;i++)
      for (int j=1;j<nCams;j++)
        dataCams.push_back( camTimePair(0,i, j,i) );// cam 0 to j at t1, then at t-1, then at t2, ..
    //////////////////////////
    //need to know ids of dsgX; in map from .. to .. -> id
    // to select homos from dsgX

    // encoded current state:
    std::vector<camTimePair> smoothCams;// need: 00,11; 11,00 10,01;01,10
    for (int j=1;j<nCams;j++)
    {
      smoothCams.push_back( camTimePair(0,0, j,1) );//fw
      smoothCams.push_back( camTimePair(j,1, 0,0) );//bw
      smoothCams.push_back( camTimePair(j,0, 0,1) );//fw
      smoothCams.push_back( camTimePair(0,1, j,0) );//bw
    }
    // then fw:
    if (nTimeSteps > 2)
      for (int j=1;j<nCams;j++)
      {
        smoothCams.push_back( camTimePair(0,2, j,0) );
        smoothCams.push_back( camTimePair(j,2, 0,0) );
      }
      if (nTimeSteps > 3)// bw
        for (int j=1;j<nCams;j++)
        {
          smoothCams.push_back( camTimePair(0,3, j,1) );
          smoothCams.push_back( camTimePair(j,3, 0,1) );
        }
        ///////////////////////////////
        //  std::srand ( unsigned ( std::time(0) ) );	
        std::srand ( 555 ); // for debugging
        printf("rotJump %.2f \n", rotJump);
        printf("depthJump %.2f \n", depthJump);
        printf("dispMax %.2f \n", dispMax);
        printf("motMax %.2f \n", motMax);
        printf("occThresh %.2f \n", occThresh);
        printf("outOfBoundsThresh %.2f \n", outOfBoundsThresh);
        printf("data thresh %.2f \n", thresh);
        printf("motion smooth theta %.2f \n", theta);
        printf("general smoothness lambda %.3f \n", lambda);
        printf("segWeight %.2f \n", segWeight);
        printf("ot:%.3f, tempPotts:%.3f autoIntern: %.1f\n", ot, tempPotts, doAuto );

        int nSolutions = max( mxGetM(prhs[11]), mxGetN(prhs[11]) );// these are the moving planes ot go with the centers
        int nCenters   = max( mxGetM(prhs[10]), mxGetN(prhs[10]) );  // these are the expansion centers; the size must be equal to: proposals prhs[23]
        int nNormals   = max( mxGetM(prhs[ 8]), mxGetN(prhs[ 8]) );
        int nStep      = min( mxGetM(prhs[ 8]), mxGetN(prhs[ 8]) );

        int patchX = patchSize;
        int patchY = patchSize;
#ifdef _DEBUG
        patchX = patchSize;
        patchY = patchSize;
#endif

        printf("gridSize %1.f \n",  gridSize );
        printf("patchSize %1.f \n", patchSize );
        printf("nSteps     %d \n", nStep );
        printf("nNormals   %d \n", nNormals );
        printf("nCenters   %d \n", nCenters );
        printf("nSolutions %d \n", nSolutions);

//        printf("MUST BE EQUAL - nSolutions: %d, nCenters: %d \n", nSolutions, nCenters);

        // whatever, the variable is not really used
        int nSegments = 8000;//(int) (mxGetM(Segments) * mxGetN(Segments));// its a cell assigning a segment to its pixel involved

        //      printf ("nSegments %d, nProposals %d, nNormals %d \n", nSegments, nSolutions, nNormals);
        //  mexEvalString("drawnow");

        std::vector<bool> oobOutsideAll(M*N, false);

        std::vector<int> seg2Prop(nSolutions); // also just read from segMap
        for(int i=0;i<nSolutions;i++)
          if (proposals)
            seg2Prop[i] = proposals[i];//i%10;
          else
            seg2Prop[i] = i;//%10;


        ///////////////////////////////////////////////////////
        std::vector<Scalar> origins(3*nSolutions, 0); // a hack for origin centric rotations

        EvalEnergyStore<Scalar> enStore;
        enStore.prepareNormals( normals, nNormals, nStep );
        enStore.setRotTra( rotations, nNormals );

        ////////////////////////////////////////////
        // all homos to all involved other views !

        // indeed like this: star config a homo from canonical to all other views:
        std::vector< genHomoG<Scalar> > homos;
        for (int j=0; j< nTimeSteps; j++)
        {
          for(int i=0;i<nCams;i++)
          {
            // unused form 0 to 0
            if (i==0 && j==0) {homos.push_back( genHomoG<Scalar>() ); continue;}

            if (j==2 && rotations_glob != NULL)
              enStore.setGlobRotTra( rotations_glob );
            if (j==3 && rotations_glob_t1t2  != NULL)
              enStore.setGlobRotTra( rotations_glob_t1t2 );

            int id = j*nCams+i;
            if (i==0)  // cases: j!=0 i==0(timestep!=0, camera canonical)
              homos.push_back( genHomoG<Scalar> ( K_[0] ) );
            else       // cases: j!=0 i!=0(timestep!=0, camera not canonical)
              homos.push_back( genHomoG<Scalar> ( K_[0], MC_[i-1], mc_[i-1], K_[i]) );

            homos[id].setNormals( enStore.getNormals() );
            if (j==0) continue;
            homos[id].setRotTra(  enStore.getTra(), enStore.getRot() );
            homos[id].setCenters( &(origins[0]) );// could put all 0's here -> Rt w.r.t. origin

            if (j==2)
              homos[id].do_inverse_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
            if (j==3)
              homos[id].do_double_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
          }
        }

        std::vector< int >  dsgId    (nCams*nCams*nTimeSteps*nTimeSteps, -1);
        std::vector< Datasegments<Scalar>::P4i > dsgIdInv;
        std::vector< bool > dataOn   (nCams*nCams*nTimeSteps*nTimeSteps, false );
        //      std::vector< bool > smoothOn (nCams*nCams*nTimeSteps*nTimeSteps, false );// cam time pair already a dsgX

        for( int i=0;i<dataCams.size();i++)
        {
          int wtf = getQuadId( dataCams[i][0], dataCams[i][1], dataCams[i][2], dataCams[i][3], nCams, nTimeSteps);
          dataOn [  getQuadId( dataCams[i][0], dataCams[i][1], dataCams[i][2], dataCams[i][3], nCams, nTimeSteps) ] = true;
        }
        // ? start with those receiving the homos ? then 
        // use a helper vector: a> warp==#'s b> use for data term, also which needs info of which
        std::vector < Datasegments<Scalar>* > dsgX;//(nCams*nTimeSteps);
        // store ids of initialisers and data dsg's !

        std::vector < int > dsg_init;// from canonical to all other cams
        std::vector < int > dsg_data;// to call warping
        // the rest of the smooth configurations, not already somewhere else
        std::vector < int > dsg_smooth;// to call for getting homos - can not contain anything from 0,0 to somewhere - this is already covered!!!!!
        //    std::vector < Datasegments<Scalar>::P2 > dsg_initCamTime;

        for (int j=0; j<nTimeSteps; j++)
        {
          for(int i=0;i<nCams;i++)
          {
            if (i==0 && j==0) continue;

            if ( dataOn [ getQuadId( 0, 0, i, j, nCams, nTimeSteps) ] )
              dsg_data.push_back(dsgX.size());

            dsg_init.push_back(dsgX.size());
            dsgId[ getQuadId( 0,0,i,j, nCams, nTimeSteps) ] = dsgX.size();// !!! store where the config is located !!!

            // pushback copies which is hard if pointer are involved - sigh (*oobVex[ getQuadId( 0, 0, i, j, nCams, nTimeSteps)] )
            dsgX.push_back( new Datasegments<Scalar> (N, M, thresh, outOfBoundsThresh, imgs[0][0], imgs[j][i], segImg[0][0], nNormals, segImg[j][i], nNormals, oobOutsideAll,  oobOutsideAll) );

            dsgX.back()->setGridSize( gridSize );
            dsgX.back()->setTimeCam( Datasegments<Scalar>::P4i(0,0,i,j) );
            //         dsgX.back()->setCentersKs( K_[0], K_[i], centers[0][0],  centers[j][i] );
            dsgX.back()->setOccPen ( occThresh );
            //         dsg_initCamTime.push_back( Datasegments<Scalar>::P2(i,j) );
            if ( !dataOn [ getQuadId( 0, 0, i, j, nCams, nTimeSteps) ] )
              dsgX.back()->setOccBuffering(false);

            if (j==0)
              if (i >0 && mc_[i-1][0] > 0)
                dsgX.back()->setMaxDisp( -dispMax );
              else
                dsgX.back()->setMaxDisp(  dispMax );
            else
            {
              if (i==0)
                dsgX.back()->setMaxMot( motMax );// never called anyway
              else
                dsgX.back()->setMaxMot( motMax+dispMax );// never called anyway
            }

            if ( autoScores[0][0] != NULL && autoScores[j][i] != NULL )
              dsgX.back()->set_AutoScoresPix ( autoScores[0][0], autoScores[j][i]);
          }
        }

        for (int k = 0; k< dataCams.size(); k++ )
        {
          int c1 = dataCams[k][0];
          int c2 = dataCams[k][2];
          int t1 = dataCams[k][1];
          int t2 = dataCams[k][3];

          if ( dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] >=0)
            continue;
          dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] = dsgX.size();

          dsg_data.push_back(dsgX.size());  //(*oobVex[ getQuadId( c1, t1, c2, t2, nCams, nTimeSteps)] )
          dsgX.push_back( new Datasegments<Scalar> (N, M, thresh, outOfBoundsThresh, imgs[t1][c1], imgs[t2][c2], segImg[t1][c1], nNormals, segImg[t2][c2], nNormals, oobOutsideAll, oobOutsideAll) );

          //          dsgX.back()->setCentersKs( K_[c1], K_[c2], centers[t1][c1],  centers[t2][c2] );
          dsgX.back()->setGridSize( gridSize );
          dsgX.back()->setTimeCam( Datasegments<Scalar>::P4i(c1,t1,c2,t2) );
          dsgX.back()->setOccPen ( occThresh );

          if (t1==t2)
          {
            Scalar tx_c1 = c1>0 ? mc_[c1-1][0] : 0;
            Scalar tx_c2 = c2>0 ? mc_[c2-1][0] : 0;
            if ( tx_c2 < tx_c1 )
              dsgX.back()->setMaxDisp( dispMax ); // TODO: check
            else
              dsgX.back()->setMaxDisp( -dispMax ); // TODO : check
          }
          else
          {
            if (c1==c2)
              dsgX.back()->setMaxMot( motMax );// never called anyway
            else
              dsgX.back()->setMaxMot( motMax+dispMax );// never called anyway
          }

          if ( autoScores[t1][c1] != NULL &&  autoScores[t2][c2] != NULL )
            dsgX.back()->set_AutoScoresPix ( autoScores[t1][c1], autoScores[t2][c2]);
        }

        for (int k = 0; k< smoothCams.size(); k++ )
        {
          int c1 = smoothCams[k][0];
          int c2 = smoothCams[k][2];
          int t1 = smoothCams[k][1];
          int t2 = smoothCams[k][3];

          if ( dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] >=0) // done
            continue;
          dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] = dsgX.size();

          dsg_smooth.push_back(dsgX.size());  //(*oobVex[ getQuadId( c1, t1, c2, t2, nCams, nTimeSteps)] )
          dsgX.push_back( new Datasegments<Scalar> (N, M, thresh, outOfBoundsThresh, imgs[t1][c1], imgs[t2][c2], segImg[t1][c1], nNormals, segImg[t2][c2], nNormals, oobOutsideAll, oobOutsideAll) );

          // never is actually
          if ( !dataOn [ getQuadId( c1, t1, c2, t2, nCams, nTimeSteps) ] )
            dsgX.back()->setOccBuffering(false);

          //          dsgX.back()->setCentersKs( K_[c1], K_[c2], centers[t1][c1],  centers[t2][c2] );
          dsgX.back()->setGridSize( gridSize );
          dsgX.back()->setTimeCam( Datasegments<Scalar>::P4i(c1,t1,c2,t2) );
          dsgX.back()->setOccPen ( occThresh );

          if (t1==t2)
            dsgX.back()->setMaxDisp( dispMax ); // !! fails ?!
          else
          {
            if (c1==c2)
              dsgX.back()->setMaxMot( motMax );// never called anyway
            else
              dsgX.back()->setMaxMot( motMax+dispMax );// never called anyway
          }

          if ( autoScores[t1][c1] != NULL &&  autoScores[t2][c2] != NULL )
            dsgX.back()->set_AutoScores ( autoScores[t1][c1], autoScores[t2][c2]);
        }

        // for each proposal, center its (new) box after joining similar props together
        std::vector<dataElem<Scalar>::P2i> patchXY( nSolutions, dataElem<Scalar>::P2i(patchX,patchY) );
        //  std::vector<dataElem<Scalar>::P3>  propCenters;
        std::vector<Scalar>  propCenters;// store new centers here, redirect pointer to this later
        dsgX[0]->setCentersKs( K_[0], K_[0], centers,  centers, nCenters );
        //        dsgX[0]->init();// init after K's and centers are given
        dsgX[0]->reorganizeProposals ( seg2Prop, nSolutions, centers, propCenters, patchX, patchY, patchXY, false );
        centers = (Scalar*) (&(propCenters[0]));nCenters = propCenters.size()/3;

        int maxVariables = (2*patchX-5)*(2*patchY-5);
#ifdef _writeEnergies_
        printf("maxVariables before: %d - ", maxVariables);
#endif
        for (int i =0; i < patchXY.size(); i++)
          maxVariables = max( ((2*patchXY[i][0]-5)*(2*patchXY[i][1]-5) ) , maxVariables); 
#ifdef _writeEnergies_
        printf("maxVariables after : %d\n", maxVariables);
#endif
        for (int i =0;i<dsgX.size(); i++)
        {
          Datasegments<Scalar>::P4i ct_ct = dsgX[i]->getTimeCam();
          dsgX[i]->setCentersKs( K_[ct_ct[0]], K_[ct_ct[2]], centers,  centers, nCenters );
          //          dsgX[i]->init();// init after K's and centers are given
        }

        // push homographies:
        for (int i = 0;i < nNormals; i++)
        {
          for (int k = 0;k < dsg_init.size(); k++)
            dsgX[dsg_init[k]]->preInitPixel( i, &homos[dsg_init[k]+1], enStore.getNormals() );

          for (int k = 0;k < dsg_data.size(); k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
            int id1= dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
            int id2= dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
            dsgX[dsg_data[k]]->preinitPixelFromView( i, dsgX[id1]->getElement(), dsgX[id2]->getElement() );
          }
          for (int k = 0;k < dsg_smooth.size(); k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_smooth[k]]->getTimeCam();
            if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
            if (ct_ct[2] == 0 && ct_ct[3] ==0 ) continue;
            int id1 = dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
            int id2 = dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
            dsgX[dsg_smooth[k]]->preinitPixelFromView( i, dsgX[id1]->getElement(), dsgX[id2]->getElement() );
          }
        }

#pragma omp parallel for
        for (int k=0; k<dsg_data.size();k++)
        {
          Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
          dsgX[ dsg_data[k] ]->buildFromCurrentSolutionPixel( );
        }

        std::vector< EvalEnergyFullFramePixel<Scalar> > evalSmooth(smoothCams.size(), EvalEnergyFullFramePixel<Scalar> (N,M));

        for(int sc=0;sc<smoothCams.size();sc++)
        {
          camTimePair ct = smoothCams[sc];
          evalSmooth[sc].setRotJump  (rotJump);
          evalSmooth[sc].setDepthJump(depthJump);
          evalSmooth[sc].setGamma ( gamma );
          evalSmooth[sc].setRotWeight( theta );
          evalSmooth[sc].prepare( N*M, nSegments  );
          evalSmooth[sc].setSegImg( N, M, segImg[ct[1]][ct[0]], nSegments );
          evalSmooth[sc].setEdgeScale( segWeight ) ;

          // jump at 0.5 pixel -> needs to set jump smaller
          if(ct[0] == 0) // from cam 0 to .. 
            evalSmooth[sc].set2dMotionMatrix ( K_[0], MC_[std::max(0,ct[2]-1)], mc_[std::max(0,ct[2]-1)], pixJump ); 
          else  // from other cam so .. must not point to cam 0 - for me it is like this inv assume that K is equal for all cams
            evalSmooth[sc].set2dMotionMatrix_inv ( K_[ ct[0] ], MC_[ std::max(0,ct[0]-1) ], mc_[std::max(0,ct[0]-1)], pixJump );
          evalSmooth[sc].setRotTra( enStore.getRot(), enStore.getTra() );

          if( ct[0]==0 && ct[1]==0 )
            evalSmooth[sc].setNormals( enStore.getNormals() );
          else
          {
            int id = dsgId[ getQuadId( 0,0,ct[0],ct[1], nCams, nTimeSteps) ];
            evalSmooth[sc].setNormals( dsgX[id]->getvNormals() );
          }

          evalSmooth[sc].setHalfEdges(heX[ct[1]][ct[0]], heY[ct[1]][ct[0]]);
          if ( halfEdgeiXY != NULL)
            evalSmooth[sc].setCrossHalfEdges(heXY[ct[1]][ct[0]], heYX[ct[1]][ct[0]]);
          if ( halfEdgeiXY != NULL)
            evalSmooth[sc].setCrossHalfEdges(heXY[ct[1]][ct[0]], heYX[ct[1]][ct[0]]);
        }

        Scalar scaleGeo = 0.25;
        Scalar scaleMot = 1.0;
        int mvpDist = 1;// smaller

        for (int k=0; k<dsgX.size();k++)
        {
          Datasegments<Scalar>::P4i ct_ct = dsgX[k]->getTimeCam();
          dsgX[ k ]->setPotts(tempPotts);
//          dsgX[ k ]->setMVPDistance(mvpDist);// also not used
          dsgX[ k ]->setOccThresh(ot);
          if (ct_ct[0] == ct_ct[2] )
            dsgX[ k ]->setScale3D(scaleGeo);//some heavy bullshit, not used
          else
            dsgX[ k ]->setScale3D(scaleMot);//some heavy bullshit, not used
        }

#ifdef _MotionChecking_
        MotionCheck<Scalar> mChecker( NULL, nSegments, nNormals );
        mChecker.setPenalty(8. * thresh); // in dependence on data
        mChecker.setNormals( enStore.getNormals() );
        mChecker.setRotTra( enStore.getRot(), enStore.getTra() );
        mChecker.setWidthHeight(N,M);
        mChecker.setK( K_[0] );
#endif

        Scalar OLD_ENERGY = std::numeric_limits<Scalar>::max();
        Scalar ENERGY(0);

        // smaller maxValue: faster ! e.g. smoScale = 4 or so
        int maxValue = 1 << 28;// does lower values here lead to better performance in energy/overall ???
        int non_sub = 0;
        int non_sol = 0;
        int old_non_sol = 0;

        int multiplierVar = dataCams.size()*(16./4.);

        int maxNodes = multiplierVar*maxVariables; // how much larger can those become???
        int numEdges = 6*maxNodes;// 8 neighboorhood, double counting

#ifdef _NO_OPENMP
        int numThreads = 1;
#else
        // one less than the computer has threads -- otherwise system freezes
        int numThreads = min (_max_num_threads_used_, omp_get_max_threads());
#ifdef _DEBUG
        numThreads = 1;
#endif
        omp_set_num_threads( numThreads );
#endif

        printf("Threads used: %d\n", numThreads);

        std::vector< std::pair<int,int> > pidOrd;

        std::vector< GraphCutContainer<Scalar> > gs ( 4 );
        for (int i=0;i<gs.size();i++)
          gs[i].create(maxNodes, numEdges);

        for (int k=0; k<dsgX.size();k++)
          dsgX[k]->initElemStack(numThreads);

        std::clock_t fullStart(std::clock());
        // fuse all proposals in a row

        // nSegments replaced by nNprmals
        std::vector<int> order(nSolutions, 0);
        for ( int i = 0;i < nSolutions;i++ )
          order[i]=i;

        std::clock_t startP3(std::clock());
        // segments where mCheck fails
        int badCount(0), goodCount(0);

#ifdef _DEBUG
        int nIterations = 1;
#else
        int nIterations = 1;
#endif
        int maxVars = 0;
        int maxEdges = 0;

        std::vector<int> numNonSubs(N*M,0);
        int maxNodesUsed(0), maxEdgesUsed(0), nRuns(0), avNodesUsed(0), avEdgesUsed(0);

        std::vector<int> unknowns=order;
        for (int its=0; its<nIterations ; its++)
        {
          if(its > 0)
#pragma omp parallel for
            for (int k=0; k<dsg_data.size();k++)
              dsgX[ dsg_data[k] ]->buildFromCurrentSolutionPixel( );

          if ( doAuto > 0.)
            for (int k=0; k<dsg_data.size();k++)
              dsgX[ dsg_data[k] ]->updateAutoScoresPix();

#ifdef _writeEnergies_
          std::vector<Scalar> dataEnergies(dsg_data.size(),0);
          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            Scalar endE = dsgX[ dsg_data[k] ]->getEnergyPixel();
            dataEnergies[k] = endE;
            printf("Start Energy:%.2f time:%d cam:%d <-> time:%d cam:%d\n", endE, ct_ct[1], ct_ct[0], ct_ct[3], ct_ct[2]);
          }
          Scalar sumDataEMulti (std::accumulate(dataEnergies.begin(),dataEnergies.end(),0.) ), sumSmoothEMulti(0);

          sumSmoothEMulti = getSmoothnessEnergy(evalSmooth, smoothCams, dsgId, nCams, nTimeSteps, dsgX, lambda );
          printf("ENERGY: %.2f = sum(Data+smooth) :%.2f + %.2f \n", sumDataEMulti + lambda * sumSmoothEMulti, sumDataEMulti, lambda * sumSmoothEMulti);

          ENERGY = sumDataEMulti + lambda * sumSmoothEMulti;
          printf("------------ Start Energy is %.2f --------------\n", ENERGY);
          printf("Start Energies (full=data) %.2f = ( ", ENERGY);
          for (int k=0; k<dataEnergies.size();k++)
            printf("%.2f ", dataEnergies[k]);
          printf(" )\n");
#endif

          order = unknowns;
          unknowns.clear();
          // full sweep over all proposal:
          random_shuffle( order.begin(), order.end() );

          dsgX[0]->sortInitNewPixel ( order, seg2Prop );
          printf("props to investigate: %d\n", order.size());

#ifdef _DEBUG
          order.resize(50);
#endif
#ifdef _shortCut_forTesting_
          order.clear();
#endif

          std::vector<bool> processed (order.size(), false);
          int sstop(0);

#pragma omp parallel
          {
            while (!sstop)
            {
#ifdef _NO_OPENMP
              const int pid =0;
#else
              const int pid = omp_get_thread_num();
#endif

              int propId(-1), cid(-1), varNr(0), ord(0), nVars(0), nEdges(0);
              std::vector< vector< Binary<Scalar> > > allBinaries( dsg_data.size() );
              std::vector< vector< Unary<Scalar> >  > allUnaries(  dsg_data.size() );

              Scalar worstDataScore(0);
#pragma omp critical
              {
                //        printf("enter critical pid: %d", pid);//mexEvalString("drawnow");

                for(; ord < processed.size(); ord++)
                {
                  if ( !processed[ord] )
                  {
                    propId = order[ord];// propId - the id of the moving plane - this is in general not related to the center !
                    cid    = propId % nCenters;//%nSegments

                    bool okToRun=true;
                    for( int k=0;k<dsg_init.size();k++)
                      okToRun = okToRun && dsgX[dsg_init[k]]->testNewPixel( cid, seg2Prop[propId], patchXY[propId][0], patchXY[propId][1] );
                    if (okToRun) break;// okToRun == true now

                  }
                }

                // end it now
                if (cid == -1)
                {
                  sstop = 1;
#pragma omp flush(sstop)
                }

                if ( pid >= numThreads)
                {
                  printf("Warning: process id:%d > number of threads:%d \n", pid, numThreads);
                  ord = order.size();
#ifdef WIN32
                  Sleep( 20 );
#endif
                  varNr =0;
                }

                if (ord < order.size() )
                {
                  processed[ord] =true;

                  // divided regularly
                  // printf("Thread: %d/%d on %d \n", pid, numThreads, ord);
                  pidOrd.push_back (std::pair<int,int> (pid, ord) );

                  varNr=0;
                  for (int k=0;k<dsg_init.size();k++)
                    varNr = dsgX[dsg_init[k]]->initNewPixel( cid, seg2Prop[propId], patchXY[propId][0], patchXY[propId][1], 1, varNr );
                  varNr += dsgX[dsg_init[0]]->getElemVarsPixel(); // + elements in first box

                  //                printf("varNr:%d ", varNr);
                }
                if (varNr > 0)
                {
                  for (int k = 0;k < dsg_data.size(); k++)
                  {
                    Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
                    if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
                    int id1= dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
                    int id2= dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
                    dsgX[dsg_data[k]]->initNewPixelFromView( cid, seg2Prop[propId], dsgX[id1]->getElement(), dsgX[id2]->getElement());
                  }
                  for (int k = 0;k < dsg_smooth.size(); k++)
                  {
                    Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_smooth[k]]->getTimeCam();
                    if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
                    int id1= dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
                    int id2= dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
                    dsgX[dsg_smooth[k]]->initNewPixelFromView( cid, seg2Prop[propId], dsgX[id1]->getElement(), dsgX[id2]->getElement());
                  }
                  for (int k = 0;k < dsg_data.size(); k++)
                  {
                    Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
                    dsgX[dsg_data[k]]->createWarpsPixel( );
                    dsgX[dsg_data[k]]->expandRegionPixel( seg2Prop[propId] );

                    allBinaries[k] = dsgX[dsg_data[k]]->getBinaries();
                    allUnaries[k]  = dsgX[dsg_data[k]]->getUnaries();

                    for (int kx=0; kx< allUnaries[k].size();kx++)
                      worstDataScore += max( allUnaries[k][kx].node[0], allUnaries[k][kx].node[1] );
                    for (int kx=0; kx< allBinaries[k].size();kx++)
                      worstDataScore += max(max( max( allBinaries[k][kx].edge[0], allBinaries[k][kx].edge[1] ),allBinaries[k][kx].edge[2]), allBinaries[k][kx].edge[3] );
                  }

                  for (int k = 0;k < dsgX.size(); k++)
                    dsgX[k]->storeElem(pid);

                  nVars = varNr ;//expPropVars[ propId ];
                }//var > 0
              }// omp critical

              // get rid of this critical here !! done
              //#pragma omp critical
              if (nVars > 0)
              {
                std::vector<Scalar> F00;
                std::vector<Scalar> F01;
                std::vector<Scalar> F10;
                std::vector<Scalar> F11;
                std::vector<Scalar> F0(nVars,0);
                std::vector<Scalar> F1(nVars,0);
                std::vector<int> Idk;
                std::vector<int> Idl;

                // now that is an art:
                Scalar worstSmoothScore(0);

                for(int sc=0;sc<evalSmooth.size();sc++)
                {
                  camTimePair ct = smoothCams[sc];
                  // a) find the boxes:
                  //
                  int qid_10 = getQuadId( ct[0], ct[1], ct[0], ct[3], nCams, nTimeSteps);
                  int qid_01 = getQuadId( ct[0], ct[1], ct[2], ct[1], nCams, nTimeSteps);
                  int qid_11 = getQuadId( ct[0], ct[1], ct[2], ct[3], nCams, nTimeSteps);

                  int qid10 = dsgId[ qid_10 ];// saem cam other time
                  int qid01 = dsgId[ qid_01 ];// same time other cam
                  int qid11 = dsgId[ qid_11 ];
                  // b) get whether fw or bw direction of dsg
                  //

                  Datasegments<Scalar>::P4i dummy1 = getQuadIdInv( qid_01, nCams, nTimeSteps );
                  Datasegments<Scalar>::P4i dummy2 = getQuadIdInv( qid_10, nCams, nTimeSteps );
                  Datasegments<Scalar>::P4i dummy3 = getQuadIdInv( qid_11, nCams, nTimeSteps );
                  bool i01 = getQuadIdInv( qid_01, nCams, nTimeSteps ) == Datasegments<Scalar>::P4i( ct[0], ct[1], ct[2], ct[1]);
                  bool i10 = getQuadIdInv( qid_10, nCams, nTimeSteps ) == Datasegments<Scalar>::P4i( ct[0], ct[1], ct[0], ct[3]);
                  bool i11 = getQuadIdInv( qid_11, nCams, nTimeSteps ) == Datasegments<Scalar>::P4i( ct[0], ct[1], ct[2], ct[3]);

                  dataElem<Scalar> de10 = dsgX[qid10]->getStoredElem(pid);
                  dataElem<Scalar> de01 = dsgX[qid01]->getStoredElem(pid);
                  dataElem<Scalar> de11 = dsgX[qid11]->getStoredElem(pid);

                  const std::vector<int>& loc2glob = i01 ? de01.seg2varFW : de01.seg2varBW;

                  int maxVar=0;
                  for (int kk=0;kk<loc2glob.size();kk++)
                    maxVar = std::max( maxVar, loc2glob[kk]);
                  assert (maxVar <= F0.size() );

                  worstSmoothScore += evalSmooth[sc].compute_score_combiDepth_consistent( seg2Prop[propId], 
                    i01 ? de01.boxFW : de01.boxBW, i01 ? de01.seg2varFW : de01.seg2varBW,// i01 ? de01.varsFW : de01.varsBW, 
                    i01 ? dsgX[qid01]->getvHoms() : dsgX[qid01]->getviHoms(), 
                    i10 ? dsgX[qid10]->getvHoms() : dsgX[qid10]->getviHoms(), 
                    i11 ? dsgX[qid11]->getvHoms() : dsgX[qid11]->getviHoms(), 
                    F00, F01, F10, F11, F0, F1, Idk, Idl );
                }

                int nEdgesSmooth =0;
                nEdgesSmooth += F00.size();
                nEdges = nEdgesSmooth;
                for (int kk=0;kk<allBinaries.size();kk++)
                  nEdges += allBinaries[kk].size();

                maxNodesUsed = max( maxNodesUsed, nVars);
                maxEdgesUsed = max( maxEdgesUsed, nEdges);
                nRuns ++;
                avNodesUsed += nVars;
                avEdgesUsed += nEdges;


                assert (maxNodes >= nVars && numEdges >= nEdges );
                if (maxNodes < nVars || numEdges < nEdges ) // insufficient memory:
                {
                  printf("\nWarning violation: %5d,%5d edges: %6d,%6d: \n\n\n\n", maxNodes, nVars, numEdges, nEdges );
                  //mexEvalString("drawnow");
                }

                gs[pid].setup(nVars);

#ifdef _MotionChecking_
                dataElem<Scalar> de1 = dsgX[0]->getStoredElem(pid);//->getElement ();
                std::vector< Datasegments<Scalar>::P2 > mScoresUn( de1.varsFW, Datasegments<Scalar>::P2(0.,0.) );
                mChecker.computeGeometryConsistencyPerPixelLocal( N, M, seg2Prop[propId], de1.boxFW, de1.seg2varFW, segImg[0][0], dsgX[0]->getvNormals( ), mScoresUn, MotionCheck<Scalar>::P3(mc_[0]), dispMax );
#endif

                Scalar scale = (Scalar)(maxValue) / (lambda * worstSmoothScore + worstDataScore);
                Scalar scaleSmo = scale * lambda;

#ifdef _MotionChecking_
                for (int k=0; k< mScoresUn.size();k++)
                  gs[pid].AddUnaryTerm( k, scale*(mScoresUn[k])[0], scale*(mScoresUn[k])[1] );
#endif

                for (int k=0;k<allUnaries.size(); k++)
                  addUnariesBinaries( (gs[pid]),  allUnaries[k],  allBinaries[k], scale, non_sub );

                /// smoothness:
                addSmoothnessConstraints( (gs[pid]), F0, F1, F00, F01, F10, F11, Idk, Idl, 0, F0.size(), scaleSmo, non_sub, nVars);
                ////////////////////////
                if (pid >= numThreads) // THOU SHALL NOT PASS
                {
                  printf("Warning: pid: %d not in range of # numThreads: %d", pid, numThreads);
                  varNr = 0;
                }

                gs[pid].finalize();
                gs[pid].solve();

#ifdef _qpboVersion_				  
                // NEW: improve moves:
                gs[pid].resolve(false, false);
#endif		  

                std::vector<int> solution01 ( nVars, 0 );
                std::vector<int> noSolution ( nVars, 0 );
                for (int j = 0; j < nVars;j++)
                {
                  int x = (gs[pid]).GetLabel(j);
                  if (x<0)
                  {          
                    noSolution[j]=1;
                    non_sol++;
                  }

                  if (x==1)
                    solution01[j] = 1;
                }

#ifdef check_correctness
                Scalar diff = dsgX[0]->getEnergyPixelLocal( solution01, propId ) ;
                if (diff > 0)
                  printf("energy diff is %7.2f\n", diff);
#endif

                for (int k=0; k<dsg_data.size();k++)
                  dsgX[dsg_data[k]]->updateSolutionsPixel( seg2Prop[propId], solution01, pid );
                for (int k=0; k<dsgX.size();k++)
                  dsgX[k]->invalidateElem( pid );

                if ( (non_sol - old_non_sol) > nVars/40 )
                  unknowns.push_back(propId);
                old_non_sol = non_sol ;

                if (maxVars  < nVars)  maxVars  = nVars;
                if (maxEdges < nEdges) maxEdges = nEdges;

                ///////////////
#ifdef __Energy_Check__Data__
                Scalar ENERGY2 = dsg.getEnergyPixel( );
                if (ENERGYOLD < ENERGY2)
                  printf("propId %4d ns %4d, Energies new, gain: %08.2f, %05.2f \n",  propId, non_sol - old_non_sol, ENERGY2, ENERGYOLD-ENERGY2 );
                mexEvalString("drawnow");
                ENERGYOLD = ENERGY2;
                old_non_sol = non_sol ;
#endif
                (gs[pid]).reset();
              }//nVars>0
            }//ompcritical but should be while ! stop
            // once only:
#ifdef _DEBUG
            OLD_ENERGY = ENERGY;
#endif
          }//omp parallel - but it is while !stop
        }// nits
        std::clock_t stopP3(std::clock());
        printf("QPBO took %.2f\n", double(stopP3-startP3)/CLOCKS_PER_SEC );

#ifdef _writeEnergies_
        printf("ords processed : %d\n", pidOrd.size() );

#pragma omp parallel for
        for (int k=0; k<dsg_data.size();k++)
        {
          Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
          dsgX[ dsg_data[k] ]->buildFromCurrentSolutionPixel( );
        }
#endif

#ifdef _writeEnergies_
        std::vector<Scalar> dataEnergiesNew( dsg_data.size(), 0 );
        for (int k=0; k<dsg_data.size(); k++)
        {
          Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
          Scalar endE = dsgX[ dsg_data[k] ]->getEnergyPixel();
          dataEnergiesNew[k] = endE;
          printf("End Energy:%.2f time:%d cam:%d <-> time:%d cam:%d\n", endE, ct_ct[1], ct_ct[0], ct_ct[3], ct_ct[2]);
        }
        Scalar sumDataEMulti = std::accumulate(dataEnergiesNew.begin(),dataEnergiesNew.end(),0.);
        printf("ENERGY: %.2f  NON submodularity detected %d (%.2f pct), no solutions: %d (%.2f pct); \n", sumDataEMulti, non_sub, double(100*non_sub)/avEdgesUsed, non_sol, double(100*non_sol)/avNodesUsed);
        printf("av/max Nodes :%d/%d, av/max Edges: %d/%d \n", avNodesUsed/nRuns, maxNodesUsed, avEdgesUsed/nRuns, maxEdgesUsed);
        Scalar allSmoothScore = getSmoothnessEnergy(evalSmooth, smoothCams, dsgId, nCams, nTimeSteps, dsgX, lambda );
        printf("ENERGY: %.2f = sum(Data+smooth) :%.2f + %.2f \n", sumDataEMulti + lambda * allSmoothScore, sumDataEMulti, lambda * allSmoothScore);
#endif

        // debugging - i need the following outputs:
        // where are the imppenalties? and where the non-subs?!
        if (nlhs > 3)
        {
#ifndef _writeEnergies_
          Scalar allSmoothScore = getSmoothnessEnergy(evalSmooth, smoothCams, dsgId, nCams, nTimeSteps, dsgX, lambda );
          std::vector<Scalar> dataEnergiesNew( dsg_data.size(), 0 );
          for (int k=0; k<dsg_data.size(); k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            Scalar endE = dsgX[ dsg_data[k] ]->getEnergyPixel();
            dataEnergiesNew[k] = endE;
          }
          Scalar sumDataEMulti = std::accumulate(dataEnergiesNew.begin(),dataEnergiesNew.end(),0.);
#endif
         mwSize dims[2];
         dims[0] = 1;
         dims[1] = 6;
         plhs[3] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
         Scalar* outputI = (Scalar*) mxGetPr(plhs[3]);
         outputI[0] = sumDataEMulti + lambda * allSmoothScore; 
         outputI[1] = non_sub; outputI[2] = non_sol; 
         outputI[3] = avNodesUsed; outputI[4] = avEdgesUsed; outputI[5] = nRuns;
        }

        // -- put all in one !!!
        if (nlhs > 1)
        {
          // new:
          std::vector<Scalar> depthMap2(N*M);
          mwSize dims[3];
          dims[0] = M;
          dims[1] = N;
          dims[2] = nCams*nTimeSteps;//dsg_data.size();
          plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
          Scalar* outputI = (Scalar*) mxGetPr(plhs[1]);
          memset( outputI, 0, sizeof(double)*dims[2]*dims[1]*dims[0] );

          for (int k=0;k<dsg_data.size();k++)
          {
            // should be used to sort things out
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            int id1 = ct_ct[0] + ct_ct[1] * nCams;
            int id2 = ct_ct[2] + ct_ct[3] * nCams;

            if (outputI[id1 * N*M] == 0)
            {
              dsgX[ dsg_data[k] ]->buildDepthMap( depthMap2 );
              for (int i=0;i < N*M; i++)
                outputI[i + id1*N*M] = depthMap2[i];
            }
            if (outputI[id2 * N*M] == 0)
            {
              dsgX[ dsg_data[k] ]->buildDepthMap2( depthMap2 );
              for (int i=0;i < N*M; i++)
                outputI[i + id2 * N*M] = depthMap2[i];
            }
          }
        }
        if (nlhs > 2)
        {
          // new
          std::vector<int> dataEnergyMap(N*M);
          std::vector<int> dataEnergyMap2(N*M);
          mwSize dims[3];
          dims[0] = M;
          dims[1] = N;
          dims[2] = dsg_data.size()*2;
          plhs[2] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
          int* outputI = (int*) mxGetPr(plhs[2]);

          for (int k=0;k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            dsgX[ dsg_data[k] ]->getEnergyMapPixel( dataEnergyMap, dataEnergyMap2 );
#ifdef _writeEnergies_
            printf("Data plot order: ct_ct: %d, %d -> %d, %d \n", ct_ct[0], ct_ct[1], ct_ct[2], ct_ct[3] );
#endif
            // (segImg[ct_ct[1]][ct_ct[0]])
            for (int i=0;i < N*M; i++)
              outputI[i+    2*k*N*M] = dataEnergyMap [i];
            for (int i=0;i < N*M; i++)
              outputI[i+N*M+2*k*N*M] = dataEnergyMap2[i];
          }
        }

        if (nlhs > 0)
        {
          mwSize dims[3];
          dims[0] = M;
          dims[1] = N;
          dims[2] = nCams*nTimeSteps;
          plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
          Scalar * output  = (Scalar*) mxGetPr(plhs[0]);
          for( int k=0;k<nTimeSteps;k++ )
            for(int j=0;j<nCams;j++ )
            {
              for ( int i = 0; i < N*M; i++ )
                output[ i + (k*nCams+j)*N*M] = (Scalar) (segImg[k][j] [i]);
              free( segImg[k][j] );
            }
        }

        for( int i=0;i<dsgX.size();i++)
          delete dsgX[i];
}
#endif // __VCSF_PIX_GC__
