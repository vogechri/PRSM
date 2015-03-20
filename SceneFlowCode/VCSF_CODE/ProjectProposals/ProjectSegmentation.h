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

#ifndef __PROJECTSEGMENTATION__
#define __PROJECTSEGMENTATION__

#include <stdio.h>
#include <string.h>

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include "EvalEnergyStore.h"
#include "genHom.h"
#include "OcclusionMappingBufferVC.h"

#include <algorithm>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

#include <ctime>

using namespace std;
using namespace Math;

// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void ProjectSegmentation( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

  size_t N,M;
  if (nrhs < 6)
    mexErrMsgTxt("Must have at least 6 input arguments");

  int* segImg(NULL), *origSegImg(NULL);
  Scalar* Kr_(NULL);

  Scalar* p2d_    = (Scalar*)  mxGetPr(prhs[0]); // K_l^-1
  Scalar* Kl_     = (Scalar*)  mxGetPr(prhs[1]); // K_l

  // optional: always 0,0,0 if not needed
  Scalar* mc_       =  (Scalar*) mxGetPr (prhs[2]);
  Scalar* MC_       =  (Scalar*) mxGetPr (prhs[3]);

  const mxArray* Segments  = (prhs[6]);
  Scalar* normals   = (Scalar*) mxGetPr(prhs[ 7]);
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[8]);
  Scalar* centers   =  (Scalar*) mxGetPr (prhs[9]);

  Scalar* centers2  =  (Scalar*) mxGetPr (prhs[10]);

//  const mxArray* edges     = (prhs[11]);
  const mxArray* edges2    = (prhs[12]);

       origSegImg       =  (int*) mxGetPr(prhs[13]);
  int* origSegImg2      =  (int*) mxGetPr(prhs[14]);

  int* proposalSolution =  (int*) mxGetPr(prhs[15]);// as many as segments map from segments to proposal set

  if (nrhs > 16)
    Kr_             = (Scalar*)  mxGetPr(prhs[16]); // K_r
  else
    Kr_ = Kl_;

  int view;
  if (nrhs > 17)
    view             = (Scalar) *mxGetPr(prhs[17]); // view
  else
    view = 0;

  int secondView    = -1;
  if (nrhs > 18)
    secondView      = (Scalar) *mxGetPr(prhs[18]); // 2nd view

  int* proposalSolution2 = NULL;
  if (nrhs > 19)
    proposalSolution2 =  (int*) mxGetPr(prhs[19]);// as many as segments map from segments to proposal set


  int thirdView    = -1;
  if (nrhs > 20)
    thirdView      = (Scalar) *mxGetPr(prhs[20]); // 2nd view

  int* proposalSolution3 = NULL;
  if (nrhs > 21)
    proposalSolution3 =  (int*) mxGetPr(prhs[21]);// as many as segments map from segments to proposal set

  Scalar* rotations_glob(NULL);
  if (nrhs > 22)
    rotations_glob =  (Scalar*) mxGetPr (prhs[22]);

  //////////////////////////

  Scalar* output, *output2 = NULL;
  Scalar* output3, *output4 = NULL;
  Scalar occThresh = 1.0; // not 0 :)

  M           = (size_t) mxGetM(prhs[13]);
  N           = (size_t) mxGetN(prhs[13]);

  printf("View: %d\n", view);

  int nNormals = max( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );
  int nStep    = min( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );
  int nSegments  = (int) (mxGetM(prhs[15]) * mxGetN(prhs[15]));
  int nSegments2 = (int) (mxGetM(edges2) * mxGetN(edges2));

  mwSize dims[2];dims[0] = N;dims[1] = M;
  // Create the output array
  {
    // segment map later
    plhs[0]      =  mxCreateDoubleMatrix(1, nSegments2, mxREAL);
    output  = (Scalar*) mxGetPr(plhs[0]);
    memset(output, 0, sizeof(Scalar) * nSegments2);
  }

  printf ("nSegments %d, nNormals %d \n", nSegments, nNormals);
  ///////////////////////////////////////////////////////

  std::vector<int> solution(nSegments);
  for (int i=0;i< nSegments; i++)
    solution[i] = proposalSolution[i];

  std::vector<int> solution2;// that the output

  std::vector<int> solutionB ( nSegments2 );
  if ( proposalSolution2 != NULL )
  {
    for (int i=0;i< nSegments2; i++)
      solutionB[i] = proposalSolution2[i];
  }
  std::vector<int> solutionC ( nSegments2 );
  if ( proposalSolution3 != NULL )
  {
    for (int i=0;i< nSegments2; i++)
      solutionC[i] = proposalSolution3[i];
  }

  EvalEnergyStore<Scalar> enStore;
  enStore.prepareNormals( normals, nNormals, nStep );
  enStore.setRotTra( rotations, nNormals );

  if (rotations_glob  != NULL)
    enStore.setGlobRotTra( rotations_glob );

  genHomoG<Scalar> gHom1 (Kr_, MC_, mc_);// t,t, left right
  genHomoG<Scalar> gHom2 (Kl_);          // t,t+1 left
  genHomoG<Scalar> gHom3 (Kr_, MC_, mc_);// t,t+1 left right

  gHom1.setNormals( enStore.getNormals() );
  gHom2.setNormals( enStore.getNormals() );
  gHom2.setRotTra(  enStore.getTra(), enStore.getRot() );
  
  gHom3.setNormals( enStore.getNormals() );
  gHom3.setRotTra(  enStore.getTra(), enStore.getRot() );

//  gHom1.setCenters( centers ); // not needed and used to distinguish
  gHom2.setCenters( centers );
  gHom3.setCenters( centers );
  //

  genHomoG<Scalar> gHom1Li (Kl_);
  genHomoG<Scalar> gHom1Ri (Kr_, MC_, mc_);// l t,r t-1 left right

  gHom1Li.setNormals( enStore.getNormals() );
  gHom1Li.setRotTra(  enStore.getTra(), enStore.getRot() );
  gHom1Li.setCenters( centers );// could put all 0's here -> Rt w.r.t. origin
  gHom1Ri.setNormals( enStore.getNormals() );
  gHom1Ri.setRotTra(  enStore.getTra(), enStore.getRot() );
  gHom1Ri.setCenters( centers );// could put all 0's here -> Rt w.r.t. origin
  if (rotations_glob  != NULL)
  {
    gHom1Li.do_inverse_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
    gHom1Ri.do_inverse_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
  }

  int patchSize = 25;
  bool doOcclusionStuff =false;
  OcclusionMapBufferVC<Scalar> testOPix((2*patchSize+1)*(2*patchSize+1),  doOcclusionStuff);

  testOPix.SetOccPenalty(occThresh);
  testOPix.setSegImg( N,M, origSegImg, (int) (mxGetM(Segments) * mxGetN(Segments)) );
  testOPix.setP2d( p2d_ );

 if ( proposalSolution2 != NULL && secondView>-1 && secondView < 4)
 {
   if ( secondView ==0)
     testOPix.setSecondTestView( &gHom1, solutionB);
   if ( secondView ==1)
     testOPix.setSecondTestView( &gHom2, solutionB);
   if ( secondView ==2)
     testOPix.setSecondTestView( &gHom3, solutionB);
   if ( secondView ==3)
     testOPix.setSecondTestView( &gHom1Li, solutionB);
   if ( secondView ==4)
     testOPix.setSecondTestView( &gHom1Ri, solutionB);
 }

 if ( proposalSolution3 != NULL && thirdView>-1 && thirdView < 4)
 {
   if ( thirdView ==0)
     testOPix.setThirdTestView( &gHom1, solutionC);
   if ( thirdView ==1)
     testOPix.setThirdTestView( &gHom2, solutionC);
   if ( thirdView ==2)
     testOPix.setThirdTestView( &gHom3, solutionC);
   if ( secondView ==3)
     testOPix.setThirdTestView( &gHom1Li, solutionC);
   if ( secondView ==4)
     testOPix.setThirdTestView( &gHom1Ri, solutionC);
 }

  if (view==0)
    testOPix.extendSegmentation( &gHom1, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kr_  );
  if (view==1)
    testOPix.extendSegmentation( &gHom2, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kl_  );
  if (view==2)
    testOPix.extendSegmentation( &gHom3, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kr_  );
  if (view==3)
    testOPix.extendSegmentation( &gHom1Li, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kl_  );
  if (view==4)
    testOPix.extendSegmentation( &gHom1Ri, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kr_  );

  if (view!=0 && view!=1 && view!=2 && view!=3 && view!=4 )
    testOPix.extendSegmentation( &gHom1, enStore.getNormals(), edges2, centers, nSegments, centers2, solution, origSegImg, origSegImg2, solution2, nSegments2, Kl_, Kr_  );


  for (int i = 0 ; i < nSegments2; i++)
    output[i] = solution2[i];
}
#endif// __PROJECTSEGMENTATION__