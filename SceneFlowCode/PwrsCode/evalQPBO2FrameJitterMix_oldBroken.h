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

/// here we allow proposals to be moved to their neighbours
/// the selection is taken randomly, we pick a segment, distribute the motion
/// to the neighbours and fix these
/// finally we stop if no segment can be picked any more
/// we do this N times
/// also possible Growing seeds like in patchMatch

// idea to do the mixing: there are n segments
// randomly we generate n new proposals, combining normal and motion from neighbors
// this delivers the new proposal indices N+1..2N
// these are evaluated as before, interstingly the ids are fix:
// 0: 1:N , 1: N+1..2N
// the winners are moved into the first N places 
// REPEAT until satisfied
#ifndef __QPBO_JIT_2FRAME__
#define __QPBO_JIT_2FRAME__
////////////////////////////////////////////////////////////////
// Combine existing moving planes fom neighbor segments  //
////////////////////////////////////////////////////////////////
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "AccumData.h"
#include "DataDefinitions.h"

#include "EvalEnergyFull2Frame.h"

#include "EvalEnergyStore.h"

#include <algorithm>
#include <list>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T
#include "../QPBO-v1.3.src/QPBO.h"

#include <ctime>

using namespace std;
using namespace Math;

#define __Small_Improvement__ 0.501

// also use cross view
//#define _use_5th_view

// shrink the store to nSegments entries
template<typename Scalar>
void getStoreConsistency(int nSegments, std::vector<int>& currentSolution, EvalEnergyStore<Scalar>& enStore )
{
  std::vector<typename EvalEnergyStore<Scalar>::P3>* norStored = enStore.manipulateNormals();
  std::vector<typename EvalEnergyStore<Scalar>::P3>* traStored = enStore.manipulateTra( );
  std::vector<typename EvalEnergyStore<Scalar>::M3>* rotStored = enStore.manipulateRot( );

  std::vector<int> valids(norStored->size(), 0);
  for(int i=0;i< currentSolution.size(); i++)
    valids[currentSolution[i]] = 1;

  std::vector< int > replacements(norStored->size(), 0);

  int firstIn = 0;
  for(int i=nSegments;i < valids.size(); i++)
    if (valids[i] ) // actually used: replace all currentsolution entries == i by firstIn
    {
      while ( firstIn<currentSolution.size() && valids[firstIn] ) 
        firstIn++;
      if (firstIn >  currentSolution.size())
        break;
      // replace:
      (*norStored)[firstIn] = (*norStored)[i];
      (*traStored)[firstIn] = (*traStored)[i];
      (*rotStored)[firstIn] = (*rotStored)[i];

      replacements[i] = firstIn;
      valids[firstIn] = 1;
    }

  norStored->resize(nSegments);
  traStored->resize(nSegments);
  rotStored->resize(nSegments);

  // replace entries in current solution
  for(int i=0;i< currentSolution.size(); i++)
    if (currentSolution[i] >= nSegments) //replace
    {
      currentSolution[i] = replacements[ currentSolution[i] ];
    }

}

//////////////////////////////////////////////////////////////////////////////////////
// now use preselected stuff : segments which were turned in the last iteration should be further expanded - right, so do them first
// newFlipped are the parts where a change occured, so more likely to be cahnged again, process these first
// LATER: use sampling? varying the current solution just a bit ? or gradient descent on ..
template<typename Scalar>
std::vector<int> generateNeighProposals(int nSegments, const mxArray* edges, const std::vector<int>& currentSolution, 
                                        std::vector<int>& newFlipped, EvalEnergyStore<Scalar>& enStore )
{
  std::vector<int> propIds(nSegments, 0);
  std::vector<int> valids (nSegments, 1);
  std::vector<int> result (nSegments, 0);

  std::vector<typename EvalEnergyStore<Scalar>::P3>* norStored = enStore.manipulateNormals();
  std::vector<typename EvalEnergyStore<Scalar>::P3>* traStored = enStore.manipulateTra( );
  std::vector<typename EvalEnergyStore<Scalar>::M3>* rotStored = enStore.manipulateRot( );

  norStored->resize(nSegments);
  traStored->resize(nSegments);
  rotStored->resize(nSegments);

  for (int i=0;i<nSegments;i++)
  {
    propIds[i] = i;
    result [i] = currentSolution[i];
  }
  ///////////
  std::random_shuffle( propIds.begin(), propIds.end() );

  int fsz = newFlipped.size();
  if (fsz > 0)
  {
    std::random_shuffle( newFlipped.begin(), newFlipped.end() );
    newFlipped.resize( propIds.size() + fsz );
    std::copy( propIds.begin(), propIds.end(), newFlipped.begin() + fsz);
    propIds = newFlipped;
  }

  // first come first served
  int countNew = nSegments;
  for (int i = 0; i < propIds.size() ;i++)
  {
    int segi = propIds[i];

    if (!valids[segi])
      continue;

    int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges, segi) )/5;
    Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges, segi) );

    // use a different proposals for the middle (not the same anyway)
    // a little more randomized by setting this also different from the current solution
    result[segi] = currentSolution[(edge [5*0 ])-1];
    valids[segi] = 0;

    int genNew = rand() % 3;
    switch (genNew)
    {
    case 1: // here simply keep the normal, xchange the rotation
      {
      typename EvalEnergyStore<Scalar>::M3 rotation    = (*rotStored)[ currentSolution[segi] ];
      typename EvalEnergyStore<Scalar>::P3 translation = (*traStored)[ currentSolution[segi] ];

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
          norStored->push_back( (*norStored)[currentSolution[id]] );
          traStored->push_back( translation );
          rotStored->push_back( rotation );
          result[id] = countNew++;
          valids[id] = 0;
        }
      }
      break;
      }
    case 2: // here simply keep the rotation, xchange the normal
      {
      typename EvalEnergyStore<Scalar>::P3 normal = (*norStored)[ currentSolution[segi] ];

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
          norStored->push_back( normal );
          traStored->push_back( (*traStored)[currentSolution[id]] );
          rotStored->push_back( (*rotStored)[currentSolution[id]] );
          result[id] = countNew++;
          valids[id] = 0;
        }
      }
      break;
      }

    default:
    case 0:
      {
      /// distribute current proposal to neighbours: grow
      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
         result[id] = currentSolution[segi];
         valids[id] = 0;
        }
      }
      break;
      }
    }
  }
  return result;
}


template<typename Scalar>
std::vector<int> generateNeighProposals_refit(int nSegments, const mxArray* edges, const std::vector<int>& currentSolution, 
                                              std::vector<int>& newFlipped, EvalEnergyStore<Scalar>& enStore, Scalar* centers, Scalar* K_ )
{
  std::vector<int> propIds(nSegments, 0);
  std::vector<int> valids (nSegments, 1);
  std::vector<int> result (nSegments, 0);

  typename EvalEnergyStore<Scalar>::M3 K(K_);

  std::vector<typename EvalEnergyStore<Scalar>::P3>* norStored = enStore.manipulateNormals();
  std::vector<typename EvalEnergyStore<Scalar>::P3>* traStored = enStore.manipulateTra( );
  std::vector<typename EvalEnergyStore<Scalar>::M3>* rotStored = enStore.manipulateRot( );

  norStored->resize(nSegments);
  traStored->resize(nSegments);
  rotStored->resize(nSegments);

  for (int i=0;i<nSegments;i++)
  {
    propIds[i] = i;
    result [i] = currentSolution[i];
  }
  ///////////
  std::random_shuffle( propIds.begin(), propIds.end() );

  int fsz = newFlipped.size();
  if (fsz > 0)
  {
    std::random_shuffle( newFlipped.begin(), newFlipped.end() );
    newFlipped.resize( propIds.size() + fsz );
    std::copy( propIds.begin(), propIds.end(), newFlipped.begin() + fsz);
    propIds = newFlipped;
  }

  // first come first served
  int countNew = nSegments;
  for (int i = 0; i < propIds.size() ;i++)
  {
    int segi = propIds[i];

    if (!valids[segi])
      continue;

    int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges, segi) )/5;
    Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges, segi) );

    // use a different proposals for the middle (not the same anyway)
    // a little more randomized by setting this also different from the current solution
    result[segi] = currentSolution[(edge [5*0 ])-1];
    valids[segi] = 0;

    int genNew = rand() % 4;//3 for testing
    switch (genNew)
    {
     case 3: // here simply keep the normal, xchange the rotation
      {
      typename EvalEnergyStore<Scalar>::M3 rotation    = (*rotStored)[ currentSolution[segi] ];
      typename EvalEnergyStore<Scalar>::P3 translation = (*traStored)[ currentSolution[segi] ];

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
          norStored->push_back( (*norStored)[currentSolution[id]] );
          traStored->push_back( translation );
          rotStored->push_back( rotation );
          result[id] = countNew++;
          valids[id] = 0;
        }
      }
      break;
      }
    case 1: // here refit the normal, xchange the rotation
      {
      typename EvalEnergyStore<Scalar>::M3 rotation    = (*rotStored)[ currentSolution[segi] ];
      typename EvalEnergyStore<Scalar>::P3 translation = (*traStored)[ currentSolution[segi] ];

      typename EvalEnergyStore<Scalar>::P3 kt       = K * translation; // from new

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
          typename EvalEnergyStore<Scalar>::P3 center ( &centers[3 * id] );
          typename EvalEnergyStore<Scalar>::P3 n_start  = (*norStored)[currentSolution[id]];

          typename EvalEnergyStore<Scalar>::P3 c_3d     = center / (center|n_start);
          Scalar dStart = c_3d[2];
          if ( dStart != dStart )
          {
            printf ("Fail generateNeighProposals_refit %f \n", dStart);
//            mexEvalString("drawnow");
          }

          typename EvalEnergyStore<Scalar>::P3 c_2d     = K * c_3d;  c_2d /= dStart;
          /// 
          typename EvalEnergyStore<Scalar>::P3 c_2dG    = K * ( (*traStored)[ currentSolution[id] ] + c_3d );c_2dG /= c_2dG[2];
          if ( c_2dG[2] != c_2dG[2] )
          {
            printf ("Fail generateNeighProposals_refit %f \n", c_2dG[2]);
//            mexEvalString("drawnow");
          }

          typename EvalEnergyStore<Scalar>::P3 u = c_2d-c_2dG;
          typename EvalEnergyStore<Scalar>::P3 v = kt - c_2dG * kt[2];// can be 0
          Scalar s = -(u[0]*v[0] + u[1]*v[1]) / (v[0]*v[0]+v[1]*v[1]);
          if ( fabs(v[0]*v[0]+v[1]*v[1]) < 0.0000001 )
            s=1;
          if ( s != s )
          {
            printf ("Fail generateNeighProposals_refit %f \n", s);
//            mexEvalString("drawnow");
          }

          typename EvalEnergyStore<Scalar>::P3 N_scale = n_start * (s*dStart);

          if ( fabs(s*dStart) < 0.0000001 )
            N_scale = n_start;

          norStored->push_back( N_scale ); // same as before
          traStored->push_back( translation );
          rotStored->push_back( rotation );
          result[id] = countNew++;
          valids[id] = 0;
        }
      }
      break;
      }
    case 2: // here simply keep the rotation, xchange the normals
      {
//      EvalEnergyStore<Scalar>::M3 rotation    = (*rotStored)[ currentSolution[segi] ];
      typename EvalEnergyStore<Scalar>::P3 normal = (*norStored)[ currentSolution[segi] ];

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
          norStored->push_back( normal );
          traStored->push_back( (*traStored)[currentSolution[id]] );
          rotStored->push_back( (*rotStored)[currentSolution[id]] );
          result[id] = countNew++;
          valids[id] = 0;
        }
      }
      break;
      }

    default:
    case 0:
      {
      /// distribute current proposal to neighbours: growing
      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;
        if ( currentSolution[segi] != result[id] && valids[id] ) // leave untouched if already proposed here
        {
         result[id] = currentSolution[segi];
         valids[id] = 0;
        }
      }
      break;
      }
    }
  }
  return result;
}


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

  evalSmooth.compute_score_Fuse00( current);
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

  score+= lambda*evalSmooth.compute_score_Fuse00( current);

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
  score = lambda * evalSmooth.compute_score_Fuse00(current);

  return score;
}


// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void evalQPBO2FramesNew( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

#ifndef NDEBUG
  const int nIts = 3;
  int jitterIterations = 3;
#else
  int jitterIterations = 80;
  const int nIts = 80;
#endif

  size_t N,M;
  size_t elements;
  char buffer[256];

  Scalar thresh = 0.1;
  Scalar lambda = 0.02;
  Scalar theta  = 100.0;
  Scalar dispMax = 32.;
  Scalar pixJump = 1.0;
  Scalar rotJump(20.);//   setRotJump (rotJump)
  Scalar depthJump(20.);//   setRotJump (rotJump)
  Scalar occThresh(0.8);
  Scalar motMax = 350; // the maximal motion in the scene

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
  Scalar* K_      = (Scalar*)  mxGetPr(prhs[5]); // K_l
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

  if (nrhs > 32)
    jitterIterations =  (int) (*mxGetPr(prhs[32]));

  //////////////////////////

  Scalar outOfBoundsThresh = occThresh;
  std::srand ( unsigned ( std::time(0) ) );

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
    plhs[0]      =  mxCreateDoubleMatrix(nSegments, 1, mxREAL); //mxCreateDoubleMatrix(N, M, mxREAL); // not double  but ...
    output  = (Scalar*) mxGetPr(plhs[0]);
  }

  ///////////////////////////////////////////////////////
  std::srand ( unsigned ( std::time(0) ) );

  EvalEnergyStore<Scalar> enStore;
  enStore.prepareNormals( normals, nNormals, nStep );
  enStore.setRotTra( rotations, nNormals );

//  printf("generateHoms\n");
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

  std::vector<int> currentSolution( nSegments, 0 );
  std::vector<int> oldSolution( nSegments, 0 );

  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      currentSolution[i]=i;

  aFD.computeVecScores (Segments, nNormals, segImg, currentSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
  std::vector<Scalar> currentData = aFD.getAllScores();
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  // now smoothness cost setup:

  EvalEnergyFullFrame<Scalar> evalSmooth;
  evalSmooth.setRotJump  (rotJump);
  evalSmooth.setDepthJump(depthJump);

  evalSmooth.setGamma ( 1.0 );
  evalSmooth.setRotWeight( theta );
  evalSmooth.set2dMotionMatrix ( K_, Kr_, MC_, mc_, pixJump ); // jump at 0.5 pixel -> needs to set jump smaller

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
  Scalar ENERGY(0);

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

  QPBO<int> q(nSegments, numEdges); // max number of nodes & edges

  ENERGY = getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda);
#ifdef __QUIET__
  printf("------------ Start Energy is %f --------------\n", ENERGY);
  printf("Start Energies (full=data+smooth) %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, currentData ), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda) );

  std::clock_t fullStart(std::clock());
#endif
  // fuse all proposals in a row
  std::vector<int> newFlipped;

  for(int its=0;its <jitterIterations; its++)
  {
    std::vector<int> newPSet = generateNeighProposals_refit(nSegments, edges, currentSolution, newFlipped, enStore, centers, K_);

    newFlipped.clear();

    enStore.flipNormals();
#ifdef __do_depth_only__
    aFD.computeVecScoresDepthOnly (Segments, nNormals, segImg, oobStore, oobOutsideR);// uses warp_noOmp_patchBased
#else
    aFD.computeVecScores (Segments, nNormals, segImg, newPSet, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
#endif
    std::vector<Scalar> newData = aFD.getAllScores();
    enStore.flipNormals();
    ///////////////////////

  Scalar worstDataScore(0);
  Scalar bestDataScore(0);

  for (int i=0; i< nSegments;i++)
  {
    worstDataScore += max( currentData[i], newData [i] );
    bestDataScore  += min( currentData[i], newData [i] );
  }

#ifdef __DEBUG_ENERGIES__
  for ( int i = 0;i < nSegments;i++ )
    oldSolution[i]=currentSolution[i];
#endif

      q.AddNode(nSegments);
      Scalar worstSmoothScore = evalSmooth.compute_score_Fuse( currentSolution, newPSet );

      Scalar scale = (Scalar)(maxValue) / (lambda * worstSmoothScore + worstDataScore);

      Scalar scaleSmo = scale * lambda;
      // add unary data terms:
      for ( int j = 0; j< nSegments; j++ )
      {
        q.AddUnaryTerm(j, scale*(currentData[j]) , 
                          scale*(newData    [j]) );
      }

      for( int j=0; j < idk.size(); j++ )
      {
        if (f00[j]+f11[j] > f10[j]+f01[j])
          non_sub++;
        q.AddPairwiseTerm(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);
      }

      q.Solve();
      q.ComputeWeakPersistencies();

//#define useProbe
#ifndef useProbe
      int allNew = 0;
      for (int j = 0; j < nSegments;j++)
      {
        int x = q.GetLabel(j);
        if (x<0)
          non_sol++;

        if (x==1) {
          allNew++;
          currentSolution[j] = newPSet[j];
          currentData[j]     = newData[j];
          newFlipped.push_back(j);
        };
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
      int allNew = 0;
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
        }
        if (x==1) {
          currentSolution[j] = newPSet[j];
          newFlipped.push_back(j);
          allNew++;
#ifdef _DEBUG
          new1old0[j]=1;
          newListed.push_back(j);
#endif
        };
      }

      int doProbe = 1;int impIts =0;int allImpIts = 0;
      const int maxProbe =3;//3
      if (doProbe && doImprove)
      {
        int nVars = nSegments;
        for (int j = 0; j < nVars;j++)
        {
          int x = q.GetLabel(j);
          if (x>=0)
            q.SetLabel(j,x);
        }

        QPBO<int>::ProbeOptions options;
        int* mapping = (int*) malloc( sizeof(int) * nVars );
        q.Probe( mapping, options );

        int* mapping2 = (int*) malloc( sizeof(int) * nVars );
        for (int i=1;i < maxProbe; i++)
        {
          q.Probe( mapping2, options );
          q.MergeMappings( nVars, mapping, mapping2 );
        }
        free(mapping2);
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
          int pos = noSolutions[j];
          int x = q.GetLabel( mapping[pos]/2 );
          if( (x>=0) && ((x+ mapping[pos]) % 2) )
          {
            currentSolution[pos] = newPSet[pos];
            newFlipped.push_back(pos);
            allImpIts++;
            impIts++;
            allNew++;
          }
        }
        free( mapping );
      }
#endif

      if (0 && allNew > 0)
      {
        enStore.flipNormals();
#ifdef __do_depth_only__
        aFD.computeVecScoresDepthOnly (Segments, nNormals, segImg, oobStore, oobOutsideR);// uses warp_noOmp_patchBased
#else
        aFD.computeVecScores (Segments, nNormals, segImg, currentSolution, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);//, oobOutsideRR, oobOutsideLR);// uses warp_noOmp_patchBased
#endif
        enStore.flipNormals();
        newData = aFD.getAllScores();


        Scalar Last_Energy = ENERGY;
        //      OLD_ENERGY = ENERGY;
        ENERGY = getEnergy ( currentSolution, nSegments, newData, evalSmooth, theta, lambda);

        // reverse all if went wrong
        if (Last_Energy < ENERGY)
        {
          printf("reversing the results ENERGY%f, Last_Energy%f OLD_ENERGY%f \n", ENERGY, Last_Energy, OLD_ENERGY );
          std::copy( oldSolution.begin(), oldSolution.end(), currentSolution.begin() );

          ENERGY = Last_Energy;
          OLD_ENERGY = Last_Energy;
        }
        else
        {
          OLD_ENERGY = Last_Energy;
          std::copy( newData.begin(), newData.end(), currentData.begin() );

          printf("improved the solution  ENERGY%f, Last_Energy%f OLD_ENERGY%f \n", ENERGY, Last_Energy, OLD_ENERGY  );
        }

#ifdef __QUIET__
// expensive to blabber:
        if (OLD_ENERGY != ENERGY)
          printf("It: %d, -", its);
#endif
      }

      getStoreConsistency( nSegments, currentSolution, enStore );
      q.Reset();
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

  enStore.flipNormals();
  printf("Energies (full=data+smooth) comp %f = (%f,%f) 0:N %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, currentData ), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda), getEnergy ( currentSolution, nSegments, dataScores1, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, dataScores1), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda) );

  printf("QPBO took %f\n", double(std::clock()-fullStart)/CLOCKS_PER_SEC);
  printf("ENERGY: %f  NON submodularity detected %d, no solutions: %d\n", ENERGY, non_sub, non_sol);
#endif
  for (int i=0;i < currentSolution.size(); i++)
    output[i] =currentSolution[i];

  enStore.flipNormals();
  const std::vector<EvalEnergyStore<Scalar>::P3>* norStored = enStore.getNormals();
  const std::vector<EvalEnergyStore<Scalar>::P3>* traStored = enStore.getTra( );
  const std::vector<EvalEnergyStore<Scalar>::M3>* rotStored = enStore.getRot( );

  if (nlhs > 1)
  {
    plhs[1]           =  mxCreateDoubleMatrix(3, (*norStored).size(), mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[1]);
    for (int i=0;i < (*norStored).size(); i++)
    {
      EvalEnergyStore<Scalar>::P3 No = (*norStored)[i];
      outputDI[3*i]   = No[0];
      outputDI[3*i+1] = No[1];
      outputDI[3*i+2] = No[2];
    }
  }
  if (nlhs > 2)
  {
    mwSize dims[3] = {4,4,1};
    dims[2] = (*traStored).size();
    plhs[2]           =  mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[2]);
    for (int i=0;i < (*traStored).size(); i++)
    {
      EvalEnergyStore<Scalar>::M3 R = (*rotStored)[i];
      EvalEnergyStore<Scalar>::P3 T = (*traStored)[i];
      outputDI[16*i  ]  = R(0,0);
      outputDI[16*i+1 ] = R(1,0);
      outputDI[16*i+2 ] = R(2,0);
      outputDI[16*i+3 ] = Scalar(0);
      outputDI[16*i+4 ] = R(0,1);
      outputDI[16*i+5 ] = R(1,1);
      outputDI[16*i+6 ] = R(2,1);
      outputDI[16*i+7 ] = Scalar(0);
      outputDI[16*i+8 ] = R(0,2);
      outputDI[16*i+9 ] = R(1,2);
      outputDI[16*i+10] = R(2,2);
      outputDI[16*i+11] = Scalar(0);
      outputDI[16*i+12] = T[0];
      outputDI[16*i+13] = T[1];
      outputDI[16*i+14] = T[2];
      outputDI[16*i+15] = Scalar(1);
    }
  }
}
#undef __do_depth_only__
#undef useProbe
#undef __Small_Improvement__
#endif 
