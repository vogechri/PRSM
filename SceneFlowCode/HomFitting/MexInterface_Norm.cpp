/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */


//#define _NO_OPENMP

#include <string.h>
#include <iostream>
#include <vector>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h"

//#define __doStereo__

#ifdef __doStereo__
#include "get_Stereo.cpp"
#else
#include "getDerivativeHomoRotL1.cpp"
#include "getDerivativeHomoNormL1.cpp"
#include "getDerivativeHomoNormL1_ALLSeg.cpp"
//#include "getDerivativeHomoRotL1_ALLSeg.cpp"
#include "getDerivativeHomoRotL1_AllSeg.cpp"
//#include "getDerivativeHomoRotL1_ALLSeg_rightonly.cpp"
#endif

#define zero 1
//#define DebugOut
//#define DebugOut2

//mex mex_sum_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
// no openmp but gcc 4.1
// mex CC=g++-4.1 CXX=g++-4.1 CFLAGS='-fPIC' LD=g++-4.1 -O -I/home/christop/CPP -lm -output Smooth_mex MexInterface.cpp

// mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp
// mex CC=g++-4.2 CXX=g++-4.2 CFLAGS='-fPIC -fopenmp' LD=g++-4.2 -O -I/home/christop/CPP -lm -lgomp -output Smooth_mex MexInterface.cpp
void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{

  int m, j;
  double *weights;

  // check: only one input and one output argument
//  if (nrhs < 10)
//    mexErrMsgTxt("Must have 9 input arguments:\n iterations I; size of data matrix NxM N,M; double: theta, Edge-Weights Wu; Edge-Weights Wd; Righthand Side B; Vertex Scales S; Scales S2; Init Solution D)");

//  if (nrhs < 13)
 //   mexErrMsgTxt("Must have 12 input arguments:\n iterations I; size of data matrix NxM N,M; double: theta, Edge-Weights Wu; Edge-Weights Wd; Righthand Side B; Vertex Scales S; Scales S2; Init Solution D)");


#ifdef __doStereo__
  if (!mxIsDouble(prhs[5]))
  {
    mexErrMsgTxt("Fifth input argument must be a double.");
  }

    if (mxIsDouble(prhs[3]))
          mexErrMsgTxt("First four input argument must be a uint-8.");
#endif

  // WORKING
  // only return the derivative, now with left-right and up-down weghting, because of boundary
  //TVL1Minimize_lrud_mex
//  get_FdF    ( nlhs, plhs, nrhs, prhs );

  // LATEST lrud_Mix_Mex Mixed robust functions here
  //get_FdF_Mix ( nlhs, plhs, nrhs, prhs );


#ifdef __doStereo__
    get_Stereo ( nlhs, plhs, nrhs, prhs );
#else
//    get_L1HomoNAllSeg( nlhs, plhs, nrhs, prhs );
    get_L1HomoN ( nlhs, plhs, nrhs, prhs );

//    get_L1HomoAllSeg ( nlhs, plhs, nrhs, prhs );
//    get_L1Homo ( nlhs, plhs, nrhs, prhs );


//    get_L1HomoAllSeg_right( nlhs, plhs, nrhs, prhs ); / not working good
#endif

  // return the derivative and the preconditioner pi: A*pi*(pi^-1*x) - (pi^-1*b) => A*pi*x^ - b^ to solve. To receive x x = pi*x^
  // this implies b=b(x) and therefore b^ can be computed as b(pi^-1*x). OR
  // pi^-1*x = (A*pi)^-1 * (pi^-1* b) = pi^-1*A^-1 * (pi^-1* b)
  //
  //
  //
// NOT WORKING
  //  get_FdFPCG ( nlhs, plhs, nrhs, prhs );


  // WORKING: full procedure with reparameterization (hope was speedup)both entities in depth
//  get_FdF_Sub ( nlhs, plhs, nrhs, prhs );

  // WORKING: full procedure without reparameterization, both entities in disparity
//  get_FdF_NoSub ( nlhs, plhs, nrhs, prhs );

  // WORKING: new tryout, b already in disparity, one in disparity(b) one in depth(x)
// used in smoothing depth while fixing the disparity
  //get_FdF_newSub ( nlhs, plhs, nrhs, prhs );
}


