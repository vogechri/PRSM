/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */

//#define _NO_OPENMP

//#include <string.h>
#include <iostream>
#include "mex.h"

#include "Interpolate.cpp"
#include "Interpolate_2ndDerivatives.cpp"
#include "nrutil.h"
#include "nrutil.c"

//#include "bcuint.cpp"
//#include "bcucof.cpp"

#define zero 1
//#define DebugOut
//#define DebugOut2


// no openmp but gcc 4.1
// mex CC=g++-4.1 CXX=g++-4.1 CFLAGS='-fPIC' LD=g++-4.1 -O -I/home/christop/CPP -lm -output Smooth_mex MexInterface.cpp

// mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp
// mex CC=g++-4.2 CXX=g++-4.2 CFLAGS='-fPIC -fopenmp' LD=g++-4.2 -O -I/home/christop/CPP -lm -lgomp -output Smooth_mex MexInterface.cpp
void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{

//  int m, n, i, e, j, it;
//  int N,M,I;

//  int elements, nsubs;
//  double *output;
//  double *temp, *weights;
//  mwIndex id, idv;

  char buffer[256];

  /* check: only one input and one output argument */
  if (nrhs < 3)
    mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM N,M; Image I; X-Indices X and Y-Indices Y\n)");

//  if (!mxIsDouble(prhs[3]))
//  {
//    mexErrMsgTxt("First three input argument must be a double.");
//  }

  bool isDouble(1);
  bool isSingle(1);
  bool isInt32(1);
  bool isInt16(1);

  bool isSingleLast2(0);
  bool isDoubleLast2(0);

  for(int i = 0;i<3;i++)
    if (!mxIsDouble(prhs[i]))
      isDouble = 0;

  for(int i = 0;i<3;i++)
    if (!mxIsSingle(prhs[i]))
      isSingle = 0;

  if (!mxIsInt32(prhs[0]))
    isInt32 = 0;
  else
  {
    isSingleLast2=true;
    isDoubleLast2=true;
    for( int i = 1; i < 3; i++ )
    {
      if ( !mxIsSingle(prhs[i]) )
        isSingleLast2 = 0;
      if ( !mxIsDouble(prhs[i]) )
        isDoubleLast2 = 0;
    }
  }

  if (!mxIsInt16(prhs[0]))
    isInt16 = 0;
  else
  {
    isSingleLast2=true;
    isDoubleLast2=true;
    for( int i = 1; i < 3; i++ )
    {
      if ( !mxIsSingle(prhs[i]) )
        isSingleLast2 = 0;
      if ( !mxIsDouble(prhs[i]) )
        isDoubleLast2 = 0;
    }
  }


  bool isInt(isDoubleLast2 | isSingleLast2);

  if (!isDouble && ! isSingle && !isInt)
    mexErrMsgTxt("First three input arguments must be all double or all single or first int32, other single!");

//  if (nlhs > 4)
//   Interpolate_2ndDerivatives ( nlhs, plhs, nrhs, prhs );
//  else
  {
    // normal algorithm
    if (isDouble)
      Interpolate ( nlhs, plhs, nrhs, prhs );
    if ( isSingle)
      Interpolate_single ( nlhs, plhs, nrhs, prhs );
//    if ( isInt16 )
//      Interpolate_int16 ( nlhs, plhs, nrhs, prhs );
//    if (isSingleLast2 & isInt32)
//      Interpolate_int32 ( nlhs, plhs, nrhs, prhs );
//    if (isDoubleLast2 & isInt32)
//      Interpolate_int32_d2 ( nlhs, plhs, nrhs, prhs );

  }  
  //  Interpolate_2ndDerivatives( nlhs, plhs, nrhs, prhs );

}


