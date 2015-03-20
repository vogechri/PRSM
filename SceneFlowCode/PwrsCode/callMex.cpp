/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */


#define _NO_OPENMP
#include <iostream>
#include "mex.h"
#include <mat.h>
//#include "engine.h"

#define __QUIET__

//#define occlusionVersion
#define __perSegmentVersion__

#ifndef occlusionVersion

//#include "evalQPBO2FrameJitterMix.h"         // without occlusion handling
#include "evalQPBO2FrameCompressedOcc.h"       // without occlusion handling

//#include "evalQPBO2FrameNew.h"         // without occlusion handling
#else
//#include "evalQPBO2FrameNewMultiOcc.h"  // with occlusion handling
#include "SuperPixelSeg3.h"  // with occlusion handling
#endif

#define zero 1
//#define DebugOut
//#define DebugOut2

//mex mex_sum_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
// no openmp but gcc 4.1
// mex CC=g++-4.1 CXX=g++-4.1 CFLAGS='-fPIC' LD=g++-4.1 -O -I/home/christop/CPP -lm -output Smooth_mex MexInterface.cpp

// mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp
// mex CC=g++-4.2 CXX=g++-4.2 CFLAGS='-fPIC -fopenmp' LD=g++-4.2 -O -I/home/christop/CPP -lm -lgomp -output Smooth_mex MexInterface.cpp

typedef void (*mexFunction_t) (int nargout, mxArray *pargout[], int nargin, const mxArray *pargin[]);

int main ( int argc, const char *argv[])
{

//  int m, n, i, e, j, it;
//  int N,M,I;

//  int elements, nsubs;
//  double *output;
//  double *temp, *weights;

//  char buffer[256];

  /* check: only one input and one output argument */
//  if (nrhs < 10)
//    mexErrMsgTxt("Must have 9 input arguments:\n iterations I; size of data matrix NxM N,M; double: theta, Edge-Weights Wu; Edge-Weights Wd; Righthand Side B; Vertex Scales S; Scales S2; Init Solution D)");

//  if (nrhs < 13)
 //   mexErrMsgTxt("Must have 12 input arguments:\n iterations I; size of data matrix NxM N,M; double: theta, Edge-Weights Wu; Edge-Weights Wd; Righthand Side B; Vertex Scales S; Scales S2; Init Solution D)");


//  Engine *ep(NULL);
  char buf[1024];
  /*
  if (!(ep = engOpen("matlab - nodisplay")))
  {
    fprintf(stderr, "Can not open matlab engine\n");
    return -1;
  }
  engOutputBuffer(ep, buf, 1023);
  */


  MATFile * pmat;
//  pmat = matOpen("testSmooth.mat","r");
  
//  engEvalString(ep, "load input.mat");
//  mxArray *arg1a = matGetVariable( pmat, "d");
//  mxArray *arg2a = matGetVariable( pmat, "N");
//  mxArray *arg3a = matGetVariable( pmat, "M");

//  mxArray *arg1 = engGetVariable(ep, "d");
//  mxArray *arg2 = engGetVariable(ep, "N");
//  mxArray *arg3 = engGetVariable(ep, "M");

#ifdef __perSegmentVersion__

#ifndef occlusionVersion
  pmat = matOpen("BugSeg_117_10.mat","r");
  
//  engEvalString(ep, "load input.mat");
  mxArray *arg1 = matGetVariable( pmat, "I1");
  mxArray *arg2 = matGetVariable( pmat, "I2");
  mxArray *arg3 = matGetVariable( pmat, "I3");
  mxArray *arg4 = matGetVariable( pmat, "I4");
  mxArray *arg5 = matGetVariable( pmat, "pp2d");
  mxArray *arg6 = matGetVariable( pmat, "Kr");
  mxArray *arg7 = matGetVariable( pmat, "segIds");
  mxArray *arg8 = matGetVariable( pmat, "nprop");
  mxArray *arg9 = matGetVariable( pmat, "dt");

  mxArray *arg10 = matGetVariable( pmat, "ds");
  mxArray *arg11 = matGetVariable( pmat, "ts");
  mxArray *arg12 = matGetVariable( pmat, "segEdges");
  mxArray *arg13 = matGetVariable( pmat, "segWeights");
  mxArray *arg14 = matGetVariable( pmat, "Tr");
  mxArray *arg15 = matGetVariable( pmat, "rtprop");
  mxArray *arg16 = matGetVariable( pmat, "centers2D");
  mxArray *arg17 = matGetVariable( pmat, "Rr");
  mxArray *arg18 = matGetVariable( pmat, "tj");
  mxArray *arg19 = matGetVariable( pmat, "dj");

  mxArray *arg20 = matGetVariable( pmat, "dD");
  mxArray *arg21 = matGetVariable( pmat, "oob");
  mxArray *arg22 = matGetVariable( pmat, "segimg");
  mxArray *arg23 = matGetVariable( pmat, "ioobs");
  mxArray *arg24 = matGetVariable( pmat, "ioobsFlow");
  mxArray *arg25 = matGetVariable( pmat, "ioobsFlow2");
  mxArray *arg26 = matGetVariable( pmat, "maxMot");
  mxArray *arg27 = matGetVariable( pmat, "wx");
  mxArray *arg28 = matGetVariable( pmat, "wy");
  mxArray *arg29 = matGetVariable( pmat, "wxy");
  mxArray *arg30 = matGetVariable( pmat, "wixy");
  mxArray *arg31 = matGetVariable( pmat, "segImgFlip");

//  'I1','I2','I3','I4', 'pp2d', 'Kr', 'segIds', 'nprop', 'dt', 'ds', 'tS',
// 'segEdges', 'segWeights', 'Tr', 'rtprop', 'centers2D', 'Rr', 'tj', 'dj',
// 'dD', 'oob', 'segimg', 'ioobs', 'ioobsFlow', 'ioobsFlow2', 'maxMot', 'wx', 'wy', 'wxy', 'wixz', 'segImgFlip');



  mxArray *pargout[2] = {0};//,0};
  const mxArray *pargin[31] = { arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9, arg10, 
                               arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, 
                               arg21, arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30, arg31};

  int nlhs = 1;
  int nrhs = 31;
#else
  pmat = matOpen("image123OccFull.mat","r");
  
//  engEvalString(ep, "load input.mat");
  mxArray *arg1 = matGetVariable( pmat, "I1");
  mxArray *arg2 = matGetVariable( pmat, "I2");
  mxArray *arg3 = matGetVariable( pmat, "I3");
  mxArray *arg4 = matGetVariable( pmat, "I4");
  mxArray *arg5 = matGetVariable( pmat, "pp2d");
  mxArray *arg6 = matGetVariable( pmat, "Kr");
  mxArray *arg7 = matGetVariable( pmat, "segIds");
  mxArray *arg8 = matGetVariable( pmat, "nprop");
  mxArray *arg9 = matGetVariable( pmat, "dt");

  mxArray *arg10 = matGetVariable( pmat, "ds");
  mxArray *arg11 = matGetVariable( pmat, "tS");
  mxArray *arg12 = matGetVariable( pmat, "segEdges");
  mxArray *arg13 = matGetVariable( pmat, "segWeights");
  mxArray *arg14 = matGetVariable( pmat, "Tr");
  mxArray *arg15 = matGetVariable( pmat, "rtprop");
  mxArray *arg16 = matGetVariable( pmat, "centers2D");

  mxArray *arg17 = matGetVariable( pmat, "Rr");
  mxArray *arg18 = matGetVariable( pmat, "tj");
  mxArray *arg19 = matGetVariable( pmat, "dj");

  mxArray *arg20 = matGetVariable( pmat, "dD");
  mxArray *arg21 = matGetVariable( pmat, "oob");
  mxArray *arg22 = matGetVariable( pmat, "occ");

  mxArray *arg23 = matGetVariable( pmat, "segimg");

  mxArray *arg24 = matGetVariable( pmat, "occList");
  mxArray *arg25 = matGetVariable( pmat, "occListL");
  mxArray *arg26 = matGetVariable( pmat, "occListR");

  mxArray *arg27 = matGetVariable( pmat, "Kl");
  mxArray *arg28 = matGetVariable( pmat, "segImgFlip");

  mxArray *arg29 = matGetVariable( pmat, "ioobs");
  mxArray *arg30 = matGetVariable( pmat, "ioobsFlow");
  mxArray *arg31 = matGetVariable( pmat, "ioobsFlow");
  mxArray *arg32 = matGetVariable( pmat, "maxMot");

  mxArray *arg33 = matGetVariable( pmat, "wx");
  mxArray *arg34 = matGetVariable( pmat, "wy");
  mxArray *arg35 = matGetVariable( pmat, "wxy");
  mxArray *arg36 = matGetVariable( pmat, "wixy");

  mxArray *pargout[2] = {0};//,0};
  const mxArray *pargin[36] = { arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9, arg10, 
                               arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, 
                               arg21, arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30,
                               arg31, arg32, arg33, arg34, arg35, arg36 };

  int nlhs = 1;
  int nrhs = 36;

#endif

#else
  pmat = matOpen("image163Pix.mat","r");
  
//  engEvalString(ep, "load input.mat");
  mxArray *arg1 = matGetVariable( pmat, "I1");
  mxArray *arg2 = matGetVariable( pmat, "I2");
  mxArray *arg3 = matGetVariable( pmat, "I3");
  mxArray *arg4 = matGetVariable( pmat, "I4");
  mxArray *arg5 = matGetVariable( pmat, "pp2d");
  mxArray *arg6 = matGetVariable( pmat, "Kr");
  mxArray *arg7 = matGetVariable( pmat, "segIds");
  mxArray *arg8 = matGetVariable( pmat, "nprop");
  mxArray *arg9 = matGetVariable( pmat, "dt");

  mxArray *arg10 = matGetVariable( pmat, "ds");
  mxArray *arg11 = matGetVariable( pmat, "tS");
  mxArray *arg12 = matGetVariable( pmat, "segEdges");
  mxArray *arg13 = matGetVariable( pmat, "segWeights");
  mxArray *arg14 = matGetVariable( pmat, "Tr");
  mxArray *arg15 = matGetVariable( pmat, "rtprop");
  mxArray *arg16 = matGetVariable( pmat, "centers2D");
  mxArray *arg17 = matGetVariable( pmat, "Rr");

  mxArray *arg18 = matGetVariable( pmat, "homos");
  mxArray *arg19 = matGetVariable( pmat, "viewNormals");

  mxArray *arg20 = matGetVariable( pmat, "tj");
  mxArray *arg21 = matGetVariable( pmat, "dj");

  mxArray *arg22 = matGetVariable( pmat, "dD");
  mxArray *arg23 = matGetVariable( pmat, "oob");
  mxArray *arg24 = matGetVariable( pmat, "segImgFlip");

  mxArray *arg25 = matGetVariable( pmat, "wx");
  mxArray *arg26 = matGetVariable( pmat, "wy");
  mxArray *arg27 = matGetVariable( pmat, "wxy");
  mxArray *arg28 = matGetVariable( pmat, "wixy");

  mxArray *arg29 = matGetVariable( pmat, "Kl");
  mxArray *arg30 = matGetVariable( pmat, "ps");
  mxArray *arg31 = matGetVariable( pmat, "maxMot");
  mxArray *arg32 = matGetVariable( pmat, "its");

  mxArray *arg33 = matGetVariable( pmat, "occW");

  mxArray *arg34 = matGetVariable( pmat, "ioobs");
  mxArray *arg35 = matGetVariable( pmat, "ioobsFlow");
  mxArray *arg36 = matGetVariable( pmat, "ioobsFlow");

  mxArray *pargout[2] = {0};//,0};
  const mxArray *pargin[36] = { arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9, arg10, 
                               arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, 
                               arg21, arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30,
                               arg31, arg32, arg33,arg34,arg35, arg36};
  int nlhs = 1;
  int nrhs = 36;
  
//save( 'segImgFlip', 'wx', 'wy', 'wxy', 'wixy', 'Kl', 'ps', 'maxMot', 'occW' ,'ioobs', 'ioobsFlow', 'ioobsFlow2');
#endif
  //if (!mxIsDouble(prhs[3]))
  //{
  //  mexErrMsgTxt("First three input argument must be a double.");
  //}

  // new per triangle with self weighting: best!
  evalQPBO2FramesNew ( nlhs, pargout, nrhs, pargin );
//  SuperPix2Frames( nlhs, pargout, nrhs, pargin );
  int notLast = 1;
//  engClose(ep);
}


