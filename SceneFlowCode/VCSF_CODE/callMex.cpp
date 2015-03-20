/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */

#ifdef _DEBUG
#define _NO_OPENMP
#endif

#include <iostream>
#include "mex.h"
#include <mat.h>

#define __perSegmentVersion__
//#define ProjectSeg

#define _testLocalReplacement_

#ifndef ProjectSeg

#ifdef __perSegmentVersion__
#include "VCSF_GC.h"
#else
#include "VCSF_Pix_GC.h"
#endif

#else
#include "ProjectSegmentation.h" // with growing
#endif

typedef void (*mexFunction_t) (int nargout, mxArray *pargout[], int nargin, const mxArray *pargin[]);

int main ( int argc, const char *argv[])
{

  MATFile * pmat;

#ifdef ProjectSeg

  const char **dir;
  const char *name;
  int	  ndir;
  mxArray *pa;

//  pmat = matOpen("projectSeg127_10.mat","r");
  pmat = matOpen("projectSeg2nd_127_10.mat","r");
  dir = (const char **)matGetDir(pmat, &ndir);

  mxFree(dir);

  /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
  if (matClose(pmat) != 0) {
    return(1);
  }

//  pmat = matOpen("projectSeg127_10.mat","r");
  pmat = matOpen("projectSeg2nd_127_10.mat","r");

  mxArray *pargout[1] = {0};//,0};
  const mxArray *pargin[20];
  for (int i=0; i < ndir; i++) 
  {
    pargin[i] = matGetNextVariable(pmat, &name);
    if (pargin[i] == NULL) 
    {
	    return(1);
    }
  }

  int nlhs = 1;
  int nrhs = ndir;

#else

#ifdef __perSegmentVersion__

#ifndef occlusionVersion
//  pmat = matOpen("perSeg014_10.mat","r");
  pmat = matOpen("LargeperSeg014_10.mat","r");
  pmat = matOpen("per3FSeg138_10.mat","r");

  pmat = matOpen("perSeg127_10.mat","r");
  pmat = matOpen("perSeg131_10.mat","r");

  pmat = matOpen("perSegAuto141_10.mat","r");
  pmat = matOpen("perSegAuto191_10.mat","r");

  pmat = matOpen("perSegNew047_10.mat","r");//2067443
//  pmat = matOpen("perSegAuto047_10.mat","r");

  pmat = matOpen("perSegNew3f_144_10.mat","r");
//  pmat = matOpen("perSegFine3f_116_10.mat","r");

  pmat = matOpen("simple2F_002_10.mat","r");

  pmat = matOpen("simpleBug2F_002_10.mat","r");

  pmat = matOpen("3frameFail_002_10.mat","r");

#ifdef _testLocalReplacement_
  pmat = matOpen("replace_121_10.mat","r");
#endif

  const char **dir;
  const char *name;
  int	  ndir;
//  mxArray *pa;
  
  dir = (const char **)matGetDir(pmat, &ndir);

  mxFree(dir);

  /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
  if (matClose(pmat) != 0) {
    return(1);
  }

//  pmat = matOpen("perSeg014_10.mat", "r");
//  pmat = matOpen("LargeperSeg014_10.mat","r");
  pmat = matOpen("per3FSeg138_10.mat","r");

  pmat = matOpen("perSeg127_10.mat","r");
  pmat = matOpen("perSeg131_10.mat","r");
  pmat = matOpen("perSegAuto141_10.mat","r");
  pmat = matOpen("perSegAuto191_10.mat","r");
  pmat = matOpen("perSegNew047_10.mat","r");
//  pmat = matOpen("perSegAuto047_10.mat","r");
  pmat = matOpen("perSegNew3f_144_10.mat","r");
//  pmat = matOpen("perSegFine3f_116_10.mat","r");

  pmat = matOpen("simple2F_002_10.mat","r");

  pmat = matOpen("simpleBug2F_002_10.mat","r");

  pmat = matOpen("3frameFail_002_10.mat","r");

#ifdef _testLocalReplacement_
  pmat = matOpen("replace_121_10.mat","r");
#endif

  if (pmat == NULL) {
    return(1);
  }

//  const mxArray *argin[38];
  const mxArray *pargin[43];

  for (int i=0; i < ndir; i++) 
  {
    mxArray* temp = matGetNextVariable(pmat, &name);
    int number = atoi( &name[1] );
    pargin[number-1] = temp;
    if (pargin[number-1] == NULL) 
	    return(1);
  }


//  engEvalString(ep, "load input.mat");

  mxArray *pargout[2] = {0};//,0};
//  const mxArray *pargin[31] = { arg1,  arg2,  arg3,  arg4,  arg5,  arg6,  arg7,  arg8,  arg9, arg10, 
//                               arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20, 
//                               arg21, arg22, arg23, arg24, arg25, arg26, arg27, arg28, arg29, arg30, arg31};

  int nlhs = 1;
  int nrhs = ndir;
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
//  pmat = matOpen("image163Pix.mat","r");
  const char **dir;
  const char *name;
  int	  ndir;
  mxArray *pa;

  pmat = matOpen("per3FPix158_10.mat","r");
  pmat = matOpen("perPixNew010_10.mat","r");
//  pmat = matOpen("perPixNew3f_144_10.mat","r");

  dir = (const char **)matGetDir(pmat, &ndir);

  mxFree(dir);

  /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
  if (matClose(pmat) != 0) {
    return(1);
  }

  pmat = matOpen("per3FPix158_10.mat","r");
  pmat = matOpen("perPixNew010_10.mat","r");
//  pmat = matOpen("perPixNew3f_144_10.mat","r");

  mxArray *pargout[1] = {0};//,0};
  const mxArray *pargin[57];
  for (int i=0; i < ndir; i++) 
  {
    pargin[i] = matGetNextVariable(pmat, &name);
    if (pargin[i] == NULL) 
    {
	    return(1);
    }
  }

  int nlhs = 3;
  int nrhs = ndir;



//save( 'segImgFlip', 'wx', 'wy', 'wxy', 'wixy', 'Kl', 'ps', 'maxMot', 'occW' ,'ioobs', 'ioobsFlow', 'ioobsFlow2');
#endif

#endif

#ifdef ProjectSeg

  ProjectSegmentation( nlhs, pargout, nrhs, pargin );
//   generateLinearMotionProposals( nlhs, plhs, nrhs, prhs );

#else
  run_VCSF ( nlhs, pargout, nrhs, pargin );

#endif
  printf("end");
  //  engClose(ep);
}