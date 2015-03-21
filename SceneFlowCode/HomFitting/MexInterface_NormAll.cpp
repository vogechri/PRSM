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

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{

#ifdef __doStereo__
  if (!mxIsDouble(prhs[5]))
  {
    mexErrMsgTxt("Fifth input argument must be a double.");
  }

    if (mxIsDouble(prhs[3]))
          mexErrMsgTxt("First four input argument must be a uint-8.");
#endif

#ifdef __doStereo__
    get_Stereo ( nlhs, plhs, nrhs, prhs );
#else
    get_L1HomoNAllSeg( nlhs, plhs, nrhs, prhs );
#endif
}


