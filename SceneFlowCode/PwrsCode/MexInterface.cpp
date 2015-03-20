/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2012
 *      Author: christoph vogel
 */

#include <iostream>  // iostream must be before mex -> else redefinition
#include "mex.h"

#ifdef _OCC_
#define _occHandlingOn_
#ifdef _SEG_
// gather statistics:
//#include "evalQPBO2FrameSegOccCounting.h"    // without occlusion handling
#include "evalQPBO2FrameSegOcc.h"    // without occlusion handling
#elif defined (_JIT_)
#include "evalQPBO2FrameJitterMixOcc.h"     // without occlusion handling
#elif defined (_FUSE_)
#include "evalQPBO2FrameFuseOcc.h"          // without occlusion handling
#elif defined (_PIX_)
#include "SuperPixelSegOcc.h"
#endif

#else

#ifdef _SEG_
#include "evalQPBO2FrameSeg.h"    // without occlusion handling
#elif defined (_JIT_)
#include "evalQPBO2FrameJitterMix.h"     // without occlusion handling
#elif defined (_FUSE_)
#include "evalQPBO2FrameFuse.h"          // without occlusion handling
#elif defined (_PIX_)
#include "SuperPixelSeg.h"
#endif

#endif

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
#if defined (_PIX_)
  SuperPix2Frames( nlhs, plhs, nrhs, prhs );
#else
  evalQPBO2FramesNew( nlhs, plhs, nrhs, prhs );
#endif
}
