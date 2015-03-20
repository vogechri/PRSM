/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */

#ifdef _DEBUG
//#define _NO_OPENMP
#endif

//#define _QUIET_

#include <iostream>  // iostream must be before mex -> else redefinition
#include "mex.h"

// either one of these here 
//#define _SEG_
//#define _PIX_

//#define _testLocalReplacement_

#ifdef _PIX_
#include "VCSF_Pix_GC.h"
#else
#include "VCSF_GC.h"
#endif

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  run_VCSF( nlhs, plhs, nrhs, prhs );
}
