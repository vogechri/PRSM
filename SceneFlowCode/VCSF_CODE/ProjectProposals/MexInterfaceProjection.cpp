/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */

#include <iostream>  // iostream must be before mex -> else redefinition
#include "mex.h"

#ifdef _projectMVPs_
#define _returnProjectionOnly_
#define _useNoCenters_
#include "linearMotionProposals.h"
#endif

#ifdef _projectMVPs_grow_
#define _useNoCenters_
//#define _returnProjectionOnly_
#include "linearMotionProposals.h"
#endif

#ifdef _projectMVP_labels_
#include "ProjectSegmentation.h" // with growing !
#endif

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{

  // project and grow:
#ifdef _projectMVP_labels_
  ProjectSegmentation( nlhs, plhs, nrhs, prhs );
#endif

#ifdef _projectMVPs_grow_
//  ProjectSegmentation( nlhs, plhs, nrhs, prhs );
  generateLinearMotionProposals( nlhs, plhs, nrhs, prhs );
#endif

  // pure projection:
#ifdef _projectMVPs_
  generateLinearMotionProposals( nlhs, plhs, nrhs, prhs );
#endif
}
