/*
 * MexInterface.cpp
 *
 *  Created on: Feb 12, 2010
 *      Author: christoph vogel
 */

#include <iostream>  // iostream must be before mex -> else redefinition
#include "mex.h"

#include "FullStepMatlab.h"

//mex mex_sum_openmp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
// no openmp but gcc 4.1
// mex CC=g++-4.1 CXX=g++-4.1 CFLAGS='-fPIC' LD=g++-4.1 -O -I/home/christop/CPP -lm -output Smooth_mex MexInterface.cpp

// mex CXX=g++-4.2 CXX=g++-4.2 LD=g++-4.2 -I/home/christop/CPP -lm -output mexTest IRLS_mex.cpp
// mex CC=g++-4.2 CXX=g++-4.2 CFLAGS='-fPIC -fopenmp' LD=g++-4.2 -O -I/home/christop/CPP -lm -lgomp -output Smooth_mex MexInterface.cpp
void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  growGrid ( nlhs, plhs, nrhs, prhs );
}
