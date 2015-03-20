/*
 * get_FdF_simple.cpp
 *
 * compute laplace F and dF both in matrix form and in a simple way (no scale2)
 *
 *  Created on: Feb 12, 2010
 *      Author: christop
 */

#include "SparseMatrixComp.h"
#include "mex.h"

//#define USE_INVALIDS

void get_FdF_simple( int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[])
{
int N,M,I;

int elements, nsubs;
double *output1, *output2, *output3;
//mwIndex id, idv;

const double epsilon = 0.0000001;
const double rEpsilon = 1./epsilon;

char buffer[256];

  I     = (int) *mxGetPr(prhs[0]);
  N     = (int) *mxGetPr(prhs[1]);
  M     = (int) *mxGetPr(prhs[2]);

if (nlhs != 1 && nlhs != 2 && nlhs != 3)
  mexErrMsgTxt("Must have 1 or 2 or 3 output arguments");

//Edge-Weights W; Righthand Side B; Vertex Scales S; Init Solution D
//if (nrhs !=U+7)
//  mexErrMsgTxt("Must have amount of unknowns (U)+7 input arguments");

/* prevent you from passing a sparse matrix,
a string matrix, or a complex array. mxIsComplex
is used to determine if there is an imaginary
part of the mxArray. mxIsClass is used to determine
if the mxArray belongs to a particular class */
if ( mxIsComplex(prhs[3])|| mxIsClass(prhs[3],"sparse") || mxIsChar(prhs[3]) )
  mexErrMsgTxt("Input must be real, full, and nonstring");

// Get the number of elements in the input argument
elements=mxGetNumberOfElements(prhs[3]);

sprintf (buffer, "The amount of elements in the matrix is not equal to the amount given by the specified dimensions!\n all: %d product %d   N: %d, M: %d\n", elements, N*M, N, M);
if (elements != N*M)
  mexErrMsgTxt(buffer);

// Get the number of dimensions in array
nsubs=mxGetNumberOfDimensions(prhs[3]);

// Jacobi:
double *weightsH(NULL), *weightsV(NULL), *b1, *b2, *b3, *x1, *x2, *x3;//, *dualW, *scales;
int *validMap(NULL), *invalidMap(NULL);
validMap = (int*) mxGetPr(prhs[3]);

// could be done as well BUT NOT NEEDED CURRENTLY

#ifdef USE_INVALIDS
invalidMap = (int*) mxGetPr (prhs[4]);
#define ADD_TO 1
#else
#ifdef USE_WEIGHTS
weightsH = mxGetPr(prhs[4]); // edge weights
weightsV = mxGetPr(prhs[5]); // edge weights
#  define ADD_TO 2
#else
#  define ADD_TO 0
#endif
#endif


// confusing non-sense: b is not b but always 0
// the b is used as a initial value for x !!!

// x1 is never used as initial solution - just the masked part is filled with x1 

//dualW    = mxGetPr(prhs[8]);
b1        = mxGetPr(prhs[4+ADD_TO]);
x1        = mxGetPr(prhs[5+ADD_TO]);

int rhs = 1;

if (nrhs >= 8+ADD_TO)
  rhs = 2;
if (nrhs >= 10+ADD_TO)
  rhs = 3;

//printf("RHS: %d", nrhs);

b2 = b1;
b3 = b1;
x2 = x1;
x3 = x1;

if (rhs > 1)
{
  b2        = mxGetPr(prhs[6+ADD_TO]);
  x2        = mxGetPr(prhs[7+ADD_TO]);
}
if (rhs > 2)
{
  b3        = mxGetPr(prhs[8+ADD_TO]);
  x3        = mxGetPr(prhs[9+ADD_TO]);
}

// Create the output array(s)
plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
output1 = mxGetPr(plhs[0]);
memcpy(output1, x1, sizeof(double)*N*M);

if (nlhs >= 2)
{
  plhs[1] = mxCreateDoubleMatrix(N, M, mxREAL);
  output2 = mxGetPr(plhs[1]);
  memcpy(output2, x2, sizeof(double)*N*M);
}
if (nlhs >= 3)
{
  plhs[2] = mxCreateDoubleMatrix(N, M, mxREAL);
  output3 = mxGetPr(plhs[2]);
  memcpy(output3, x3, sizeof(double)*N*M);
}

//printf("N %d, M %d -- nlhs %d, rhs %d \n",N, M, nlhs, rhs);

if (nlhs == 1)
{
  SolveLaplace<1> slp( N, M, b1, validMap, output1);
  if (weightsH != 0 )
    slp.setWeights( weightsH, weightsV );
  if (invalidMap != NULL)
    slp.setInvalid(invalidMap);

  slp.setup_Matrix( I );
}
if (nlhs == 2)
{
  SolveLaplace<2> slp( N, M, b1, b2, validMap, output1, output2);
  if (weightsH != 0 )
    slp.setWeights( weightsH, weightsV );
  if (invalidMap != NULL)
    slp.setInvalid(invalidMap);
  slp.setup_Matrix( I );
}
if (nlhs == 3)
{
  SolveLaplace<3> slp( N, M, b1, b2, b3, validMap, output1, output2, output3);
  if (weightsH != 0 )
    slp.setWeights( weightsH, weightsV );
  if (invalidMap != NULL)
    slp.setInvalid(invalidMap);
  slp.setup_Matrix( I );
}
}

#undef USE_INVALIDS
