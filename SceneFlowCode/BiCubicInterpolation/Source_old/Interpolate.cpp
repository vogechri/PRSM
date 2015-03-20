/*
 * Interpolate.cpp
 *
 *  Created on: Mar 1, 2010
 *      Author: christop
 */

#include <math.h>
#include "bcuint.cpp"

#include "BicubicIp.h"
#include "BilinearIp.h"

#define NRANSI
#include "nrutil.h"
#undef NRANSI

//#define DebugOut

#define use_class

/*
void Interpolate_int32_d2( int nlhs, mxArray *plhs[],
                           int nrhs, const mxArray *prhs[])
{
  int N(0),M(0);

  int elementsX, elementsY, *I, *outputF;
  double *X, *Y;

  char buffer[256];

  if (nrhs < 3)
    mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM Image I; X-Indices X and Y-Indices Y\n)");

  N = (int) mxGetN(prhs[0]);
  M = (int) mxGetM(prhs[0]);

  I = (int*) mxGetPr(prhs[0]);
  X = (double*) mxGetPr(prhs[1]);
  Y = (double*) mxGetPr(prhs[2]);

  if ((nlhs !=3) && (nlhs !=1))
    mexErrMsgTxt("Must have 1 or 3 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY] )");

  elementsX=(int)mxGetNumberOfElements(prhs[1]);
  elementsY=(int)mxGetNumberOfElements(prhs[2]);

  mwSize dims[2];
  dims[0] = elementsX;
  dims[1] = 1;

  sprintf (buffer, "There is some Error %d elementsX, %d elementsY.\n", elementsX, elementsY);

  if (elementsX != elementsY)
    mexErrMsgTxt(buffer);

  //printf("Elements: img(%d,%d) val: %d \n", N,M, elementsX);
//  plhs[0]  = mxCreateNumericArray(2, dims, mxINT16_CLASS, mxREAL); // the function value
  plhs[0]  = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL); // the function value
  outputF  = (int*) mxGetPr(plhs[0]);

  //BiCubicIp<float> bcIP(N, M, I);
  BiLinearIPInt<double, int> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  //bcIP.interpolate_intOmp( elementsX, X, Y, outputF );
  bcIP.interpolate_int( elementsX, X, Y, outputF );
}

void Interpolate_int16( int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{
  int N(0),M(0);

  int elementsX, elementsY;
  short *I, *outputF;
  double *X, *Y;

  char buffer[256];

  if (nrhs < 3)
    mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM Image I; X-Indices X and Y-Indices Y\n)");

  N = (int) mxGetN(prhs[0]);
  M = (int) mxGetM(prhs[0]);

  I = (short*) mxGetPr(prhs[0]);
  X = (double*) mxGetPr(prhs[1]);
  Y = (double*) mxGetPr(prhs[2]);

  if ((nlhs !=3) && (nlhs !=1))
    mexErrMsgTxt("Must have 1 or 3 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY] )");

  elementsX=(int)mxGetNumberOfElements(prhs[1]);
  elementsY=(int)mxGetNumberOfElements(prhs[2]);

  mwSize dims[2];
  dims[0] = elementsX;
  dims[1] = 1;

  sprintf (buffer, "There is some Error %d elementsX, %d elementsY.\n", elementsX, elementsY);

  if (elementsX != elementsY)
    mexErrMsgTxt(buffer);

  //printf("Elements: img(%d,%d) val: %d \n", N,M, elementsX);
//  plhs[0]  = mxCreateNumericArray(2, dims, mxINT16_CLASS, mxREAL); // the function value
  plhs[0]  = mxCreateNumericArray(2, dims, mxINT16_CLASS, mxREAL); // the function value
  outputF  = (short*) mxGetPr(plhs[0]);

  //BiCubicIp<float> bcIP(N, M, I);
  BiLinearIPInt<double, short> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  //bcIP.interpolate_intOmp( elementsX, X, Y, outputF );
  bcIP.interpolate_int16( elementsX, X, Y, outputF );
}

  
  void Interpolate_int32 ( int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[])
{
  int N(0),M(0);

  int elementsX, elementsY, *I, *outputF;
  float *X, *Y;

  char buffer[256];

  if (nrhs < 3)
    mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM Image I; X-Indices X and Y-Indices Y\n)");

  N = (int) mxGetN(prhs[0]);
  M = (int) mxGetM(prhs[0]);

  I = (int*) mxGetPr(prhs[0]);
  X = (float*) mxGetPr(prhs[1]);
  Y = (float*) mxGetPr(prhs[2]);

  if ((nlhs !=3) && (nlhs !=1))
    mexErrMsgTxt("Must have 1 or 3 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY] )");

  elementsX=(int)mxGetNumberOfElements(prhs[1]);
  elementsY=(int)mxGetNumberOfElements(prhs[2]);

  mwSize dims[2];
  dims[0] = elementsX;
  dims[1] = 1;

  sprintf (buffer, "There is some Error %d elementsX, %d elementsY.\n", elementsX, elementsY);

  if (elementsX != elementsY)
    mexErrMsgTxt(buffer);

  //printf("Elements: img(%d,%d) val: %d \n", N,M, elementsX);
//  plhs[0]  = mxCreateNumericArray(2, dims, mxINT16_CLASS, mxREAL); // the function value
  plhs[0]  = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL); // the function value
  outputF  = (int*) mxGetPr(plhs[0]);

  //BiCubicIp<float> bcIP(N, M, I);
  BiLinearIPInt<float, int> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  bcIP.interpolate_intOmp( elementsX, X, Y, outputF );
  //bcIP.interpolate_int( elementsX, X, Y, outputF );
}
*/
void Interpolate_single ( int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
  int N(0),M(0);

int elementsX, elementsY;
float *outputF, *outputdFX, *outputdFY, *I, *X, *Y;

char buffer[256];

/* check: only one input and one output argument */
if (nrhs < 3)
  mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM Image I; X-Indices X and Y-Indices Y\n)");

N = (int) mxGetN(prhs[0]);
M = (int) mxGetM(prhs[0]);

I = (float*) mxGetPr(prhs[0]);
X = (float*) mxGetPr(prhs[1]);
Y = (float*) mxGetPr(prhs[2]);

if ((nlhs !=3) && (nlhs !=1))
  mexErrMsgTxt("Must have 1 or 3 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY] )");

/* Get the number of elements in the input argument */
elementsX=(int)mxGetNumberOfElements(prhs[1]);
elementsY=(int)mxGetNumberOfElements(prhs[2]);

mwSize dims[2];
dims[0] = elementsX;
dims[1] = 1;

sprintf (buffer, "There is some Error %d elementsX, %d elementsY.\n", elementsX, elementsY);

if (elementsX != elementsY)
  mexErrMsgTxt(buffer);

//printf("Elements: img(%d,%d) val: %d \n", N,M, elementsX);
plhs[0]  = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,mxREAL); // the function value
outputF  = (float*) mxGetPr(plhs[0]);

if (nlhs == 3)
{
  plhs[1]  = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,mxREAL);
//  plhs[1]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative
  outputdFX = (float*) mxGetPr(plhs[1]);
  plhs[2]  = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,mxREAL);
//  plhs[2]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative
  outputdFY = (float*) mxGetPr(plhs[2]);
}


if (nlhs == 3)
{
  BiCubicIp<float> bcIP(N, M, I);
//  BiLinearIp<float> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  bcIP.interpolate( elementsX, X, Y, outputF, outputdFX, outputdFY );
}
else
{
  //BiCubicIp<float> bcIP(N, M, I);
  BiLinearIP<float> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  bcIP.interpolate_Omp( elementsX, X, Y, outputF );
}
}

void Interpolate ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
int N(0),M(0);

int elementsX, elementsY;
double *outputF, *outputdFX, *outputdFY, *I, *X, *Y, *outputdFYY, *outputdFXX, *outputdFXY;

char buffer[256];

/* check: only one input and one output argument */
if (nrhs < 3)
  mexErrMsgTxt("Must have 3 input arguments:\n size of data matrix NxM Image I; X-Indices X and Y-Indices Y\n)");


//N = *mxGetPr(prhs[0]);
//M = *mxGetPr(prhs[1]);
N = (int) mxGetN(prhs[0]);
M = (int) mxGetM(prhs[0]);

//printf("Elements: img(%d,%d) \n", N,M);


I = mxGetPr(prhs[0]);
X = mxGetPr(prhs[1]);
Y = mxGetPr(prhs[2]);

if ((nlhs !=3) && (nlhs !=1) && (nlhs !=6))
  mexErrMsgTxt("Must have 1,3 or 6 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY][dFXX/dYY/dXY] )");

//if (nrhs !=5)
//  mexErrMsgTxt("Must have amount of unknowns 5 input arguments");

/* Create the output array */
/*
for (int i=0;i<U;i++)
{
  plhs[i] = mxCreateDoubleMatrix(N, M, mxREAL);
  output  = mxGetPr(plhs[i]);
  // Populate the output
  memcpy(output, mxGetPr(prhs[i+7]), sizeof(double)*N*M);
}
*/

#ifdef use_class
/* prevent you from passing a sparse matrix,
a string matrix, or a complex array. mxIsComplex
is used to determine if there is an imaginary
part of the mxArray. mxIsClass is used to determine
if the mxArray belongs to a particular class */
//if ( mxIsComplex(prhs[5])|| mxIsClass(prhs[5],"sparse") || mxIsChar(prhs[5]) )
//  mexErrMsgTxt("Input must be real, full, and nonstring");

/* Get the number of elements in the input argument */
elementsX=(int)mxGetNumberOfElements(prhs[1]);
elementsY=(int)mxGetNumberOfElements(prhs[2]);

sprintf (buffer, "There is some Error %d elementsX, %d elementsY.\n", elementsX, elementsY);

if (elementsX != elementsY)
  mexErrMsgTxt(buffer);

//printf("Elements: img(%d,%d) val: %d \n", N,M, elementsX);

if (N*M == elementsX)
{
  plhs[0]  = mxCreateDoubleMatrix(M, N, mxREAL); // the function value
  outputF  = mxGetPr(plhs[0]);

  if (nlhs >= 3)
  {
    plhs[1]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative
    outputdFX = mxGetPr(plhs[1]);
    plhs[2]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative
    outputdFY = mxGetPr(plhs[2]);
  }
  if (nlhs == 6)
  {
    plhs[3]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative dxx
    outputdFXX = mxGetPr(plhs[3]);
    plhs[4]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative dyy
    outputdFYY = mxGetPr(plhs[4]);
    plhs[5]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative dxy
    outputdFXY = mxGetPr(plhs[5]);
  }
}
else
{
  plhs[0]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the function value
  outputF  = mxGetPr(plhs[0]);

  if (nlhs >= 3)
  {
    plhs[1]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative
    outputdFX = mxGetPr(plhs[1]);
    plhs[2]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative
    outputdFY = mxGetPr(plhs[2]);
  }
  if (nlhs == 6)
  {
    plhs[3]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative dxx
    outputdFXX = mxGetPr(plhs[3]);
    plhs[4]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative dyy
    outputdFYY = mxGetPr(plhs[4]);
    plhs[5]  = mxCreateDoubleMatrix(elementsX, 1, mxREAL); // the derivative dxy
    outputdFXY = mxGetPr(plhs[5]);
  }
}


if (nlhs == 6)
{
  BiCubicIp<double> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  bcIP.interpolate_2nd( elementsX, X, Y, outputF, outputdFX, outputdFY, outputdFXX, outputdFYY, outputdFXY );
}
if (nlhs == 3)
{
  BiCubicIp<double> bcIP(N, M, I);
  bcIP.approximateImageGradients();
  bcIP.interpolate( elementsX, X, Y, outputF, outputdFX, outputdFY );
}
//else
if (nlhs == 1)
{
  BiLinearIP<double> bcIP(N, M, I);
//  BiCubicIp<double> bcIP(N, M, I);
  bcIP.approximateImageGradients();
//  bcIP.interpolate_Omp( elementsX, X, Y, outputF );
//  bcIP.interpolate_Omp_unRolled( elementsX, X, Y, outputF );

  bcIP.interpolate_noOmp( elementsX, X, Y, outputF );
//  bcIP.interpolate_noOmp_SSE( elementsX, X, Y, outputF );
  
}
#else
int m, n, i, e, j, nsubs;

plhs[0]  = mxCreateDoubleMatrix(M, N, mxREAL); // the function value
outputF  = mxGetPr(plhs[0]);
plhs[1]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative
outputdFX = mxGetPr(plhs[1]);
plhs[2]  = mxCreateDoubleMatrix(M, N, mxREAL); // the derivative
outputdFY = mxGetPr(plhs[2]);


/* Get the number of dimensions in array */
//nsubs=mxGetNumberOfDimensions(prhs[5]);
//msubs=mxGetNumberOfDimensions(prhs[6]);

// Start by computing derivatives:

//double 
//std::vector<double> Idx(N*M),Idy(N*M),Idxy(N*M);
double *Idx,*Idy,*Idxy;
Idx  = (double*) mxMalloc( N*M*sizeof(double) );
Idy  = (double*) mxMalloc( N*M*sizeof(double) );
Idxy = (double*) mxMalloc( N*M*sizeof(double) );

/////////////////////////////////////////////////////////
// Start OF APPROXIMATING IMAGE GRADIENTS AS PRE-STEP (those are used further below for bi-cubic stuff)
/////////////////////////////////////////////////////////
#pragma omp for //num_threads(omp_get_num_procs())
for(int n=0;n<N;n++)
{
  for(int m=0;m<M;m++)
  {
    int id = n*M+m; // m is the row index, that is determines y, id +1 is row below current one
//    int id_test  = m*N+n;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f\n", n, m, I[id]);
#endif


  if ( (m>0) && (m<M-1) )
    Idy[id]  = (I[id+1]-I[id-1])/2.0;
  if ( (n>0) && (n<N-1) )
    Idx[id]  = (I[id+M]-I[id-M])/2.0;

  if ( (m>0) && (m<M-1) && (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M+1]-I[id+M-1] - I[id-M+1] + I[id-M-1]) /4.0;

  }
}

#ifdef DebugOut
//  mexPrintf("First Row\n");
#endif

// along first row ______
//                 ...
for(int n=0;n<N;n++)
{
  int m = 0;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id+1], I[id]);
#endif

  Idy[id]  = (I[id+1]-I[id])/1.0;

  if ( (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M+1]-I[id+M] - I[id-M+1] + I[id-M]) /2.0;
  else
  {
    if ( (n>0) ) // corner upper right
      Idxy[id] = (I[id+1]-I[id] - I[id-M+1] + I[id-M]) /2.0;
    else // n== 0 // corner upper left
      Idxy[id] = (I[id+M+1]-I[id+M] - I[id+1] + I[id]) /2.0;
  }
}

#ifdef DebugOut
//  mexPrintf("Last Row\n");
#endif
//              ...
// along last row ______
for(int n=0;n<N;n++)
{
  int m = M-1;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, M-1, I[id], I[id-1]);
#endif

  Idy[id]  = (I[id]-I[id-1])/1.0;

  if ( (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M]-I[id+M-1] - I[id-M] + I[id-M-1]) /2.0;
  else
  {
    if ( (n>0) )// corner lower right
      Idxy[id] = (I[id]-I[id-1] - I[id-M] + I[id-M-1]) /2.0;
    else// corner lower left
      Idxy[id] = (I[id+M]-I[id+M-1] - I[id] + I[id-1]) /2.0;
  }
}

#ifdef DebugOut
//  mexPrintf("First Column\n");
#endif
// along first column |...
for(int m=0;m<M;m++)
{
  int n = 0;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id+M], I[id]);
#endif

  Idx[id]  = (I[id+M]-I[id])/1.0;

  if ( (m>0) && (m<M-1) )
    Idxy[id] = (I[id+M+1]-I[id+M-1] - I[id+1] + I[id-1]) /2.0;
  else
  {
    if ( (m>0) )// corner lower left
      Idxy[id] = (I[id+M]-I[id+M-1] - I[id] + I[id-1]) /2.0;
    else// corner upper left
      Idxy[id] = (I[id+M+1]-I[id+M] - I[id+1] + I[id]) /2.0;
  }
}

#ifdef DebugOut
 // mexPrintf("Last Column\n");
#endif
// along last column ...|
for(int m=0;m<M;m++)
{
  int n = N-1;
  int id = n*M+m;

#ifdef DebugOut
 // mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id], I[id-M]);
#endif

  Idx[id]  = (I[id]-I[id-M])/1.0;

  if ( (m>0) && (m<M-1) )
    Idxy[id] = (I[id+1]-I[id-1] - I[id-M+1] + I[id-M-1]) /2.0;
  else
  {
    if ( (m>0) ) // corner lower right
      Idxy[id] = (I[id]-I[id-1] - I[id-M] + I[id-M-1]) /2.0;
    else// corner upper right
      Idxy[id] = (I[id+1]-I[id] - I[id-M+1] + I[id-M]) /2.0;
  }
}
/////////////////////////////////////////////////////////
// END OF APPROXIMATING THE IMMAGE GRADIENTS AS PRE-STEP
/////////////////////////////////////////////////////////
// works
#ifdef DebugOut
for(int n=0;n<N;n++)
{
  for(int m=0;m<M;m++)
  {
    int id = n*M+m; // m is the row index, that is determines y, id +1 is row below current one

    mexPrintf("n=%d m=%d = ( %f %f %f )\n", n, m, Idx[id], Idy[id], Idxy[id]);
  }
}
#endif

#ifdef DebugOut
int nProcs = omp_get_num_procs();
int numThreads = omp_get_num_threads();
int mThreads = omp_get_max_threads();

printf("nProcs: %d, nThreads: %d, mThreads: %d", nProcs, numThreads, mThreads);

mexPrintf("Start Interpolation\n");
#endif

#ifndef _OLD_STYLE_
#pragma omp parallel for schedule (static)
for(int id=0;id<elementsY;id++)
{
  double **cc;
  cc=matrix(1,4,1,4);

  #ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif


#else
#pragma omp parallel for schedule (static)
for(int n=0;n<N;n++)
{
  double **cc;
  cc=matrix(1,4,1,4);

#ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif

  for(int m=0;m<M;m++) // row index is y coordinate
  {
//    int id_inv = n*N+m;
    int id = n*M+m;

//    int id  = m*N+n;
#endif
    int intx, inty;
    double x,y;
    double t,u;

    // c++ style indices 
    double idx = X[id] -1.0;
    double idy = Y[id] -1.0;

    if ( (idx >= double(N-1)) )
      idx = idx - 0.000001;

    if ( (idy >= double(M-1)) )
      idy = idy - 0.000001;

    if ( (idx < double(0)) )
      idx = idx + 0.000001;

    if ( (idy < double(0)) )
      idy = idy + 0.000001;

    x = floor( idx );
    y = floor( idy );

    intx = int(x);
    inty = int(y);

    if ( (intx >= N-1) || (intx < 0) || (inty >= M-1) || (inty < 0) )
    {
      if (intx < 0)   intx = 0;
      if (intx > N-1) intx = N-1;
      if (inty < 0)   inty = 0;
      if (inty > M-1) inty = M-1;
      int myid = intx*M+inty;
      outputF[id]  = I[myid];
      outputdFX[id] = 0.0;
      outputdFY[id] = 0.0;
      continue;
    }

#ifdef DebugOut
    mexPrintf("n=%d m=%d : %d %d  %f - %f\n", n, m, intx, inty, X[id], Y[id]);
    int myid = intx*M+inty+1;
    double a  = I[myid];
    double b = Idx[myid];
    double c = Idy[myid];
    double d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty+1, a, b, c, d);
     myid = intx*M+inty+1+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty+1, a, b, c, d);
    myid = intx*M+inty+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty, a, b, c, d);
     myid = intx*M+inty;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty, a, b, c, d);
#endif

    // first get t and u and the 
//    t = idx-x;
//    u = idy-y;

    double f[5], y1[5], y2[5], y12[5];

    int k = 2, add = 1;
    // start at lower left and go counterclockwise
    int idCC = intx*M + inty + 1;

    k = 4;
    f[k]   =    I[idCC];
    y1[k]  =  Idx[idCC];
    y2[k]  =  Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 3;
    idCC   = intx*M + inty + 1 + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 2;
    idCC   = intx*M + inty + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 1;
    idCC   = intx*M + inty;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

/*
    {
    double d1,d2,**c;
    c=matrix(1,4,1,4);

    d1=idx-x;
    d2=idy-y;
    bcucof(f,y1,y2,y12,d1,d2, c);
    for (i=1;i<=4;i++)
      for (j=1;j<=4;j++)
        mexPrintf(" %f ", c[i][j]);
    mexPrintf("\n");
    }
    */

    bcuint(f, y1, y2, y12, x, x+1.0, y, y+1.0, idx, idy, &(outputF[id]), &(outputdFX[id]), &(outputdFY[id]), cc );
#ifdef _OLD_STYLE_
  }
#endif
  free_matrix(cc,1,4,1,4);
}

mxFree(Idx);
mxFree(Idy);
mxFree(Idxy);
#endif
};

#undef use_class
