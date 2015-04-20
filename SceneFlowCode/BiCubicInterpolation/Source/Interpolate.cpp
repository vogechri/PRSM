/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

This software is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this software.
*/

#include <math.h>
#include "bcuint.cpp"

#include "BicubicIp.h"
#include "BilinearIp.h"

#define NRANSI
#include "nrutil.h"
#undef NRANSI

#define use_class

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
  outputdFX = (float*) mxGetPr(plhs[1]);
  plhs[2]  = mxCreateNumericArray(2, dims, mxSINGLE_CLASS,mxREAL);
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

N = (int) mxGetN(prhs[0]);
M = (int) mxGetM(prhs[0]);

I = mxGetPr(prhs[0]);
X = mxGetPr(prhs[1]);
Y = mxGetPr(prhs[2]);

if ((nlhs !=3) && (nlhs !=1) && (nlhs !=6))
  mexErrMsgTxt("Must have 1,3 or 6 output arguments: ( F(x) [and dF(x)/dX and dF(x)/dY][dFXX/dYY/dXY] )");

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


//if (nlhs == 6)
//{
//  BiCubicIp<double> bcIP(N, M, I);
//  bcIP.approximateImageGradients();
//  bcIP.interpolate_2nd( elementsX, X, Y, outputF, outputdFX, outputdFY, outputdFXX, outputdFYY, outputdFXY );
//}
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
  bcIP.approximateImageGradients();
  bcIP.interpolate_noOmp( elementsX, X, Y, outputF );
}
};
