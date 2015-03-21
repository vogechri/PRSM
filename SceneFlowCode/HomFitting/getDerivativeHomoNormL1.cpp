// StereoTest.cpp : Defines the entry point for the console application.
//
#include <stdio.h>

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define _matlab_output

#ifdef _matlab_output
#include "mex.h"
#endif


void get_L1HomoN ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  typedef double Scalar;
  typedef  Math::VectorT<Scalar,3> P3;
  typedef  Math::VectorT<Scalar,2> P2;
  typedef  Math::VectorT<Scalar,3> Vec3;
  typedef  Math::VectorT<Scalar,6> Vec6;
  typedef  Math::Mat3x3T< Scalar > Mat3;
  typedef  Math::Mat3x3T< Scalar > M3;

  Vec3 N ( (Scalar*) mxGetPr(prhs[0]) );
  Vec3 m  ( (Scalar*) mxGetPr(prhs[1]) );
  Scalar* KP  = (Scalar*) mxGetPr(prhs[2]); // the points of left image
  Scalar* MKP = (Scalar*) mxGetPr(prhs[3]); // the fixed part of right points
  Scalar* q2d = (Scalar*) mxGetPr(prhs[4]); // the points to be matched of right image
  Scalar sigma= (Scalar) *mxGetPr(prhs[5]); // the sigma of the lorentzian

#ifdef _matlab_output___

#ifdef _matlab_output
if (nrhs != 5)
  mexErrMsgTxt("Must have 5 input argument");

if (nlhs < 1)
  mexErrMsgTxt("Must have 1 or 2  output arguments");
#endif
#endif

 int m_q1 =  mxGetM((prhs[2]));
 int m_p  =  mxGetM((prhs[4]));
 int n_q1 =  mxGetN((prhs[2]));
 int n_p  =  mxGetN((prhs[4]));


 Scalar F(0);
 Vec3  dF(0,0,0);

 for (int i =0;i < n_p; i++)
 {
   P3 pi  = P3 (&KP[3*i]);
   P3 mpi = P3 (&MKP[3*i]);
   P2 qi  = P2 (&q2d[2*i]);

   Scalar iDepth = N|pi;
   P3 qhi        = mpi - m*iDepth;

  Scalar rx = qi[0] - qhi[0]/qhi[2];
  Scalar ry = qi[1] - qhi[1]/qhi[2];

  Scalar Z1 = m[0]*qhi[2] - qhi[0] * m[2];
  Scalar Z2 = m[1]*qhi[2] - qhi[1] * m[2];

  Scalar r2 = rx*rx + ry*ry;
  F += log(1+ 0.5 * r2/sigma );

  Scalar nen = qhi[2]*qhi[2] * (sigma+r2 * 0.5);
  Scalar dFx = rx * Z1 / nen;
  Scalar dFy = ry * Z2 / nen;

  dF += pi * (dFx+dFy);

//  printf("i:%d iDepth%f, rx%f , ry%f, (dFx+dFy)%f sig: %f, pi:%f,%f, qi:%f,%f\n",i, iDepth, rx, ry, (dFx+dFy), sigma, pi[0], pi[1], qi[0], qi[1]);
 }

 plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
 Scalar* outF = mxGetPr( plhs[0]) ;
 *outF = F;

 if (nlhs > 0)
 {
   plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
   Scalar* outdF = mxGetPr( plhs[1]) ;
   outdF[0] = dF[0];
   outdF[1] = dF[1];
   outdF[2] = dF[2];
 }


#ifndef _matlab_output

//printf("The end\n");
//printf("Not yet\n");

#endif
}
