// StereoTest.cpp : Defines the entry point for the console application.
//

//#include <cv.h>
//#include <highgui.h>
//#include <cxcore.h>
#include <stdio.h>

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define _matlab_output

#ifdef _matlab_output
#include "mex.h"
#endif

//#define _l1norming_

void get_L1HomoNAllSeg ( int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{
  typedef double Scalar;
  typedef  Math::VectorT<Scalar,3> P3;
  typedef  Math::VectorT<Scalar,2> P2;
  typedef  Math::VectorT<Scalar,3> Vec3;
  typedef  Math::VectorT<Scalar,6> Vec6;
  typedef  Math::Mat3x3T< Scalar > Mat3;
  typedef  Math::Mat3x3T< Scalar > M3;

  Scalar* Ns ( (Scalar*) mxGetPr(prhs[0]) );

  const mxArray* Segments  = (prhs[1]);
  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  Vec3     m ( (Scalar*) mxGetPr(prhs[2]) );
  Scalar* KP  = (Scalar*) mxGetPr(prhs[3]); // the points of left image
  Scalar* MKP = (Scalar*) mxGetPr(prhs[4]); // the fixed part of right points
  Scalar* q2d = (Scalar*) mxGetPr(prhs[5]); // the points to be matched of right image
  Scalar* sigmas= (Scalar*) mxGetPr(prhs[6]); // the sigma of the lorentzian
  int*    valids= (int*) mxGetPr(prhs[7]); // the sigma of the lorentzian

#ifdef _matlab_output___

#ifdef _matlab_output
if (nrhs != 5)
  mexErrMsgTxt("Must have 6 input argument");

if (nlhs < 1)
  mexErrMsgTxt("Must have 1 or 2  output arguments");
#endif
#endif
/*
 int m_q1 =  mxGetM((prhs[2]));
 int m_p  =  mxGetM((prhs[4]));
 int n_q1 =  mxGetN((prhs[2]));
 int n_p  =  mxGetN((prhs[4]));
 */

 plhs[0] = mxCreateDoubleMatrix(1, nSegments, mxREAL);
 Scalar* outF = mxGetPr( plhs[0]) ;

 Scalar* outdF = NULL;
 if (nlhs > 0)
 {
   plhs[1] = mxCreateDoubleMatrix(3, nSegments, mxREAL);
   outdF = mxGetPr( plhs[1]) ;
 }

#pragma omp parallel for
 for (int seg_i = 0; seg_i < nSegments ;seg_i++)
 {
   int* ids = (int*) mxGetPr(mxGetCell(Segments, seg_i));
   int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, seg_i));

   Scalar  F(0);
   Vec3   dF(0,0,0);
   Vec3    N( &Ns[3*seg_i] );

   Scalar sigma = sigmas[seg_i];

   for (int i =0;i < nIds; i++)
   {
     if (!valids[ids[i]] ) continue;

     P3 pi  = P3 (&KP[3*ids[i]]);
     P3 mpi = P3 (&MKP[3*ids[i]]);
     P2 qi  = P2 (&q2d[2*ids[i]]);

     Scalar iDepth = N|pi;
     P3 qhi        = mpi - m*iDepth;

     Scalar rx = qi[0] - qhi[0]/qhi[2];
     Scalar ry = qi[1] - qhi[1]/qhi[2];

     Scalar Z1 = m[0]*qhi[2] - qhi[0] * m[2];
     Scalar Z2 = m[1]*qhi[2] - qhi[1] * m[2];


#ifdef _l1norming_
     // L1 norm: SUCKS
     Scalar isqrtF = 1./sqrt( rx*rx + ry*ry);
     F += 1./isqrtF;

     isqrtF = isqrtF/ (qhi[2]*qhi[2]);
     dF += pi * (rx * isqrtF * Z1 + ry * isqrtF * Z2);

// Lorentzian : works already!
#else
     Scalar r2 = rx*rx + ry*ry;
     F += log(1+ 0.5 * r2/sigma );

     Scalar nen = qhi[2]*qhi[2] * (sigma+r2 * 0.5);
     Scalar dFx = rx * Z1 / nen;
     Scalar dFy = ry * Z2 / nen;

     dF += pi * (dFx+dFy);
#endif

//     if(seg_i==1694)
//       printf("i:%d iDepth%f, rx%f , ry%f, (dFx+dFy)%f sig: %f, pi:%f,%f, qi:%f,%f, mpi:%f,%f\n",i, iDepth, rx, ry, (dFx+dFy), sigma, pi[0], pi[1], qi[0], qi[1], mpi[0], mpi[1]);
   }

   outF[seg_i] = F;

   if (outdF!=NULL)
   {
     outdF[3*seg_i]   = dF[0];
     outdF[3*seg_i+1] = dF[1];
     outdF[3*seg_i+2] = dF[2];
   }
 }

#ifndef _matlab_output

//printf("The end\n");
//printf("Not yet\n");

#endif
}
