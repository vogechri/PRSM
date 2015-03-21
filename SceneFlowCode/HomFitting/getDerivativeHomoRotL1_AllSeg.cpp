// StereoTest.cpp : Defines the entry point for the console application.
//
#include <stdio.h>

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define _matlab_output

#ifdef _matlab_output
#include "mex.h"
#endif


void get_L1HomoAllSeg ( int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;
  typedef  Math::VectorT<Scalar,3> P3;
  typedef  Math::VectorT<Scalar,2> P2;
  typedef  Math::VectorT<Scalar,3> Vec3;
  typedef  Math::VectorT<Scalar,6> Vec6;
  typedef  Math::Mat3x3T< Scalar > Mat3;
  typedef  Math::Mat3x3T< Scalar > M3;

  Scalar* RTs = (Scalar*) mxGetPr(prhs[0]);

  const mxArray* Segments  = (prhs[1]);
  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  Mat3 K  ( (Scalar*) mxGetPr(prhs[2]) );
  Mat3 M  ( (Scalar*) mxGetPr(prhs[3]) );
  Vec3 m  ( (Scalar*) mxGetPr(prhs[4]) );
  //  Mat3 R  ( (Scalar*) mxGetPr(prhs[5]) );
  Scalar* Rs = (Scalar*) mxGetPr(prhs[5]);

  Scalar sigma(0.5);

  Scalar* q1 = (Scalar*) mxGetPr(prhs[6]); // the points to be matched
  Scalar* q2(NULL);

  Scalar*  p      = (Scalar*) mxGetPr(prhs[7]); // the points to be matched
  Scalar*  iDepth = (Scalar*) mxGetPr(prhs[8]); // the inverse depth of the points

  int doL1 = 1;
  if (nrhs > 9)
    doL1 = (int) *mxGetPr(prhs[9]); // the points to be matched

  Scalar* sigmas=NULL;
  if (nrhs > 10)
    sigmas= (Scalar*) mxGetPr(prhs[10]); // the sigma of the lorentzian

  if (nrhs > 11)
    q2 = (Scalar*) mxGetPr(prhs[11]); // the points to be matched


#ifdef _matlab_output___

#ifdef _matlab_output
  if (nrhs < 9 || nrhs > 11)
    mexErrMsgTxt("Must have 9, 10 or 11 input argument");

  if (nlhs < 1)
    mexErrMsgTxt("Must have 1 or 2  output arguments");
#endif
#endif
  //Edge-Weights W; Righthand Side B; Vertex Scales S; Init Solution D
  //if (nrhs !=U+7)
  //  mexErrMsgTxt("Must have amount of unknowns (U)+7 input arguments");

  /* Create the output array */
  {
#ifdef _matlab_output
    //  plhs[0] = mxCreateDoubleMatrix(N*M, 3, mxREAL);
    //  output_Points  = mxGetPr(plhs[0]);
#else
    //  output_Points = (double*) malloc(sizeof(double) * N*M);
#endif
    /* Populate the output */
    //  memcpy(output, mxGetPr(prhs[9]), sizeof(double)*N*M);

    //  plhs[1] = mxCreateDoubleMatrix(N, M, mxREAL);
    //  output_dF  = mxGetPr(plhs[1]);
    /* Populate the output */
    //  memcpy(output_dF, mxGetPr(prhs[9]), sizeof(double)*N*M);
  }

  /* prevent you from passing a sparse matrix,
  a string matrix, or a complex array. mxIsComplex
  is used to determine if there is an imaginary
  part of the mxArray. mxIsClass is used to determine
  if the mxArray belongs to a particular class */
  //if ( mxIsComplex(prhs[9])|| mxIsClass(prhs[9],"sparse") || mxIsChar(prhs[9]) )
  //  mexErrMsgTxt("Input must be real, full, and nonstring");

  /* Get the number of elements in the input argument */
#ifdef _matlab_output

  //elements=mxGetNumberOfElements(prhs[3]);

  //sprintf (buffer, "The amount of elements in the matrix is not equal to the amount given by the specified dimensions!\n all: %d product %d   N: %d, M: %d\n", elements, N*M, N, M);
  //if (elements != N*M)
  //  mexErrMsgTxt(buffer);
#endif

  plhs[0] = mxCreateDoubleMatrix(1, nSegments, mxREAL);
  Scalar* outF = mxGetPr( plhs[0]) ;

  Scalar* outdF = NULL;
  if (nlhs > 0)
  {
    plhs[1] = mxCreateDoubleMatrix(6, nSegments, mxREAL);
    outdF = mxGetPr( plhs[1]) ;
  }

  // first compute new rotation
#pragma omp parallel for
  for (int seg_i = 0; seg_i < nSegments ;seg_i++)
  {
    int* ids = (int*) mxGetPr(mxGetCell(Segments, seg_i));
    int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, seg_i));

    Scalar  F(0);
    Vec6   dF(0,0,0, 0,0,0);

    Vec6 RT(&RTs[6*seg_i]);
    Mat3 R (&Rs [9*seg_i]);

    P3 ri_x = P3 (&( RT[ 0 ] ));
    P3 ti_x = P3 (&( RT[ 3 ] ));

    Scalar sinA = ri_x.norm();
    sinA = std::min(sinA, 1.);sinA = std::max(sinA, 0.000000001);

    P3 rVec = ri_x * (1./sinA);
    Scalar alpha = asin(sinA);
    Scalar cosA = cos(alpha);
    //    Scalar cosA = sqrt(1. - sinA*sinA);// sign anyone ? (|alpha| > pi/2
    Scalar cosA_1 = Scalar(1) - cosA;

    M3 A1( cosA,   -ri_x[2],  ri_x[1],
      ri_x[2], cosA,    -ri_x[0],
      -ri_x[1], ri_x[0],  cosA );

    M3 A2( cosA_1 * rVec[0]*rVec[0], cosA_1 * rVec[0]*rVec[1] ,cosA_1 * rVec[0]*rVec[2],
      cosA_1 * rVec[0]*rVec[1], cosA_1 * rVec[1]*rVec[1], cosA_1 * rVec[2]*rVec[1],
      cosA_1 * rVec[0]*rVec[2], cosA_1 * rVec[2]*rVec[1], cosA_1 * rVec[2]*rVec[1] );

    M3 R_new = R * (A1+A2);
    //////////////////////////////////////

    M3   KR = K*R_new;
    Vec3 Kt = K*ti_x;

    Scalar sigma = 1;
    if (sigmas !=NULL) sigma=sigmas[seg_i];

    Scalar badRotPenalty = 0;
    if(sinA>0.5) // 30 degrees here
      badRotPenalty = 10000; // rotation too large - really

    for (int i =0;i < nIds; i++)
    {
      //     if (!valids[ids[i]] ) continue;

      int id = ids[i];
      P3 pi = P3 (&p[3*id]);
      P2 qi = P2 (&q1[2*id]);

      P3 krkp  = KR * pi;
      M3 KRKpx = KR*  M3 ( 0, -pi[2],  pi[1],   pi[2], 0, -pi[0],  -pi[1], pi[0],  0 );

      P3 Kdt = Kt * iDepth[id];
      M3 Kd  = K*iDepth[id];

      Vec6 amt1 = Vec6( KRKpx(0, 0), KRKpx(0, 1), KRKpx(0, 2), Kd(0, 0), Kd(0, 1), Kd(0, 2) );
      Vec6 amt2 = Vec6( KRKpx(1, 0), KRKpx(1, 1), KRKpx(1, 2), Kd(1, 0), Kd(1, 1), Kd(1, 2) );
      Vec6 amt3 = Vec6( KRKpx(2, 0), KRKpx(2, 1), KRKpx(2, 2), Kd(2, 0), Kd(2, 1), Kd(2, 2) );

      Vec3 phi = -krkp + Kdt;

      Scalar rx = phi[0]/phi[2] - qi[0];
      Scalar ry = phi[1]/phi[2] - qi[1];

//      if(seg_i==1694)
//        printf("i:%d iDepth%f, rx%f , ry%f, pi:%f,%f, qi:%f,%f, ri_x, ti_x:%f%f%f ,%f%f%f\n",i, iDepth[id], rx, ry, pi[0], pi[1], qi[0], qi[1], ri_x[0], ri_x[1], ri_x[2], ti_x[0], ti_x[1], ti_x[2] );


      if (doL1)
      {
        Scalar isqrtF = 1./(sqrt( rx*rx + ry*ry)+0.0001);
        F += 1./isqrtF;
        dF += ( rx * isqrtF / (phi[2]*phi[2]) )* ( amt1 * phi[2] - amt3 * phi[0]  ) +
          ( ry * isqrtF / (phi[2]*phi[2]) )* ( amt2 * phi[2] - amt3 * phi[1]  );
      }
      else // L2 or lorentzian
      {
        /*
        F += rx*rx + ry*ry;
        dF +=
        ( 2.*rx / (phi[2]*phi[2]) ) * ( amt1 * phi[2] - amt3 * phi[0] ) +
        ( 2.*ry / (phi[2]*phi[2]) ) * ( amt2 * phi[2] - amt3 * phi[1] );     
        */

        Scalar r2 = rx*rx + ry*ry;
        F += log(1+ 0.5 * r2/sigma );

        Scalar nen = phi[2]*phi[2] * (sigma+r2 * 0.5);

        dF +=  ( rx / nen) * ( amt1 * phi[2] - amt3 * phi[0] ) +
          ( ry / nen) * ( amt2 * phi[2] - amt3 * phi[1] );

      }
    }

    if (q2 != NULL)
    {
      M3 MR = M*R_new;
      Vec3 Mt = M*ti_x;
      for (int i =0;i < nIds; i++)
      {
        //     if (!valids[ids[i]] ) continue;

        int id = ids[i];
        P3 pi = P3 (&p[3*id]);
        P2 qi = P2 (&q2[2*id]);

        P3 krkp  = MR * pi - m * iDepth[id];
        M3 KRKpx = MR*  M3 ( 0, -pi[2],  pi[1],   pi[2], 0, -pi[0],  -pi[1], pi[0],  0 );

        P3 Kdt   = Mt * iDepth[id];
        M3 Kd    = M*iDepth[id];

        Vec6 amt1 = Vec6( KRKpx(0, 0), KRKpx(0, 1), KRKpx(0, 2), Kd(0, 0), Kd(0, 1), Kd(0, 2) );
        Vec6 amt2 = Vec6( KRKpx(1, 0), KRKpx(1, 1), KRKpx(1, 2), Kd(1, 0), Kd(1, 1), Kd(1, 2) );
        Vec6 amt3 = Vec6( KRKpx(2, 0), KRKpx(2, 1), KRKpx(2, 2), Kd(2, 0), Kd(2, 1), Kd(2, 2) );

        Vec3 phi = -krkp + Kdt;

        Scalar rx = phi[0]/phi[2] - qi[0];
        Scalar ry = phi[1]/phi[2] - qi[1];

        if (doL1)
        {
          Scalar isqrtF = 1./(sqrt( rx*rx + ry*ry)+0.0001);
          F += 1./isqrtF;
          dF +=
            ( rx * isqrtF / (phi[2]*phi[2]) ) * ( amt1 * phi[2] - amt3 * phi[0] ) +
            ( ry * isqrtF / (phi[2]*phi[2]) ) * ( amt2 * phi[2] - amt3 * phi[1] );
        }
        else // L2
        {
          /*
          F += rx*rx + ry*ry;
          dF +=
          ( 2.*rx / (phi[2]*phi[2]) ) * ( amt1 * phi[2] - amt3 * phi[0] ) +
          ( 2.*ry / (phi[2]*phi[2]) ) * ( amt2 * phi[2] - amt3 * phi[1] );  
          */

          Scalar r2 = rx*rx + ry*ry;
          F += log(1+ 0.5 * r2/sigma );

          Scalar nen = phi[2]*phi[2] * (sigma+r2 * 0.5);

          dF +=
            ( rx / nen ) * ( amt1 * phi[2] - amt3 * phi[0] ) +
            ( ry / nen ) * ( amt2 * phi[2] - amt3 * phi[1] );
        }
      }
    }
    //////////////////////

    outF[seg_i] = F+badRotPenalty;

    if (nlhs > 0)
    {
      outdF[6*seg_i+0] = dF[0];
      outdF[6*seg_i+1] = dF[1];
      outdF[6*seg_i+2] = dF[2];
      outdF[6*seg_i+3] = dF[3];
      outdF[6*seg_i+4] = dF[4];
      outdF[6*seg_i+5] = dF[5];
    }
  }

#ifndef _matlab_output

  //printf("The end\n");
  //printf("Not yet\n");

#endif
}
