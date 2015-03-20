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
#ifndef _BILINEAR_IP_H
#define _BILINEAR_IP_H

#include <math.h>
#include<vector>

//#define _NO_OPENMP

#ifndef _NO_OPENMP
#include <omp.h>
#endif

#undef DebugOut

template <typename S_>
class BiLinearIP
{
public:

  typedef S_ Scalar;

  BiLinearIP(int _N, int _M, Scalar*& _Img ) 
    : N(_N), M(_M), I(_Img), maxConstN(Scalar(_N)-1.000001), maxConstM(Scalar(_M)-1.000001)
  {};

  ~BiLinearIP() {};

  void interpolate_noOmp( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF )
  {
    for(int id=0;id<elements;id++)
    {
//      int intx, inty;
//      Scalar x,y;

      // c++ style indices not matlab
      Scalar idx = X[id] -Scalar(1.0);
      Scalar idy = Y[id] -Scalar(1.0);

      // assume it is a bit above, then we have >= true, so no interpolation
      // now assume exact match, still no interpol == with interpolation
      /*
      if ( (idx >= Scalar(N-1)) )
        idx = idx - Scalar(0.0001);

      if ( (idy >= Scalar(M-1)) )
        idy = idy - Scalar(0.0001);

      if ( (idx < Scalar(0)) )
        idx = idx + Scalar(0.0001);

      if ( (idy < Scalar(0)) )
        idy = idy + Scalar(0.0001);
      
      // cpp rounds down anyway:
      //      x = floor( idx );
      //      y = floor( idy );

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
        continue;
      }
      */
      /*
      // so ceil(0.00000000000001) = 1, ceil(N-1)=N-1
      idx = max( Scalar(0.00000001), min(idx, Scalar(N-1)) );
      idy = max( Scalar(0.00000001), min(idy, Scalar(M-1)) );


//      Scalar x1 = floor(idx);
      const Scalar x2 = ceil(idx);
//      Scalar y1 = floor(idy);
      const Scalar y2 = ceil(idy);

      const Scalar a = (x2-idx);
      const Scalar b = (y2-idy);
      const Scalar W11 = a*b;
      const Scalar W21 = (1-b)*a;
      const Scalar W12 = b*(1-a);
      const Scalar W22 = (1-W11-W21-W12);

//      int idCC   =   (int)(x2*M+y2);
      const Scalar* p22 = &(I[(int)(x2*M+y2)]);

      const Scalar& I22 = p22[0];
      const Scalar& I12 = p22[-1];
      const Scalar& I21 = p22[-M];
      const Scalar& I11 = p22[-M-1];     
      */
/*
      int idCC   =   (int)(x2*M+y2);
      const Scalar& I22 = I[idCC];
      const Scalar& I12 = I[idCC-1];
      const Scalar& I21 = I[idCC-M];
      const Scalar& I11 = I[idCC-M-1];     
*/

      idx = std::max( Scalar(0.000001), std::min(idx, maxConstN) );
      idy = std::max( Scalar(0.000001), std::min(idy, maxConstM) );
     
      const int x2 = (int)(idx);// thats floor
      const int y2 = (int)(idy);

      const Scalar a = (idx-(Scalar)(x2));
      const Scalar b = (idy-(Scalar)(y2));
      const Scalar W22 = a*b;
      const Scalar W12 = (1-b)*a;
      const Scalar W21 = b*(1-a);
      const Scalar W11 = (1-W22-W21-W12);

      const Scalar* p11 = &(I[x2*M+y2]);

      const Scalar& I11 = p11[0];     
      const Scalar& I21 = p11[1];
      const Scalar& I12 = p11[M];
      const Scalar& I22 = p11[M+1];

      outputF[id] = W11*I11 + W12*I12 + W21*I21 + W22*I22;
    }
  }


  void interpolate_Omp( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF )
  {
#pragma omp parallel for schedule (static)
    for(int id=0;id<elements;id++)
    {
//      int intx, inty;
//      Scalar x,y;

      // c++ style indices not matlab
      Scalar idx = X[id] -Scalar(1.0);
      Scalar idy = Y[id] -Scalar(1.0);
      /*
      // so ceil(0.00000000000001) = 1, ceil(N-1)=N-1
      idx = std::max( Scalar(0.00000001), std::min(idx, Scalar(N-1)) );
      idy = std::max( Scalar(0.00000001), std::min(idy, Scalar(M-1)) );

      const Scalar x2 = ceil(idx);
      const Scalar y2 = ceil(idy);

      const Scalar a = (x2-idx);
      const Scalar b = (y2-idy);
      const Scalar W11 = a*b;
      const Scalar W21 = (1-b)*a;
      const Scalar W12 = b*(1-a);
      const Scalar W22 = (1-W11-W21-W12);

      const Scalar* p22 = &(I[(int)(x2*M+y2)]);

      const Scalar& I22 = p22[0];
      const Scalar& I12 = p22[-1];
      const Scalar& I21 = p22[-M];
      const Scalar& I11 = p22[-M-1];     
      */

      idx = std::max( Scalar(0.000001), std::min(idx, maxConstN) );
      idy = std::max( Scalar(0.000001), std::min(idy, maxConstM) );
     
      const int x2 = (int)(idx);// thats floor
      const int y2 = (int)(idy);

      const Scalar a = (idx-(Scalar)(x2));
      const Scalar b = (idy-(Scalar)(y2));
      const Scalar W22 = a*b;
      const Scalar W12 = (1-b)*a;
      const Scalar W21 = b*(1-a);
      const Scalar W11 = (1-W22-W21-W12);

      const Scalar* p11 = &(I[x2*M+y2]);

      const Scalar& I11 = p11[0];     
      const Scalar& I21 = p11[1];
      const Scalar& I12 = p11[M];
      const Scalar& I22 = p11[M+1];

      outputF[id] = W11*I11 + W12*I12 + W21*I21 + W22*I22;
    }
  }

  void interpolate_Omp_unRolled( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF )
  {
    const int unrolC = 4;

    int step4 = elements/unrolC;
    int last  = unrolC*(step4);

#pragma omp parallel for schedule (static)
    for(int id=0;id<elements;id+=unrolC)
    {
      // c++ style indices not matlab
      Scalar idx1 = X[id]  ;
      Scalar idx2 = X[id+1];
      Scalar idx3 = X[id+2];
      Scalar idx4 = X[id+3];

      Scalar idy1 = Y[id]  ;
      Scalar idy2 = Y[id+1];
      Scalar idy3 = Y[id+2];
      Scalar idy4 = Y[id+3];

      idx1 = std::max( Scalar(0.00000001), std::min(idx1-Scalar(1.0), maxConstN) );
      idx2 = std::max( Scalar(0.00000001), std::min(idx2-Scalar(1.0), maxConstN) );
      idx3 = std::max( Scalar(0.00000001), std::min(idx3-Scalar(1.0), maxConstN) );
      idx4 = std::max( Scalar(0.00000001), std::min(idx4-Scalar(1.0), maxConstN) );
      idy1 = std::max( Scalar(0.00000001), std::min(idy1-Scalar(1.0), maxConstM) );
      idy2 = std::max( Scalar(0.00000001), std::min(idy2-Scalar(1.0), maxConstM) );
      idy3 = std::max( Scalar(0.00000001), std::min(idy3-Scalar(1.0), maxConstM) );
      idy4 = std::max( Scalar(0.00000001), std::min(idy4-Scalar(1.0), maxConstM) );

      const Scalar x2_1 = (Scalar)((int)(idx1));
      const Scalar x2_2 = (Scalar)((int)(idx2));
      const Scalar x2_3 = (Scalar)((int)(idx3));
      const Scalar x2_4 = (Scalar)((int)(idx4));

      const Scalar y2_1 = (Scalar)((int)(idy1));
      const Scalar y2_2 = (Scalar)((int)(idy2));
      const Scalar y2_3 = (Scalar)((int)(idy3));
      const Scalar y2_4 = (Scalar)((int)(idy4));

      // load ahead: works indeed
      const Scalar* p11_1 = &(I[(int)(x2_1*M+y2_1)]);
      const Scalar* p11_2 = &(I[(int)(x2_2*M+y2_2)]);
      const Scalar* p11_3 = &(I[(int)(x2_3*M+y2_3)]);
      const Scalar* p11_4 = &(I[(int)(x2_4*M+y2_4)]);

      const Scalar a1 = (idx1-x2_1);
      const Scalar a2 = (idx2-x2_2);
      const Scalar a3 = (idx3-x2_3);
      const Scalar a4 = (idx4-x2_4);

      const Scalar b1 = (idy1-y2_1);
      const Scalar b2 = (idy2-y2_2);
      const Scalar b3 = (idy3-y2_3);
      const Scalar b4 = (idy4-y2_4);


      const Scalar W22_1 = a1*b1;
      const Scalar W22_2 = a2*b2;
      const Scalar W22_3 = a3*b3;
      const Scalar W22_4 = a4*b4;

      const Scalar W12_1 = a1*(Scalar(1)-b1);
      const Scalar W12_2 = a2*(Scalar(1)-b2);
      const Scalar W12_3 = a3*(Scalar(1)-b3);
      const Scalar W12_4 = a4*(Scalar(1)-b4);

      const Scalar W21_1 = b1*(Scalar(1)-a1);
      const Scalar W21_2 = b2*(Scalar(1)-a2);
      const Scalar W21_3 = b3*(Scalar(1)-a3);
      const Scalar W21_4 = b4*(Scalar(1)-a4);

      const Scalar W11_1 = (Scalar(1)-W22_1-W21_1-W12_1);
      const Scalar W11_2 = (Scalar(1)-W22_2-W21_2-W12_2);
      const Scalar W11_3 = (Scalar(1)-W22_3-W21_3-W12_3);
      const Scalar W11_4 = (Scalar(1)-W22_4-W21_4-W12_4);


      outputF[id]   = W11_1*p11_1[0] + W21_1*p11_1[1] + W12_1*p11_1[M] + W22_1*p11_1[M+1];
      outputF[id+1] = W11_2*p11_2[0] + W21_2*p11_2[1] + W12_2*p11_2[M] + W22_2*p11_2[M+1];
      outputF[id+2] = W11_3*p11_3[0] + W21_3*p11_3[1] + W12_3*p11_3[M] + W22_3*p11_3[M+1];
      outputF[id+3] = W11_4*p11_4[0] + W21_4*p11_4[1] + W12_4*p11_4[M] + W22_4*p11_4[M+1];
    }

    for(int id = last;id<elements;id++)
    {
      // c++ style indices not matlab
      Scalar idx = X[id]-Scalar(1.0);
      Scalar idy = Y[id]-Scalar(1.0);

      idx = std::max( Scalar(0.000001), std::min(idx, maxConstN) );
      idy = std::max( Scalar(0.000001), std::min(idy, maxConstM) );

      const Scalar x2 = (Scalar)((int)(idx));// thats floor
      const Scalar y2 = (Scalar)((int)(idy));

      const Scalar a = (idx-x2);
      const Scalar b = (idy-y2);
      const Scalar W22 = a*b;
      const Scalar W12 = (1-b)*a;
      const Scalar W21 = b*(1-a);
      const Scalar W11 = (1-W22-W21-W12);

      const Scalar* p11 = &(I[(int)(x2*M+y2)]);

      const Scalar& I11 = p11[0];     
      const Scalar& I21 = p11[1];
      const Scalar& I12 = p11[M];
      const Scalar& I22 = p11[M+1];

      outputF[id] = W11*I11 + W12*I12 + W21*I21 + W22*I22;
    }
  }

  /// compute interpolation and store in the respective containers
  void interpolateDull( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF )
  {

#pragma omp parallel for schedule (static)
    for(int id=0;id<elements;id++)
    {
#ifdef DebugOut
#ifndef _NO_OPENMP
      int numPrcs = omp_get_num_procs();
      int th_id = omp_get_thread_num();
      printf("row %d threadid is : %d\n", n, th_id);
#endif
#endif

      int tNum(0); 
#ifndef _NO_OPENMP
      tNum = omp_get_thread_num( ); // the current thread, defines where to store the temp per thread solutuion
#endif
      int intx, inty;
      Scalar x,y;

      // c++ style indices 
      Scalar idx = X[id] -Scalar(1.0);
      Scalar idy = Y[id] -Scalar(1.0);

      if ( (idx >= Scalar(N-1)) )
        idx = idx - Scalar(0.0001);

      if ( (idy >= Scalar(M-1)) )
        idy = idy - Scalar(0.0001);

      if ( (idx < Scalar(0)) )
        idx = idx + Scalar(0.0001);

      if ( (idy < Scalar(0)) )
        idy = idy + Scalar(0.0001);

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
        continue;
      }

      Scalar x1 = floor(idx);
      Scalar x2 = ceil(idx);
      Scalar y1 = floor(idy);
      Scalar y2 = ceil(idy);

      Scalar a = (x2-idx);
      Scalar b = (y2-idy);
      Scalar W11 = a*b;
      Scalar W21 = (1-b)*a;
      Scalar W12 = b*(1-a);
      Scalar W22 = (1-a)*(1-b);

      int idCC = x1*M+y1;
      Scalar I11 = I[idCC];

      idCC   = x2*M+y1;
      Scalar I12 = I[idCC];

      idCC   =  x2*M+y2;
      Scalar I22 = I[idCC];

      idCC   =  x1*M+y2;
      Scalar I21 = I[idCC];

      outputF[id] = W11*I11 + W12*I12 + W21*I21 + W22*I22;
    }

  }

  void approximateImageGradients() {};

private:

  int N;
  int M;

  const Scalar maxConstN;
  const Scalar maxConstM;

  const Scalar* I;
};


template <typename S_, typename I_>
class BiLinearIPInt
{
public:

  typedef S_ Scalar;
  typedef I_ intType;

  BiLinearIPInt(int _N, int _M, intType*& _Img ) 
    : N(_N), M(_M), I(_Img), maxConstN(Scalar(_N)-1.000001), maxConstM(Scalar(_M)-1.000001)
  {};

  ~BiLinearIPInt() {};

// lacks clamp borders ! do outside
inline intType GetPixel(const intType* img, Scalar x, Scalar y)
{
 unsigned long Fx = (long)(x * 65536.0); // convert to Fixed16
 unsigned long Fy = (long)(y * 65536.0); // convert to Fixed16
 unsigned long px = (Fx & 0xFFFF0000) >> 16; // floor function
 unsigned long py = (Fy & 0xFFFF0000) >> 16; // floor function

 // all 32 bits, so pixel has values between 0 and 65535
 // now a one is 0x10000, mult by 1 is shift by 16 bits = 4 stellen, so after multiplying we shift by 16
 const unsigned int stride = M;

// int testx  = px >> 16;
// int testy  = py >> 16;

 // position = x2*M+y2
 const intType* p0 = img + px * stride + py; // pointer to first pixel
// const intType* p0 = img + (px>>16) * stride + (py>>16); // pointer to first pixel

 // load the four neighboring pixels
 const intType& p1 = p0[0];
 const intType& p3 = p0[1];
 const intType& p2 = p0[0 + stride];
 const intType& p4 = p0[1 + stride];

 // Calculate the weights for each pixel
 unsigned long fx = Fx & 0x0000FFFF; // frac function
 unsigned long fy = Fy & 0x0000FFFF; // frac function
 unsigned long fx1 = 0x00010000 - fx; // 1 - fx
 unsigned long fy1 = 0x00010000 - fy; // 1 - fy;

// long w1 = (fx1 * fy1) >> 16;
 intType w2 = (fx * fy1)  >> 16;
 intType w3 = (fx1 * fy ) >> 16;
 intType w4 = (fx * fy )  >> 16;
 intType w1 = 0x00010000 - (w4+w2+w3);

 //intType sumW = w1+w2+w3+w4;
 //intType temp = 0;
 //if (sumW != 0x00010000 )
 //  temp++;
 //if (w2 < 0 && w4 < 0 && w3 < 0 && w1 < 0)
 //  temp++;

 // Calculate the weighted sum of pixels (for each color channel)
 intType outr = (intType)((p1 * w1 + p2 * w2 + p3 * w3 + p4 * w4) >> 16);

 return outr;
// Pixel(outr >> 8, outg >> 8, outb >> 8, outa >> 8);
}

  void interpolate_int( int elements, Scalar*& X, Scalar *&Y, intType *&outputF )
  {
    for(int id=0;id<elements;id++)
    {
      // c++ style indices not matlab
      Scalar idx = X[id] -Scalar(1.0);
      Scalar idy = Y[id] -Scalar(1.0);

      idx = std::max( Scalar(0.00000001), std::min(idx, Scalar(N-1)) );
      idy = std::max( Scalar(0.00000001), std::min(idy, Scalar(M-1)) );

      outputF[id] = GetPixel(I, idx, idy);
    }
  }

  void interpolate_intOmp( int elements, Scalar*& X, Scalar *&Y, intType *&outputF )
  {
#pragma omp parallel for schedule (static)
    for(int id=0;id<elements;id++)
    {
      // c++ style indices not matlab
      Scalar idx = X[id] -Scalar(1.0);
      Scalar idy = Y[id] -Scalar(1.0);

      idx = std::max( Scalar(0.00000001), std::min(idx, Scalar(N-1)) );
      idy = std::max( Scalar(0.00000001), std::min(idy, Scalar(M-1)) );

      outputF[id] = GetPixel(I, idx, idy);
    }
  }

  void approximateImageGradients() {};

private:

  int N;
  int M;

  const Scalar maxConstN;
  const Scalar maxConstM;

  const intType* I;
};

#endif

