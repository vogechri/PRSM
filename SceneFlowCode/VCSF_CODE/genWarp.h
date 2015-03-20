/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
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

#ifndef __GEN_WARP_HH
#define __GEN_WARP_HH

#include "DataDefinitionsVC.h"
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "BilinearIp.h"

using namespace std;
using namespace Math;

// 7x7 census
//#define boxRadius 3 


template<typename Scalar> class genWarpS
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<int, 4>     P4i;
  
  genWarpS( int N_, int M_, Scalar*& img2_)
    : nPixel(N_*M_), nSegments(0), nNormals(0), N(N_), M(M_), dummyWarp(0), bcIP(N_, M_, img2_), 
      Idx(NULL), Idy(NULL), img2(img2_), i_out(NULL), i_in_short(NULL), i_out_short(NULL), i_inout_short(NULL)
  {if (N_*M_ > 0) {init();generateShorts(img2_);};};

  ~genWarpS(){finish();deleteShorts();};

  void reInit(int N_, int M_, Scalar*& img2_)
  {
    finish();deleteShorts();
    img2=img2_;
    if (N_*M_ > 0) {init();generateShorts(img2_);};
  }

  /// idea is to use shorts instead of double to be faster and instead of byte to be more precise (fixed point representation)
  void generateShorts(Scalar*& img2_)
  {
    // increase to prevent looking up bad memory locations
    if (i_in_short  == NULL)
      i_in_short  = (short*) malloc( (N*M+M+1)*sizeof(short) );

    if (i_out_short == NULL)
      i_out_short = (short*) malloc( N*M*sizeof(short) );
    if (i_inout_short == NULL)
      i_inout_short = (short*) malloc( N*M*sizeof(short) );

    for (int ij = 0; ij < N*M;ij++)
    {
      i_in_short [ij]    =       (short) ( img2_[ij]*Scalar(255.0) );
      i_out_short [ij]   = 128 * (short) ( img2_[ij]*Scalar(255.0) );
      i_inout_short[ij]  = 128 * (short) ( img2_[ij]*Scalar(255.0) );
    }
    // prevent bad lookups in bilinearIP
    for (int ij = (N-1)*M; ij < N*M;ij++)
      i_in_short [ij+M]  = (short) ( img2_[ij]*Scalar(255.0) );
    i_in_short [N*M+M] = (short) ( img2_[N*M]*Scalar(255.0) );

  }

  void deleteShorts()
  {
    if (i_in_short  != NULL)
       free( i_in_short );
    if (i_out_short != NULL)
       free( i_out_short );
    if (i_inout_short != NULL)
       free( i_inout_short );
    i_in_short   = NULL;
    i_out_short  = NULL;
    i_inout_short  = NULL;
  }

  void init()
  {
    if (Idx == NULL)
      Idx  = (Scalar*) malloc( N*M*sizeof(Scalar) );
    if (Idy == NULL)
      Idy  = (Scalar*) malloc( N*M*sizeof(Scalar) );

    memset(Idx,   0, N*M*sizeof(Scalar));
    memset(Idy,   0, N*M*sizeof(Scalar));

    // no no no get some meaning inside
    for (int i = 0; i < N;i++)
      for (int j = 0; j < M;j++)
      {
        Idx[i*M + j] = i+1;// my indices count at 1, 0 is already outside (matlab style)
        Idy[i*M + j] = j+1;
      }
    bcIP.approximateImageGradients();
  }

  void finish()
  {
    if (Idx != NULL)
      free(Idx);
    Idx = NULL;
    if (Idy!= NULL)
      free(Idy);
    Idy = NULL;
  }

  Scalar* getIdx( ) 
  { 
      return Idx;
  }

  Scalar* getIdy( ) 
  { 
      return Idy;
  }

  int imgSize() {return  N*M;};

  void warp_noOmp_patchBased( P4i bigbox, M3 Hom ) 
  {
    bigbox[0] = max(0, bigbox[0]-boxRadius );
    bigbox[1] = max(0, bigbox[1]-boxRadius );
    bigbox[2] = min(N, bigbox[2]+boxRadius );
    bigbox[3] = min(M, bigbox[3]+boxRadius );

    for (int i = bigbox[0]; i < bigbox[2] ;i++)
    {
      P3 p = Hom * P3( i+1, bigbox[1], 1. );
      P3 q = Hom * P3(  0., 1, 0. );
      int j = bigbox[1];
      int pix = j + i*M;

      __m128d p1 = _mm_set_pd( p[0],p[0]-q[0] );
      __m128d p2 = _mm_set_pd( p[1],p[1]-q[1] );
      __m128d p3 = _mm_set_pd( p[2],p[2]-q[2] );

      __m128d q1 = _mm_set1_pd( q[0]*2 );
      __m128d q2 = _mm_set1_pd( q[1]*2 );
      __m128d q3 = _mm_set1_pd( q[2]*2 );

      for (; j < bigbox[3]-1 ;j+=2, pix+=2)
      {
        p1 = _mm_add_pd(p1,q1);
        p2 = _mm_add_pd(p2,q2);
        p3 = _mm_add_pd(p3,q3);

        _mm_storeu_pd( &Idx[pix],  _mm_div_pd( p1,p3 )) ;
        _mm_storeu_pd( &Idy[pix],  _mm_div_pd( p2,p3 )) ;
      }
      if( j == bigbox[3]-1 )
      {
#ifdef _DEBUG
        int noWay=0;
        if (pix>N*M-1 )
          noWay =1;
#endif

        P3 p = Hom * P3( i+1, bigbox[3], 1. );
        Idx[pix] = p[0]/p[2];
        Idy[pix] = p[1]/p[2];
      }
    }

    // could replace the bilinear warping by nn interpolation for more speed
//    bcIP.interpolate_box_short_fix( i_in_short, Idx, Idy, i_out_short, bigbox[0], bigbox[2], bigbox[1], bigbox[3]);
   bcIP.interpolate_box_short( i_in_short, Idx, Idy, i_out_short, bigbox[0], bigbox[2], bigbox[1], bigbox[3]);
  }


  void warp_noOmp_patchBased( std::vector<M3> Homs, int* segImg, int width, int height )
  {
    int nPixel = width*height;
    for( int i = 0; i<width;  i++ )
      for( int j = 0; j<height; j++ )
      {
        int pix = j + i*M;
        int currentSeg = segImg[pix];
        // speedup !: Hom*iK *p3
        P3 p = Homs[currentSeg] * P3( i+1, j+1, 1. );
        Idx[pix] = p[0]/p[2];
        Idy[pix] = p[1]/p[2];
      }

    bcIP.interpolate_short( i_in_short, nPixel, Idx, Idy, i_out_short );
  }

    /// check if warping should be done, original image is not to be warp: HERE the ids are the identity
  void warp_noOmp_patchBased( std::vector< Scalar >& idx, std::vector< Scalar >& idy, std::vector< std::pair<int,int> > pixelList )
  {
    for (int ij = 0; ij < pixelList.size(); ij++)
    {
      if (pixelList[ij].first <0 || pixelList[ij].first >= N*M|| idx.size() > N*M)
        i_inout_short[ ij ]  = 0; // could project to closest border pixel or mirror or ..
      else
        i_inout_short[ ij ]  = 128 * (short) (i_in_short[ pixelList[ij].first ]); //( img2[ pixelList[ij].first ]*Scalar(255.0) );
    }
    // also in the reference image case, I must compute a warp with the identiy map, since we are loking at a part of the image
    nPixel = idx.size();
    if (nPixel >0)
    {
      memcpy( Idx, &(idx[0]), sizeof(Scalar)*idx.size());
      memcpy( Idy, &(idy[0]), sizeof(Scalar)*idy.size());
//      bcIP.interpolate_noOmp_short( nPixel, i_in_short, Idx, Idy, i_out_short);
      bcIP.interpolate_short( i_in_short, nPixel, Idx, Idy, i_out_short);
    }
  };

  short* getWarpedShortImage()
  {
    if (i_out_short != NULL)
      return i_out_short;
    else
      return NULL;
  };

  short* getOrigShortImage()
  {
    if (i_inout_short != NULL)
      return i_inout_short;
    else
      return NULL;
  };


  int getNPixel(){return nPixel;};

  /////////////////////////////////////////////////
private:

  int nPixel;
  int nSegments;
  int nNormals;
  int N;
  int M;

  int dummyWarp;

  BiLinearIP<Scalar> bcIP;

  Scalar *Idx;
  Scalar *Idy;

  Scalar  * img2;	// new 
  Scalar *i_out;

  short* i_inout_short;
  short* i_in_short;
  short* i_out_short;
};

#endif // __GEN_HOM_HH
