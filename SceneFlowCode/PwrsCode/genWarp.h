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
License along with this software.*/

#ifndef _Gen_WARP_H
#define _Gen_WARP_H

#include "DataDefinitions.h"
#include "BilinearIp.h"
#include "genHom.h"

template<typename Scalar> class genWarp
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<int, 4>     P4i;
  
  genWarp( int N_, int M_, Scalar*& img2_, Scalar* p2d_)
    : nPixel(N_*M_), nSegments(0), nNormals(0), N(N_), M(M_), dummyWarp(0), bcIP(N_, M_, img2_), 
    locHom(NULL), Idx(NULL), Idy(NULL), p2d(p2d_), img2(img2_), i_out(NULL), 
    Idx2(NULL), Idy2(NULL), i_out2(NULL), warped1(false), warped2(false), validIndices(0),
    i_in_short(NULL), i_out_short(NULL)
  {init();generateShorts(img2_);};

  genWarp( int N_, int M_, Scalar*& img2_)
    : nPixel(N_*M_), nSegments(0), nNormals(0), N(N_), M(M_), dummyWarp(0), bcIP(N_, M_, img2_), 
    locHom(NULL), Idx(NULL), Idy(NULL), p2d(NULL), img2(img2_), i_out(NULL), 
    Idx2(NULL), Idy2(NULL), i_out2(NULL), warped1(false), warped2(false), validIndices(0), 
    i_in_short(NULL), i_out_short(NULL)
  {init();generateShorts(img2_);};

  ~genWarp(){finish();deleteShorts();};

  void generateShorts(Scalar*& img2_)
  {
    if (i_in_short  == NULL)
      i_in_short  = (short*) malloc( (N*M+M+1)*sizeof(short) );
//      i_in_short  = (short*) malloc( N*M*sizeof(short) );
    if (i_out_short == NULL)
      i_out_short = (short*) malloc( N*M*sizeof(short) );

    for (int ij = 0; ij < N*M;ij++)
    {
      i_in_short [ij]  = (short) ( img2_[ij]*Scalar(255.0) );
      i_out_short[ij]  = 128 * (short) ( img2_[ij]*Scalar(255.0) );
    }
    // duplicate out of bounds area
    for (int ij = (N-1)*M; ij < N*M;ij++)
      i_in_short [ij+M]  = (short) ( img2_[ij]*Scalar(255.0) );
    i_in_short [N*M+M] = (short) ( img2_[N*M]*Scalar(255.0) );
    //////////////////////////////////////
  }

  void deleteShorts()
  {
    if (i_in_short  != NULL)
       free( i_in_short );
    if (i_out_short != NULL)
       free( i_out_short );
    i_in_short   = NULL;
    i_out_short  = NULL;
    ///// 
  }

  void init()
  {
    if (Idx == NULL)
      Idx  = (Scalar*) malloc( N*M*sizeof(Scalar) );
    if (Idy == NULL)
      Idy  = (Scalar*) malloc( N*M*sizeof(Scalar) );
    if (i_out == NULL)
      i_out= (Scalar*) malloc( N*M*sizeof(Scalar) );

    memset(Idx,   0, N*M*sizeof(Scalar));
    memset(Idy,   0, N*M*sizeof(Scalar));
    memset(i_out, 0, N*M*sizeof(Scalar));

    // no no no get some meaning inside
    for (int i = 0; i < N;i++)
      for (int j = 0; j < M;j++)
      {
        Idx[i*M + j] = i+1;// my indices count at 1, 0 is already outside (matlab style)
        Idy[i*M + j] = j+1;
      }
    validIndices = 1;

    bcIP.approximateImageGradients();
  }

  void request2ndWarp()
  {
    if (Idx2 == NULL)
      Idx2  = (Scalar*) malloc( N*M*sizeof(Scalar) );
    if (Idy2 == NULL)
      Idy2  = (Scalar*) malloc( N*M*sizeof(Scalar) );
    if (i_out2 == NULL)
      i_out2= (Scalar*) malloc( N*M*sizeof(Scalar) );
    
    memset(Idx2,   0, N*M*sizeof(Scalar));
    memset(Idy2,   0, N*M*sizeof(Scalar));
    memset(i_out2, 0, N*M*sizeof(Scalar));

    // no no no get some meaning inside
    for (int i = 0; i < N;i++)
      for (int j = 0; j < M;j++)
      {
        Idx2[i*M + j] = i+1;// my indices count at 1, 0 is already outside (matlab style)
        Idy2[i*M + j] = j+1;
      }
  }

  // if the image remains constant (the reference image)
  void setDummy() {dummyWarp = 1;};

  void finish()
  {
    if (Idx != NULL)
      free(Idx);
    Idx = NULL;
    if (Idy!= NULL)
      free(Idy);
    Idy = NULL;
    if (i_out!= NULL)
      free(i_out);
    i_out = NULL;

    /// if there is a second pair census style
    if (Idx2 != NULL)
      free(Idx2);
    Idx2 = NULL;
    if (Idy2!= NULL)
      free(Idy2);
    Idy2 = NULL;
    if (i_out2!= NULL)
      free(i_out2);
    i_out2 = NULL;
  }

  Scalar* getWarp(bool orig = false)
  {
    if (p2d == NULL)
      return img2;
    else
      return i_out;
  }

  Scalar* getIdx( ) 
  { 
    if ( (p2d == NULL || locHom == NULL) & !validIndices)
      return NULL;
    else
      return Idx;
  }

  Scalar* getIdy( ) 
  { 
    if ( (p2d == NULL || locHom == NULL) & !validIndices)
      return NULL;
    else
      return Idy;
  }

  /// to generate the homographies from the 3d moving planes
  void setHom ( genHomography<Scalar>* genHom_)  
  {
    locHom = genHom_;
    int nNormals =  (int)(locHom->getSize());
    segmentHomos.resize(nNormals, M3(1,0,0,0,1,0,0,0,1));
  };

  /// check if warping should be done, original image is not to be warp: HERE the ids are the identity
  void warp_noOmp_patchBased( std::vector< Scalar >& idx, std::vector< Scalar >& idy )
  { 
    // alos in the reference image case, I must compute a ward with the identiy map, since we are loking at a part of the image
    nPixel = idx.size();
    memcpy( Idx, &(idx[0]), sizeof(Scalar)*idx.size());
    memcpy( Idy, &(idy[0]), sizeof(Scalar)*idy.size());

    validIndices = 1;

    bcIP.interpolate_noOmp_short( nPixel, i_in_short, Idx, Idy, i_out_short);

    warped1 = true;
  };

//  void warp_noOmp_patchBased( const mxArray*& Segments, const int proposal, const P4i bbox, const int* segImg )
  void warp_noOmp_patchBased( const mxArray*& Segments, const int proposal, const int bboxX1, const int bboxY1, const int bboxX2, const int bboxY2, const int* segImg )
  {
    if ( p2d == NULL || locHom == NULL )
      return;

    P4i bigbox( bboxX1, bboxY1, bboxX2, bboxY2 );
    bigbox[0] = max(0,   bigbox[0]-boxRadius );
    bigbox[1] = max(0,   bigbox[1]-boxRadius );
    bigbox[2] = min(N, bigbox[2]+boxRadius );
    bigbox[3] = min(M, bigbox[3]+boxRadius );

//    bigbox[2] = min(N-1, bigbox[2]+boxRadius );
//    bigbox[3] = min(M-1, bigbox[3]+boxRadius );

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
    nNormals =  (int) (locHom->getSize());

    segmentProcessed.assign(segmentHomos.size(), false );

    M3 Hom;

    int lastSeg = -1;
    for (int i = bigbox[0]; i < bigbox[2] ;i++)
      for (int j = bigbox[1]; j < bigbox[3] ;j++)
      {
        int pix = j + i*M;
        int seg = segImg[pix];

        if (lastSeg != seg )
        {
          lastSeg = seg;
          if (!segmentProcessed[seg])
          {
            Hom = locHom->getHom( proposal, seg );
            segmentProcessed[seg] = true;
            segmentHomos[seg] = Hom;
          }
          else
            Hom = segmentHomos[seg];
        }
        P3 p = Hom * P3( i+1, j+1, 1. );// a lot faster !!?? Hom*iK *p3
        Idx[pix] = p[0]/p[2];
        Idy[pix] = p[1]/p[2];
      }

    bcIP.interpolate_box_short( i_in_short, Idx, Idy, i_out_short, bigbox[0], bigbox[2], bigbox[1], bigbox[3]);
  }

  void warp_noOmp_patchBased(const mxArray*& Segments, const std::vector<int>& perSegIds)
  {
    if ( p2d == NULL || locHom == NULL )
      return;

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
    nNormals =  (int)(locHom->getSize());

    for (int i = 0; i < nSegments ;i++)
    {
      int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
      int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

      M3 Hom = locHom->getHom( perSegIds[i], i );
      for (int j=0;j<nIds;j++)
      {
        P3 p = Hom * P3( ids[j]/M+1, ids[j]%M+1, 1. );// a little faster !!?? Hom*iK *p3

        Idx[ids[j]] = p[0]/p[2];
        Idy[ids[j]] = p[1]/p[2];
      }
    }

    bcIP.interpolate_short( i_in_short, N*M, Idx, Idy, i_out_short );
  };


  void warp_noOmp_patchBased(const mxArray*& Segments, int nId)
  {
    if ( p2d == NULL || locHom == NULL )
      return;

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
    nNormals =  (int)(locHom->getSize());

    for (int i = 0; i < nSegments ;i++)
    {
      int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
      int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

      M3 Hom = locHom->getHom( nId, i );
      for (int j=0;j<nIds;j++)
      {
        P3 p = Hom * P3( ids[j]/M+1, ids[j]%M+1, 1. );
        Idx[ids[j]] = p[0]/p[2];
        Idy[ids[j]] = p[1]/p[2];
      }
    }

    bcIP.interpolate_noOmp( N*M, Idx, Idy, i_out);
  };

  short* getWarpedShortImage()
  {
    if (i_out_short != NULL)
      return i_out_short;
    else
      return NULL;
  };


  Scalar* getWarpedPartImage()
  {
    if (dummyWarp && !warped1)
      return img2;
    else
      return i_out;
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

  genHomography<Scalar>* locHom;

  std::vector<bool>  segmentProcessed;
  std::vector<M3>    segmentHomos;

  Scalar *Idx;
  Scalar *Idy;
  Scalar *p2d;

  Scalar *& img2;	
  Scalar *i_out;

  Scalar *Idx2;
  Scalar *Idy2;
  Scalar *i_out2;

  short* i_in_short;
  short* i_out_short;

  bool  warped1;
  bool  warped2;

  /// are the indices (Idx, Idy) valid indices - so can these be checked for being inside of the image, etc
  int validIndices;
};
#endif