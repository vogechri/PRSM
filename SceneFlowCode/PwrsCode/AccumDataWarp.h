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

////////////////////////////////////////////////////////////////
// Cmpute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#ifndef __ACCUM_DATA_WARP_H
#define __ACCUM_DATA_WARP_H
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define add_reg_pd( ma, out ) _mm_store_sd( &out, _mm_add_pd (ma, _mm_shuffle_pd(ma,ma, 1) ))

#define midC 160

#include "genWarp.h"
#include "OcclusionMappingBuffer.h" // maybe forward definition, include in cpp file
#include "DataDefinitions.h"

template<typename Scalar> class accumulateWarp
{
public:

  typedef Math::Mat3x3T< Scalar>     M3;
  typedef Math::VectorT< Scalar, 3>  P3;
  typedef Math::VectorT<int, 4>      P4i;

  accumulateWarp( int N_, int M_, Scalar thresh_, Scalar oobThresh_ = Scalar(0.1), int censusSize = 1 )
    : N(N_), M(M_), nSegments(0), nNormals(0), locWarp1(NULL), locWarp2(NULL), 
      maxDisp(N_),  maxMot(M_*N_), thresh(thresh_), oobThresh(oobThresh_)
  {
    setCensusSize( censusSize );
  };

  ~accumulateWarp(){};

  void setRefImage(Scalar* refImg_) {refImg = refImg_;};

  // should reflect the minimal depth possible - even better car: upper half should be 10m away or so
  void setMaxDisp(Scalar maxDisp_) {maxDisp = maxDisp_;};
  void setMaxMot (Scalar maxMot_)  {maxMot  = maxMot_;};

  void setWarp1 ( genWarp<Scalar>* genWarp_)  {locWarp1 = genWarp_;};
  void setWarp2 ( genWarp<Scalar>* genWarp_)  {locWarp2 = genWarp_;};

  std::vector<Scalar>& getScores()        {return locScores;};

  std::vector<Scalar>& getScores2()       {return locScores2;};

  std::vector<Scalar>& getFreeScores()    {return locFreeScores;};

  std::vector<int>&    getFreeVariables() {return freePixel;};
  ////////////////////////////////////////////////////////////

  void computeFullScoresCensus3OMP_Box(const mxArray*& Segments, int* segImg, 
    const int bboxX1, const int bboxY1, const int bboxX2, const int bboxY2, const std::vector<bool>& occlusionsR);

  void computeFullScoresCensus3OMP_short(const mxArray*& Segments, int* segImg, const std::vector<bool>& occlusionsR);

  /// ids == global ids are used to look up the oob pixels PatchIdGenerator<Scalar>* genPatchIds // std::vector< int >* otherEdgeIds
  void computeFullScoresMappedCensusNew_short( std::vector< int >& ids, std::vector< int >& lIds, int Mstep, const std::vector<bool>& oobPixR );

  void setCensusSize(int cSize)
  {
    int censusSize = ((2*cSize+1)*(2*cSize+1)-1)/2;

    dxC.clear();
    dyC.clear();
    dxyC.clear();
    dxyW.clear();
    dxC.resize(  censusSize );
    dyC.resize(  censusSize );
    dxyC.resize( censusSize  );
    dxyW.resize( censusSize  );

    int sumAll=0;
    for (int xx = -cSize; xx <= cSize; xx++ )
      for (int yy = -cSize; yy <= cSize; yy++ )
        if (sumAll < censusSize)
        {
          dxC[sumAll] = xx;
          dyC[sumAll] = yy;
          dxyC[sumAll] = N*xx+yy;
          dxyW[sumAll] = max(abs(xx), abs(yy))-1;
          sumAll++;
        }
  }

  ////////////////////////////////////////////////////////////
private:

  void censusBoxError4_short(const int* segImg, short *img1, short *img2, const int idsj, Scalar& locScore, int& miss );

  void censusBoxError4_short_SSE(const int* segImg, short *img1, short *img2, const int idsj, __m128i& scores );

  void censusBoxError4_short(const int* segImg, short *img1, short *img2, 
                             const int idsj, const int seg, const std::vector<bool>& occluded, 
                             std::vector<Scalar>& scores, std::vector<int>& hits, std::vector<int>& miss );

  /// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
  void censusBoxError4_short_SSE(const int* segImg, short *img1, short *img2, 
                                 const int idsj, const int seg, const std::vector<bool>& occluded, 
                                 std::vector<__m128i>& intscores, std::vector<int>& hits, std::vector<int>& miss );


  Scalar getPenaltyDisp( Scalar* Idx, Scalar* Idx2, bool& eval, std::vector<bool>& occluded, int i, int gId );

  Scalar getPenaltyMot( Scalar* Idx, Scalar* Idx2, Scalar* Idy, Scalar* Idy2, bool& eval, std::vector<bool>& occluded, int i, int gId );

  /// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
  void censusBinaryNew_short( std::vector<int>& rowPtr, std::vector<__m128i>& epsis, 
                              short *img0_a, short *img0_b, int gId,
                              std::vector<Scalar>& scores, int i, const int *dd, const int *ww, 
                              const int *dx, const int *dy, int steps, short Repsi[3], int Mstep, int Nstep);

  /// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
  void censusBinaryNew_shortSSE( std::vector<int>& rowPtr, std::vector<__m128i>& epsis, 
                                 short *img0_a, short *img0_b, int gId, std::vector<Scalar>& scores, 
                                 int i, short Repsi[3], Scalar dtaPen );


  char censusScore_short (short diff_a, short diff_b, short epsi = 160)
  {
    if (
      ( (diff_a > epsi) && (diff_b < epsi) ) ||
      ( (diff_a <-epsi) && (diff_b >-epsi) ) ||
      ( (diff_a < epsi) && (diff_b > epsi) ) ||
      ( (diff_a >-epsi) && (diff_b <-epsi) ) )
      return 1;
    return 0;
  }


  int N;
  int M;

  int nSegments;
  int nNormals;

  genWarp<Scalar>* locWarp1;
  genWarp<Scalar>* locWarp2;

  std::vector<Scalar> locScores;
  std::vector<Scalar> locFreeScores; 
  std::vector<int>    freePixel;

  std::vector<Scalar> locScores2;

  Scalar maxDisp;
  Scalar maxMot;

  Scalar thresh;
  Scalar oobThresh;
  
  /// to get the brightness at the original pixel in the image
  Scalar* refImg;

/// the neighbourhood in the censuss term
  std::vector<int> dxC;
  std::vector<int> dyC;
  std::vector<int> dxyC;
  std::vector<int> dxyW; // census epsilon implied for each of the pixel involved 
};
///////////////////////////////////

#if defined(INCLUDE_TEMPLATES) && !defined(__ACCUM_DATA_WARP_CPP_)
#include "AccumDataWarp.cpp"
#endif

#endif
