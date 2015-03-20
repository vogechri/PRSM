/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

//#define _lookupVersion_

#ifndef __ACCUM_DATA_WARP_CPP_
#define __ACCUM_DATA_WARP_CPP_

#include "AccumDataWarp.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define add_reg_pd( ma, out ) _mm_store_sd( &out, _mm_add_pd (ma, _mm_shuffle_pd(ma,ma, 1) ))

#define midC 160

#include "DataDefinitions.h"
#ifndef _NO_OPENMP
#include <omp.h>
#endif


template<class Scalar>
void
accumulateWarp<Scalar>::
computeFullScoresCensus3OMP_Box(const mxArray*& Segments, int* segImg, const int bboxX1, const int bboxY1, const int bboxX2, const int bboxY2, const std::vector<bool>& occlusionsR)
{
    P4i bigbox;
    bigbox[0] = max(0, bboxX1-boxRadius );
    bigbox[1] = max(0, bboxY1-boxRadius );
    bigbox[2] = min(N, bboxX2+boxRadius ); // if census8, no outer row needed
    bigbox[3] = min(M, bboxY2+boxRadius );

    assert (locWarp1 != NULL);
    assert (locWarp2 != NULL);

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

    locScores.assign( nSegments, 0);
    freePixel.assign( nSegments, 0);
    locFreeScores.assign( nSegments, 0);

    std::vector<bool> occluded(occlusionsR);

    short *img_s1 = locWarp1->getWarpedShortImage();
    short *img_s2 = locWarp2->getWarpedShortImage();
    std::vector<__m128i> intscores(nSegments, _mm_set1_epi16(0) );

    // more sense would be to use both images and check their coordinates
    Scalar* Idx = locWarp1->getIdx( );
    Scalar* Idy = locWarp1->getIdy( ); // NEW also its usage

    // both !
    Scalar* Idx2 = locWarp2->getIdx( );
    Scalar* Idy2 = locWarp2->getIdy( ); // NEW also its usage

    // all the penalties done
    std::vector<Scalar> scores(nSegments, 0);
    std::vector<int>    hits  (nSegments, 0);
    std::vector<int>    miss  (nSegments, 0);
    std::vector<int>    fails (nSegments, 0); // occluded != non occluded

    if (maxDisp != N) // involving the stereo case with extra check of geometry not too close
    {
      for (int i = bigbox[0]; i < bigbox[2] ;i++)
        for (int j = bigbox[1]; j < bigbox[3] ;j++)
        {
          int pix = j + i*M;

          int seg = segImg[pix];
          if ( ((maxDisp>=0) && (( Idx[pix] > (i +1) ) || ( Idx[pix] < (i +1 - maxDisp) ))) ||
               ((maxDisp< 0) && (( Idx[pix] < (i +1) ) || ( Idx[pix] > (i +1 - maxDisp) ))) ) // can not work with maxDisp < 0
          {
            // if in image big fail, if supposed to be but its not : fail
            if (Idx[pix] >= 0.5 && Idx[pix] <= N+0.5)
              locScores[seg] +=	Scalar(50.*thresh); // invalid
            else if (!occlusionsR[pix])
            {
              fails[seg]++;
            }
          }
          else
          {
            if (Idx[pix] >= 0.5 && Idx[pix] <= N+0.5)
            {
              if (occlusionsR[pix]) {fails[seg]++; occluded[pix] = true; continue;}
              freePixel[seg]++;

               if ( i < boxRadius || j < boxRadius || j >= M-boxRadius ) // j >= M-boxRadius || i >= N-boxRadius
                 censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, occluded, scores, hits, miss );
               else
                 censusBoxError4_short_SSE( segImg, img_s1, img_s2, pix, seg, occluded, intscores, hits, miss );
            }
            else
            {
              if (!occlusionsR[pix]) // those guys contribute to free score - but be aware better count good pixel per segment and assigment
              {
//                freePixel[i]++;
                fails[seg]++;
              }
              occluded[pix] = true;//basically used to not look up the neighbor for census
            }
          }
        }
    }
    else // for general motion (not only horizontal)
    {
      for (int i = bigbox[0]; i < bigbox[2] ;i++)
        for (int j = bigbox[1]; j < bigbox[3] ;j++)
        {
          int pix = j + i*M;

          if (j >= M || i >= N )
            continue;

          int seg = segImg[pix];

          if ( maxMot != N*M && (Idx[pix]-Idx2[pix]) * (Idx[pix]-Idx2[pix]) + (Idy[pix]-Idy2[pix]) * (Idy[pix]-Idy2[pix]) > maxMot*maxMot)
          {
            // if in image big fail, if supposed to be but its not : fail
            if ( (                ( Idx[pix] >= 0.5 && Idx[pix]  <= N+0.5 && Idy[pix]  >= 0.5 && Idy[pix]  <= M+0.5 )
               && ( Idx2==NULL || (Idx2[pix] >= 0.5 && Idx2[pix] <= N+0.5 && Idy2[pix] >= 0.5 && Idy2[pix] <= M+0.5)) )
               || !occlusionsR[pix] )
              locScores[seg] +=	Scalar(50.*thresh); // invalidZ
            continue;
          }

          if (                 ( Idx[pix] >= 0.5 && Idx[pix]  <= N+0.5 && Idy[pix]  >= 0.5 && Idy[pix]  <= M+0.5 )
            && ( Idx2==NULL || (Idx2[pix] >= 0.5 && Idx2[pix] <= N+0.5 && Idy2[pix] >= 0.5 && Idy2[pix] <= M+0.5)) )
          {
            if (occlusionsR[pix]) {fails[seg]++; occluded[pix] = true;continue;}
            freePixel[seg]++;


               if ( i < boxRadius || j < boxRadius || j >= M-boxRadius ) // j >= M-boxRadius || i >= N-boxRadius
                 censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, occluded, scores, hits, miss );
               else
                 censusBoxError4_short_SSE( segImg, img_s1, img_s2, pix, seg, occluded, intscores, hits, miss );
          }
          else
          {
            if (!occlusionsR[pix]) 
            {
//              freePixel[i]++;
              fails[seg]++;
            }
            occluded[pix] = true;//basically used to not look up the neighbor for census
          }
        }
      }
    
    Scalar locThresh = thresh*4./24.;

    //////////////////////////////////////////
    // adding up the stuff per segement involved:
    for (int i =0;i<nSegments;i++)
    {

      __m128i score = _mm_load_si128( &intscores[ i ] );
      // kill dumb highest 16 bits: never show up anyway so uncomment
//      score = _mm_and_si128( score, _mm_set_epi16 (0,  0xffff,  0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff) );

#ifdef __SSE4_1__
      scores[i] += locThresh * (score.m128i_i16[0] + score.m128i_i16[1] + score.m128i_i16[2] + score.m128i_i16[3] +
                                score.m128i_i16[4] + score.m128i_i16[5] + score.m128i_i16[6] + score.m128i_i16[7]);
#else
    // should be done externally and once - as long as there are less than 32768 pixel in a segment
      __m128i tempA = _mm_add_epi16( _mm_shufflelo_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ), _mm_shufflehi_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ) );
      score = _mm_add_epi16( tempA, _mm_unpackhi_epi16( tempA, tempA ) );
      tempA = _mm_add_epi16(  _mm_shufflelo_epi16( score, _MM_SHUFFLE(2,2,2,2) ), score );
      scores[i] += locThresh * _mm_extract_epi16( tempA, 0 );
#endif     

      if (scores[i] <= 0 && fails[i] <= 0) continue;
      
#ifdef _globalRescaling_
      Scalar scale = Scalar(hits[i]+miss[i]) / Scalar( max(1, hits[i]) );
#else
      Scalar scale = 1.;
#endif

#ifdef _data_truncation_
      locFreeScores[i] = min( freePixel[i] * oobThresh, scores[i]*scale);//thresh * 1./24.
#else
      locFreeScores[i] = scores[i]*scale;
#endif
      locScores[i]    += locFreeScores[i] + fails[i] * oobThresh;// fails are permenent as well as too large motion penalty
    }
    //////////////////////////////////////////
  };

template<class Scalar>
void
accumulateWarp<Scalar>::
censusBoxError4_short(const int* segImg, short *img1, short *img2, 
                              const int idsj, const int seg, const std::vector<bool>& occluded, 
                              std::vector<Scalar>& scores, std::vector<int>& hits, std::vector<int>& miss )
  {
    Scalar score (0);

    /// simplistic idea: just do something here
    short  i1 = img1[ idsj ];
    short  i2 = img2[ idsj ];
    short Repsi[3] = {1.25*128, 0.5*128, 0.};

    const int steps = 24;
    const int dd[24] = { -3-3*M, -2-3*M, -1-3*M, -3*M, 1-3*M, 2-3*M, 3-3*M, -3-2*M, -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, 3-2*M, -3-M, -2-M, -1-M, -M, 1-M, 2-M, 3-M, -3, -2, -1};
    const int ww[24] = {  3,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  2,  2,  3,  3,  2,  1,  1,  1,  2,  3,  3,  2,  1 };
    const int dy[24] = { -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1 };
    const int dx[24] = { -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1,  -1, 0,  0,  0 };
    Scalar scaleThresh = 4./24.;

    int lookId = ( refImg[idsj]*Scalar(255.0) );

    Scalar locThresh = thresh * scaleThresh;

    int localOobs = 0;
    int localHits = 0;

    int cy = idsj%M;

    for (int cid=0; cid < steps; cid++)
    {
      int nid = idsj + dd[cid];

      int py = cy + dy[cid];
      if ((nid>=0) && py >=0 && py < M )
      {
        if ( censusScore_short (img1[nid]-i1, img2[nid]-i2, Repsi[ww[cid]-1]) )
          score += locThresh;
      }
      else
        localOobs++;
    }
    score *= Scalar(steps) / std::max( Scalar(1.), Scalar(steps-localOobs) );
    scores[seg] += score;
    hits[seg]+=steps*2;
  }

/// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
template<class Scalar>
void
accumulateWarp<Scalar>::  
censusBoxError4_short_SSE(const int* segImg, short *img1, short *img2, 
                              const int idsj, const int seg, const std::vector<bool>& occluded, 
                              std::vector<__m128i>& intscores, std::vector<int>& hits, std::vector<int>& miss ) //std::vector<Scalar>& scores
  {
    /// simplistic idea: just do something here
    short  i1s = img1[ idsj ];
    short  i2s = img2[ idsj ];

#ifdef _lookupVersion_
    short Repsi[3] = {1.25*128, 0.5*128, 0.};
    int lookId  = (refImg[idsj]*Scalar(255.0));
#endif

    __m128i score = _mm_set1_epi16 (0);//,0,0,0, 0,0,0,0);
    static const __m128i ones  = _mm_set1_epi16 (1);//,1,1,1, 1,1,1,1);

    __m128i i1 = _mm_set1_epi16 (i1s);//,i1s,i1s,i1s, i1s,i1s,i1s,i1s);
    __m128i i2 = _mm_set1_epi16 (i2s);//,i2s,i2s,i2s, i2s,i2s,i2s,i2s);

    const int steps = 24;
    Scalar scaleThresh = 4./24.;
    const int dd[4] = {-3-3*M, -3-2*M, -3-M, -3};

#ifdef _lookupVersion_
//#define __Census__Compare__0__
#ifdef __Census__Compare__0__
    __m128i epsis[4] = {_mm_set_epi16 (32760, 0,0,0,0,0,0,0 ), _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0), 
      _mm_set_epi16 (32760, 0, 64, midC, midC, midC, 64, 0), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, midC, 64, 0)};
#else
    __m128i epsis[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
      _mm_set_epi16 (32760, 1, 64, midC, midC, midC, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, midC, 64, 1)};
#endif

#else
    static const __m128i epsis[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
      _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, 163, 64, 1) };

    static const __m128i iepsis[4] = {_mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 ), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
      _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -32760, -32760, -32760, -32760, -163, -64, -1) };
#endif

    const int runs = 4;
    Scalar locThresh = thresh * scaleThresh;

    for (int cid=0; cid < runs; cid++)
    {
      int nid  = idsj + dd[cid];

      // this is where the speedup is !
      __m128i i1_  = _mm_loadu_si128( (__m128i*) (&img1[nid]) );
      __m128i i2_  = _mm_loadu_si128( (__m128i*) (&img2[nid]) );

      __m128i  diff_i11 = _mm_sub_epi16(i1_, i1);
      __m128i  diff_i22 = _mm_sub_epi16(i2_, i2);
#ifdef _lookupVersion_
      __m128i mdiff_i11 = _mm_sub_epi16(i1, i1_);
      __m128i mdiff_i22 = _mm_sub_epi16(i2, i2_);
#endif

#ifndef __Census__Compare__0__
      // equality case is neglected: should a<b, a>b not >= in the equation, so hurts with 0 comparison so wither set to 1 in epsis or not like this
      __m128i a = _mm_cmpgt_epi16(diff_i11, epsis[cid]);
      __m128i b = _mm_cmplt_epi16(diff_i22, epsis[cid]);
#ifdef _lookupVersion_
      __m128i c = _mm_cmpgt_epi16(mdiff_i11, epsis[cid]);
      __m128i d = _mm_cmplt_epi16(mdiff_i22, epsis[cid]);
#else
      __m128i c = _mm_cmplt_epi16(diff_i11, iepsis[cid]);// not the opposite of above == case, why use diff and mdiff though
      __m128i d = _mm_cmpgt_epi16(diff_i22, iepsis[cid]);
#endif
      __m128i locScore =
      _mm_andnot_si128(
      _mm_and_si128(
        _mm_andnot_si128(_mm_and_si128(a,b), _mm_or_si128(a,b)),
        _mm_andnot_si128(_mm_and_si128(c,d), _mm_or_si128(c,d)) ),
         ones);
#else
      __m128i locScore = 
      _mm_and_si128( ones,
      _mm_or_si128(
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(diff_i11, epsis[cid]), _mm_cmplt_epi16(diff_i22, epsis[cid])), // i1_-i1>eps && i2_-i2<eps
        _mm_and_si128(_mm_cmplt_epi16(diff_i11, epsis[cid]), _mm_cmpgt_epi16(diff_i22, epsis[cid])) ),// i1_-i1<eps && i2_-i2>eps
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(mdiff_i11, epsis[cid]), _mm_cmplt_epi16(mdiff_i22, epsis[cid])),// -(i1_-i1)>eps && -(i2_-i2)<eps
        _mm_and_si128(_mm_cmplt_epi16(mdiff_i11, epsis[cid]), _mm_cmpgt_epi16(mdiff_i22, epsis[cid])) ) ));// -(i1_-i1)<eps && -(i2_-i2)>eps
#endif
      //( (diff_a > epsi) && (diff_b < epsi) ) ||
      //( (diff_a <-epsi) && (diff_b >-epsi) ) ||
      //( (diff_a < epsi) && (diff_b > epsi) ) ||
      //( (diff_a >-epsi) && (diff_b <-epsi) ) )

      score = _mm_add_epi16( locScore, score );
    }
    // all hits anyway
#ifndef _globalRescaling_

    __m128i tmp2 = _mm_load_si128( &intscores[seg] );
    _mm_store_si128( &intscores[seg], _mm_add_epi16( score, tmp2) );
    hits[seg]   += steps*2;
#else
    score = _mm_hadd_epi16(score, score);
    score = _mm_hadd_epi16(score, score);
    score = _mm_hadd_epi16(score, score);
    scores[seg] += locThresh * score.m128i_i16[0];

    hits[seg]   += steps*2;
#endif
  }

template<class Scalar>
void
accumulateWarp<Scalar>::  
computeFullScoresCensus3OMP_short(const mxArray*& Segments, int* segImg, const std::vector<bool>& occlusionsR)
{
    assert (locWarp1 != NULL);
    assert (locWarp2 != NULL);

    Scalar scaleData = 8. / ((2*boxRadius+1)*(2*boxRadius+1)-1);

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

    locScores.assign( nSegments, 0);
    freePixel.assign( nSegments, 0);
    locFreeScores.assign( nSegments, 0);
	
    short *img_s1 = locWarp1->getWarpedShortImage();
    short *img_s2 = locWarp2->getWarpedShortImage();

    // more sense would be to use both images and check their coordinates
    Scalar* Idx = locWarp1->getIdx( );
    Scalar* Idy = locWarp1->getIdy( ); // NEW also its usage

    // both !
    Scalar* Idx2 = locWarp2->getIdx( );
    Scalar* Idy2 = locWarp2->getIdy( ); // NEW also its usage

    if (maxDisp != N) // involving the stereo case with extra check of geometry not too close
    {
      for (int i = 0; i < nSegments ;i++)
      {
        int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
        int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

    		Scalar   locScore(0);
        __m128i mLocScore = _mm_set1_epi16(0);
        int fails(0);
        int miss(0); // out of segment/bounds

        for (int j=0;j<nIds;j++)
        {
          if ( ((maxDisp>=0) && (( Idx[ids[j]] > (ids[j]/M +1) ) || ( Idx[ids[j]] < (ids[j]/M +1 - maxDisp) ))) ||
               ((maxDisp< 0) && (( Idx[ids[j]] < (ids[j]/M +1) ) || ( Idx[ids[j]] > (ids[j]/M +1 - maxDisp) ))) )
          {
            if (Idx[ids[j]] >= 0.5 && Idx[ids[j]] <= N+0.5)
              locScores[i] +=	Scalar(50.*thresh); // invalid
            else if (!occlusionsR[ids[j]])
              fails++;
          }
          else
          {
            if (Idx[ids[j]] >= 0.5 && Idx[ids[j]] <= N+0.5)
            {
              if (occlusionsR[ids[j]]) {fails++; continue;}
              freePixel[i]++;
              if ( ids[j]/M < boxRadius || ids[j]%M < boxRadius || ids[j]%M >= M-boxRadius ) // j >= M-boxRadius || i >= N-boxRadius
                censusBoxError4_short( segImg, img_s1, img_s2, ids[j], locScore, miss);
              else
                censusBoxError4_short_SSE( segImg, img_s1, img_s2, ids[j], mLocScore); // ignore the not the segment case - as before
            }
            else
            {
              if (!occlusionsR[ids[j]]) // those guys contribute to free score - but be aware better count good pixel per segment and assigment
                fails++;
            }
          }
        }

#ifdef __SSE4_1__
        locScore += scaleData*thresh * (mLocScore.m128i_i16[0] + mLocScore.m128i_i16[1] + mLocScore.m128i_i16[2] + mLocScore.m128i_i16[3] + 
                                        mLocScore.m128i_i16[4] + mLocScore.m128i_i16[5] + mLocScore.m128i_i16[6] + mLocScore.m128i_i16[7] );
#else
    // should be done externally and once - as long as there are less than 32768 pixel in a segment
      __m128i tempA = _mm_add_epi16( _mm_shufflelo_epi16 ( mLocScore, _MM_SHUFFLE(0,1,2,3) ), _mm_shufflehi_epi16 ( mLocScore, _MM_SHUFFLE(0,1,2,3) ) );
      mLocScore = _mm_add_epi16( tempA, _mm_unpackhi_epi16( tempA, tempA ) );
      tempA = _mm_add_epi16(  _mm_shufflelo_epi16( mLocScore, _MM_SHUFFLE(2,2,2,2) ), mLocScore );
      locScore += scaleData*thresh * _mm_extract_epi16( tempA, 0 );
#endif
#ifdef _data_truncation_
        locScore = min ( freePixel[i] * oobThresh, locScore * Scalar(24*nIds)/Scalar(24*nIds-miss) );
#else
        locScore *=	Scalar(24*nIds)/Scalar(24*nIds-miss);
#endif
    		locFreeScores[i] += locScore;
        locScores[i]     += locScore + fails * oobThresh;// fails are permanent as well as too large motion penalty
      }
    }
    else // for general motion (not only horizontal)
    {
      for (int i = 0; i < nSegments ;i++)
      {
        int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
        int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

        Scalar locScore(0);
        __m128i mLocScore = _mm_set1_epi16(0);
        int fails(0);
        int miss(0); // out of segment/bounds
		
        // Idx starts at 1, not at 0
        for (int j=0;j<nIds;j++)
        {
          
          if ( maxMot != N*M && (Idx[ids[j]]-Idx2[ids[j]]) * (Idx[ids[j]]-Idx2[ids[j]]) + (Idy[ids[j]]-Idy2[ids[j]]) * (Idy[ids[j]]-Idy2[ids[j]]) > maxMot*maxMot)
          {
            // if in image big fail, if supposed to be but its not : fail
            if (                 ( Idx[ids[j]] >= 0.5 && Idx[ids[j]]  <= N+0.5 && Idy[ids[j]]  >= 0.5 && Idy[ids[j]]  <= M+0.5 )
              && ( Idx2==NULL || (Idx2[ids[j]] >= 0.5 && Idx2[ids[j]] <= N+0.5 && Idy2[ids[j]] >= 0.5 && Idy2[ids[j]] <= M+0.5)) )
              locScores[i] +=	Scalar(50.*thresh); // invalid
            else if (!occlusionsR[ids[j]])
              fails++;
            continue;
          }
          
          if (                 ( Idx[ids[j]] >= 0.5 && Idx[ids[j]]  <= N+0.5 && Idy[ids[j]]  >= 0.5 && Idy[ids[j]]  <= M+0.5 )
            && ( Idx2==NULL || (Idx2[ids[j]] >= 0.5 && Idx2[ids[j]] <= N+0.5 && Idy2[ids[j]] >= 0.5 && Idy2[ids[j]] <= M+0.5)) )
          {
            if (occlusionsR[ids[j]]) {fails++; continue;}//occluded[ids[j]] = true;
            freePixel[i]++;
            if ( ids[j]/M < boxRadius || ids[j]%M < boxRadius || ids[j]%M >= M-boxRadius ) // j >= M-boxRadius || i >= N-boxRadius
                censusBoxError4_short( segImg, img_s1, img_s2, ids[j], locScore, miss);
            else
              censusBoxError4_short_SSE( segImg, img_s1, img_s2, ids[j], mLocScore ); // ignore the not the segment case - as before
          }
          else
          {
            if (!occlusionsR[ids[j]]) 
              fails++;
          }
        }
#ifdef __SSE4_1__
        locScore += scaleData*thresh * (mLocScore.m128i_i16[0] + mLocScore.m128i_i16[1] + mLocScore.m128i_i16[2] + mLocScore.m128i_i16[3] + 
                                        mLocScore.m128i_i16[4] + mLocScore.m128i_i16[5] + mLocScore.m128i_i16[6] + mLocScore.m128i_i16[7] );
#else
    // should be done externally and once - as long as there are less than 32768 pixel in a segment
      __m128i tempA = _mm_add_epi16( _mm_shufflelo_epi16 ( mLocScore, _MM_SHUFFLE(0,1,2,3) ), _mm_shufflehi_epi16 ( mLocScore, _MM_SHUFFLE(0,1,2,3) ) );
      mLocScore = _mm_add_epi16( tempA, _mm_unpackhi_epi16( tempA, tempA ) );
      tempA = _mm_add_epi16(  _mm_shufflelo_epi16( mLocScore, _MM_SHUFFLE(2,2,2,2) ), mLocScore );
      locScore += scaleData*thresh * _mm_extract_epi16( tempA, 0 );
#endif
      // WHY NOT ?
#ifdef _data_truncation_
        locScore = min ( freePixel[i] * oobThresh, locScore * Scalar(24*nIds)/Scalar(24*nIds-miss) );
#else
        locScore *=	Scalar(24*nIds)/Scalar(24*nIds-miss);
#endif

        locFreeScores[i] = locScore;// + quickFixFails * oobThresh;
        locScores[i]    += locScore + fails * oobThresh;
      }
    }
  };// FIX IT WRT min


template<class Scalar>
void
accumulateWarp<Scalar>::
censusBoxError4_short(const int* segImg, short *img1, short *img2, 
                      const int idsj, Scalar& locScore, int& miss )
{
    Scalar score (0);

    /// simplistic idea: just do something here
    short  i1 = img1[ idsj ];
    short  i2 = img2[ idsj ];
    short Repsi[3] = {1.25*128, 0.5*128, 0.};

    const int steps = 24;
    const int dd[24] = { -3-3*M, -2-3*M, -1-3*M, -3*M, 1-3*M, 2-3*M, 3-3*M, -3-2*M, -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, 3-2*M, -3-M, -2-M, -1-M, -M, 1-M, 2-M, 3-M, -3, -2, -1};
    const int ww[24] = {  3,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  2,  2,  3,  3,  2,  1,  1,  1,  2,  3,  3,  2,  1 };
    const int dy[24] = { -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1 };
    const int dx[24] = { -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1,  -1, 0,  0,  0 };
    Scalar scaleThresh = 4./24.;

    int lookId = ( refImg[idsj]*Scalar(255.0) );

    Scalar locThresh = thresh * scaleThresh;

    int localOobs = 0;
    int localHits = 0;

    int cy = idsj%M;

    for (int cid=0; cid < steps; cid++)
    {
      int nid = idsj + dd[cid];

      int py = cy + dy[cid];
      if ((nid>=0) && py >=0 && py < M )
      {
        if ( censusScore_short (img1[nid]-i1, img2[nid]-i2, Repsi[ww[cid]-1]) )
          score += locThresh;
      }
      else
        miss++;
    }
    locScore += score;
  }

  /// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
template<class Scalar>
void
accumulateWarp<Scalar>::
censusBoxError4_short_SSE(const int* segImg, short *img1, short *img2, 
                          const int idsj, __m128i& scores )
  {
    /// simplistic idea: just do something here
    short  i1s = img1[ idsj ];
    short  i2s = img2[ idsj ];
#ifdef _lookupVersion_
    short Repsi[3] = {1.25*128, 0.5*128, 0.};
    int lookId  = (refImg[idsj]*Scalar(255.0));
#endif

    __m128i score = _mm_set1_epi16 (0);//,0,0,0, 0,0,0,0);
    static const __m128i ones  = _mm_set1_epi16 (1);//,1,1,1, 1,1,1,1);

    __m128i i1 = _mm_set1_epi16 (i1s);//,i1s,i1s,i1s, i1s,i1s,i1s,i1s);
    __m128i i2 = _mm_set1_epi16 (i2s);//,i2s,i2s,i2s, i2s,i2s,i2s,i2s);

    const int steps = 24;
    const int dd[4] = {-3-3*M, -3-2*M, -3-M, -3};

#ifdef __Census__Compare__0__
    __m128i epsis[4] = {_mm_set_epi16 (32760, 0,0,0,0,0,0,0 ), _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0), 
      _mm_set_epi16 (32760, 0, 64, midC, midC, midC, 64, 0), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, midC, 64, 0)};
#else
#ifdef _lookupVersion_
    __m128i epsis[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
      _mm_set_epi16 (32760, 1, 64, midC, midC, midC, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, midC, 64, 1)};
#else
    static const __m128i epsis[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
      _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, 163, 64, 1) };

    static const __m128i iepsis[4] = {_mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 ), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
      _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -32760, -32760, -32760, -32760, -163, -64, -1) };
#endif
#endif

    const int runs = 4;
    for (int cid=0; cid < runs; cid++)
    {
      int nid  = idsj + dd[cid];

      // this is where the speedup is !
      __m128i i1_  = _mm_loadu_si128( (__m128i*) (&img1[nid]) );
      __m128i i2_  = _mm_loadu_si128( (__m128i*) (&img2[nid]) );

      __m128i  diff_i11 = _mm_sub_epi16(i1_, i1);
      __m128i  diff_i22 = _mm_sub_epi16(i2_, i2);
#ifdef _lookupVersion_
      __m128i mdiff_i11 = _mm_sub_epi16(i1, i1_);
      __m128i mdiff_i22 = _mm_sub_epi16(i2, i2_);
#endif

#ifndef __Census__Compare__0__
      // equality case is neglected: should a<b, a>b not >= in the equation, so hurts with 0 comparison so wither set to 1 in epsis or not like this
      __m128i a = _mm_cmpgt_epi16(diff_i11, epsis[cid]);
      __m128i b = _mm_cmplt_epi16(diff_i22, epsis[cid]);
#ifdef _lookupVersion_
      __m128i c = _mm_cmpgt_epi16(mdiff_i11, epsis[cid]);
      __m128i d = _mm_cmplt_epi16(mdiff_i22, epsis[cid]);
#else
      __m128i c = _mm_cmplt_epi16(diff_i11, iepsis[cid]);// not the opposite of above == case, why use diff and mdiff though
      __m128i d = _mm_cmpgt_epi16(diff_i22, iepsis[cid]);
#endif
      __m128i locScore =
      _mm_andnot_si128(
      _mm_and_si128(
        _mm_andnot_si128(_mm_and_si128(a,b), _mm_or_si128(a,b)),
        _mm_andnot_si128(_mm_and_si128(c,d), _mm_or_si128(c,d)) ),
         ones);
#else
      __m128i locScore = 
      _mm_and_si128( ones,
      _mm_or_si128(
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(diff_i11, epsis[cid]), _mm_cmplt_epi16(diff_i22, epsis[cid])), // i1_-i1>eps && i2_-i2<eps
        _mm_and_si128(_mm_cmplt_epi16(diff_i11, epsis[cid]), _mm_cmpgt_epi16(diff_i22, epsis[cid])) ),// i1_-i1<eps && i2_-i2>eps
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(mdiff_i11, epsis[cid]), _mm_cmplt_epi16(mdiff_i22, epsis[cid])),// -(i1_-i1)>eps && -(i2_-i2)<eps
        _mm_and_si128(_mm_cmplt_epi16(mdiff_i11, epsis[cid]), _mm_cmpgt_epi16(mdiff_i22, epsis[cid])) ) ));// -(i1_-i1)<eps && -(i2_-i2)>eps
#endif

      score = _mm_add_epi16( locScore, score );
    }

     scores = _mm_add_epi16( scores, score );
  };

template<class Scalar>
Scalar
accumulateWarp<Scalar>::
getPenaltyDisp( Scalar* Idx, Scalar* Idx2, bool& eval, std::vector<bool>& occluded, int i, int gId )
{
    Scalar penalty(0);eval = true;

    if ( ((maxDisp>=0) && (( Idx[i]>Idx2[i] ) || ( Idx[i] < Idx2[i] - maxDisp))) ||
         ((maxDisp< 0) && (( Idx[i]<Idx2[i] ) || ( Idx[i] > Idx2[i] - maxDisp))) ) // can not work with maxDisp < 0
    {
      if (Idx[i] >= 0.5 && Idx[i] <= N+0.5)
        penalty = Scalar(50.)*thresh; // invalid since in image
      else if (!occluded[ gId ] )     // not in image but should be
        penalty = oobThresh; // wrong oob decision
      eval = false;
    }
    else // valid disparity
    {
      if (Idx[i] >= 0.5 && Idx[i] <= N+0.5) // inside
      {
        if( occluded[ gId] ) // should be visible
        {
          penalty = oobThresh;
          eval = true;
          occluded[ gId ] = false;
        }
      }
      else
      {
        if ( !occluded[ gId ] ) // should be occluded
        {
          occluded[ gId ] = true;
          penalty = oobThresh;
        }
        eval = false;
      }
    }
    return penalty;
  }


template<class Scalar>
Scalar
accumulateWarp<Scalar>::
getPenaltyMot( Scalar* Idx, Scalar* Idx2, Scalar* Idy, Scalar* Idy2, bool& eval, std::vector<bool>& occluded, int i, int gId ) 
{
    Scalar penalty(0);eval = true;

    if ( maxMot != N*M && (Idx[i]-Idx2[i]) * (Idx[i]-Idx2[i]) + (Idy[i]-Idy2[i]) * (Idy[i]-Idy2[i]) > maxMot*maxMot)
    {
      if ( ( Idx[i]  >= 0.5 && Idx[i]  <= N+0.5 && Idy[i]  >= 0.5 && Idy[i]  <= M+0.5 ) &&
           ( Idx2[i] >= 0.5 && Idx2[i] <= N+0.5 && Idy2[i] >= 0.5 && Idy2[i] <= M+0.5 ) ) // Idx2 ==NULL, never
           penalty = Scalar(50.)*thresh; // invalid
      else if (!occluded[ gId] )
        penalty = oobThresh; // wrong oob decision
      eval = false;
    }
    else
    {
      if ( ( Idx[i]  >= 0.5 && Idx[i]  <= N+0.5 && Idy[i]  >= 0.5 && Idy[i]  <= M+0.5 ) &&
           ( Idx2[i] >= 0.5 && Idx2[i] <= N+0.5 && Idy2[i] >= 0.5 && Idy2[i] <= M+0.5 ) ) // Idx2 ==NULL, never
      {
        if( occluded[ gId ] ) // inside but supposed to be not
        {
          eval = false;
          penalty = oobThresh; 
          // new:
          occluded[ gId ] = false;
          eval = true;
        }
      }
      else
      {
        if ( !occluded[ gId ] ) // should be occluded
        {
          occluded[ gId ] = true;
          penalty = oobThresh;
        }
        eval = false;
      }
    }
    return penalty;
  }

template<class Scalar>
void
accumulateWarp<Scalar>::
computeFullScoresMappedCensusNew_short( std::vector< int >& ids, std::vector< int >& lIds, int Mstep, const std::vector<bool>& oobPixR )
  {
    assert (locWarp1 != NULL);
    assert (locWarp2 != NULL);

    const int csize = (sqrt((double)(dxyC.size()*2+1))-1)/2;
    int nPixel = locWarp1->getNPixel();
    int Nstep = nPixel / Mstep;

    std::vector<int> rowPtr(csize*2+1);
    for(int i=0;i<csize*2+1;i++)
      rowPtr[i] = (i-csize)*Nstep - csize;
    ///////////////

    const int steps1 = 8;
    const int dd1[8] = { -1-M, -M, 1-M, -1, 1, -1+M, M, 1+M};
    const int dy1[8] = { -1,  0,  1, -1, 1, -1, 0, 1};
    const int dx1[8] = { -1, -1, -1,  0, 0, 1, 1, 1};
    const int ww1[8] = { 1, 1, 1, 1, 1, 1, 1, 1};

    const int steps3 = 48;
    const int dd3[48] = { -3-3*M, -2-3*M, -1-3*M, -3*M, 1-3*M, 2-3*M, 3-3*M, -3-2*M, -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, 3-2*M, -3-M, -2-M, -1-M, -M, 1-M, 2-M, 3-M, -3, -2, -1, 1, 2, 3, -3+M, -2+M, -1+M, M, 1+M, 2+M, 3+M, -3+2*M, -2+2*M, -1+2*M, 2*M, 1+2*M, 2+2*M, 3+2*M, -3+3*M, -2+3*M, -1+3*M, 3*M, 1+3*M, 2+3*M, 3+3*M};
    const int ww3[48] = {  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, 2, 1, 1, 1, 2, 3,  3,  2, 1, 1,2,3,3,2,1,1,1,2,3,3,2,2,2,2,2,3,3,3,3,3,3,3,3};
    const int dy3[48] = { -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  1,  2,  3, -3, -2, -1,  0,  1,  2, 3, -3,-2, -1,  0,  1,  2, 3, -3, -2, -1,  0,  1,  2,  3 };
    const int dx3[48] = { -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1,  -1, 0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1, 1,  2, 2,  2,  2,  2,  2, 2,  3,  3,  3,  3,  3,  3,  3 };
//    Scalar scaleThresh = 4./48.;

    const int steps2 = 24;
    const int dd2[24] = { -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, -2-M, -1-M, -M, 1-M, 2-M, -2, -1, 1, 2, -2+M, -1+M, M, 1+M, 2+M, -2+2*M, -1+2*M, 2*M, 1+2*M, 2+2*M };
    const int ww2[24] = { 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2};
    const int dy2[24] = { -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1, 1, 2, -2, -1,  0,  1,  2, -2, -1,  0,  1,  2 };
    const int dx2[24] = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
//    Scalar scaleThresh = 4./24.;
    __m128i mm_Mask3;
    __m128i mm_Mask2;

    short Repsi[3] = {1.25*128, 0.5*128, 0.};
    int cSize= dxyC.size();
    Scalar dtaPen (thresh*2./Scalar(cSize));

    const int *dd, *ww, *dx, *dy;
    int steps=8;

    std::vector<__m128i> epsis(csize*2+1, _mm_set1_epi32 (0));
    // epsilons:
    switch(csize)
    {
    case 1: 
      // idea here: load per hand, only one run
      // setting 160 is about 1.25: 1.25*255*128 = 160
//      epsis[0] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,160,160,160 );
      epsis[0] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,midC,midC,midC );
      epsis[2] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,midC,midC,midC );
      epsis[1] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,midC,32760,midC );
      mm_Mask3 = _mm_set_epi16 (0, 0, 0, 0, 0, 0xffff, 0xffff, 0xffff);
      mm_Mask2 = _mm_set_epi16 (0, 0, 0, 0, 0, 0xffff, 0, 0xffff);
      dd = dd1;
      ww = ww1;
      dx = dx1;
      dy = dy1;
      break;
    default:
    case 2:
      // idea load all 4, last per hand -> distribute, 3 runs only
      epsis[0] = _mm_set_epi16 (32760, 32760, 32760,  64,  64,  64,  64, 64 );
      epsis[1] = _mm_set_epi16 (32760, 32760, 32760,  64,   midC,  midC,   midC, 64);
      epsis[2] = _mm_set_epi16 (32760, 32760, 32760,  64,   midC, 32760, midC, 64);
      epsis[3] = _mm_set_epi16 (32760, 32760, 32760,  64,   midC,   midC,   midC, 64);
      epsis[4] = _mm_set_epi16 (32760, 32760, 32760,  64,  64,  64,  64, 64 );
      mm_Mask3 = _mm_set_epi16 (0, 0, 0, 0, 0xffff, 0xffff, 0xffff, 0);
      mm_Mask2 = _mm_set_epi16 (0, 0, 0, 0, 0xffff, 0, 0xffff, 0);
      dd = dd2;
      ww = ww2;
      dx = dx2;
      dy = dy2;
      steps= 24;
      break;
    case 3: 
      // could distribute last on all ? or not ?
      epsis[0] = _mm_set_epi16 (32760, 0,0,0,0,0,0,0 );
      epsis[1] = _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0);
      epsis[2] = _mm_set_epi16 (32760, 0, 64, midC, midC, midC, 64, 0);
      epsis[3] = _mm_set_epi16 (32760, 0, 64, midC, 32760, midC, 64, 0);
      epsis[6] = _mm_set_epi16 (32760, 0,0,0,0,0,0,0 );
      epsis[5] = _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0);
      epsis[4] = _mm_set_epi16 (32760, 0, 64, midC, midC, midC, 64, 0);
      mm_Mask3 = _mm_set_epi16 (0, 0, 0, 0xffff, 0xffff, 0xffff, 0, 0);
      mm_Mask2 = _mm_set_epi16 (0, 0, 0, 0xffff, 0, 0xffff, 0, 0);
      dd = dd3;
      ww = ww3;
      dx = dx3;
      dy = dy3;
      steps= 48;
      break;
    }

    short *img_s1 = locWarp1->getWarpedShortImage();
    short *img_s2 = locWarp2->getWarpedShortImage();

    /////////////////////////// 

    // accumulated scores
    locScores.clear();
    locScores.resize(nPixel,0);

    // penalties
    locScores2.clear();
    locScores2.resize(nPixel,0);

    // NEW also its usage
    Scalar *img0_1 = locWarp1->getWarpedPartImage();
    Scalar *img0_2 = locWarp2->getWarpedPartImage();

    Scalar* Idx0   = locWarp1->getIdx( );
    Scalar* Idy0   = locWarp1->getIdy( );
    Scalar* Idx0_2 = locWarp2->getIdx( );
    Scalar* Idy0_2 = locWarp2->getIdy( );

    std::vector<bool> occluded0(oobPixR);

    assert ( Idx0_2 != NULL);

    Scalar locScore0(0), locScore1(0);
    int fails0(0), fails1(0);

    if (maxDisp != N)// case of bounded motion checking
    {
      for (int i = 0; i < nPixel ;i++)
      {
        if (lIds[i] < 0) continue;

        bool eval0 = true;      
        locScores2[i]  = getPenaltyDisp( Idx0, Idx0_2, eval0, occluded0, i, ids[i] ) ;
        if ( occluded0[ids[i]] ) continue;

        int lookId( refImg[ids[i]]*Scalar(255.0) );

        if (i/Nstep < csize || i/Nstep >= Mstep-csize || i%Nstep < csize || i%Nstep >= Nstep-csize )
           censusBinaryNew_short( rowPtr, epsis, img_s1, img_s2, ids[i], locScores, i, dd, ww, dx, dy, steps, Repsi, Nstep, Mstep);//, true  ) ;
        else
           censusBinaryNew_shortSSE( rowPtr, epsis, img_s1, img_s2, ids[i], locScores, i, Repsi, dtaPen ) ;//Nstep, Mstep, steps, 
      }
    }
    else // for general motion (not only horizontal)
    {
      for (int i = 0; i < nPixel ;i++)
      {
        if (lIds[i] < 0) continue;

        bool eval0 = true;
        locScores2[i]  = getPenaltyMot( Idx0, Idx0_2, Idy0, Idy0_2, eval0, occluded0, i, ids[i] ) ;
        // same for second segment at pixel

        if ( occluded0[ids[i]] ) continue;

        int lookId( refImg[ids[i]]*Scalar(255.0) );

        if (i/Nstep < csize || i/Nstep >= Mstep-csize || i%Nstep < csize || i%Nstep >= Nstep-csize )
          censusBinaryNew_short( rowPtr, epsis, img_s1, img_s2, ids[i], locScores, i, dd, ww, dx, dy, steps, Repsi, Nstep, Mstep);//, true ) ;
        else
           //censusBinaryNew_short( rowPtr, epsis, img_s1, img_s2, ids[i], locScores, i, dd, ww, dx, dy, steps, Repsi, Nstep, Mstep, false ) ;
           censusBinaryNew_shortSSE( rowPtr, epsis, img_s1, img_s2, ids[i], locScores, i, Repsi, dtaPen );//Nstep, Mstep, steps

      };
    }
  }


template<class Scalar>
void
accumulateWarp<Scalar>::
censusBinaryNew_short( std::vector<int>& rowPtr, std::vector<__m128i>& epsis, 
                       short *img0_a, short *img0_b, int gId,
                       std::vector<Scalar>& scores, int i, const int *dd, const int *ww, 
                       const int *dx, const int *dy, int steps, short Repsi[3], int Mstep, int Nstep)
  {
    // i : local position
    short  i0_a = img0_a[ i ];
    short  i0_b = img0_b[ i ];

    // dependent on gray value in center:
    int score=0;
    int hits= 0;

    int cSize= dxyC.size();
    Scalar dtaPen (thresh*2./Scalar(cSize));

    // somehow that worked, with check see below
    int yi = i/Mstep;
    int xi = i%Mstep;

    for ( int j=0; j<steps; j++ )
    {
      // somehow that worked above, with this check
      if (yi+dy[j] >= 0 && yi+dy[j] < Nstep && xi+dx[j] >= 0 && xi+dx[j] < Mstep)
      {
        int nid = i+dy[j]*Mstep+dx[j];
        if ( censusScore_short (img0_a[nid]-i0_a, img0_b[nid]-i0_b, Repsi[ww[j]-1]) )
          score++;
        hits++;
      }
    }

#ifdef _data_truncation_
    scores[i] = min( oobThresh, score * steps * (dtaPen / std::max( Scalar(1.), Scalar(hits) )));
#else
    scores[i] = score * steps * (dtaPen / std::max( Scalar(1.), Scalar(hits) ));
#endif
    return;
  }

template<class Scalar>
void
accumulateWarp<Scalar>::
censusBinaryNew_shortSSE( std::vector<int>& rowPtr, std::vector<__m128i>& epsis, 
                          short *img0_a, short *img0_b, int gId, std::vector<Scalar>& scores, 
                          int i, short Repsi[3], Scalar dtaPen )
  {
    __m128i score = _mm_set1_epi16 (0);//,0,0,0, 0,0,0,0);
    __m128i ones  = _mm_set1_epi16 (1);//,1,1,1, 1,1,1,1);

    __m128i i1 = _mm_set1_epi16 (img0_a[ i ]);
    __m128i i2 = _mm_set1_epi16 (img0_b[ i ]);

    for(int k =0;k< rowPtr.size(); k++)
    {
      int nid = i+rowPtr[k];
      __m128i i1_ = _mm_loadu_si128( (__m128i*) ( &img0_a[ nid ] ) );
      __m128i i2_ = _mm_loadu_si128( (__m128i*) ( &img0_b[ nid ] ) );

      // ignore epsi adjustment:
      __m128i  diff_i11 = _mm_sub_epi16(i1_, i1);
      __m128i  diff_i22 = _mm_sub_epi16(i2_, i2);
      __m128i mdiff_i11 = _mm_sub_epi16(i1, i1_);
      __m128i mdiff_i22 = _mm_sub_epi16(i2, i2_);

      __m128i locScore = 
      _mm_and_si128( ones,
      _mm_or_si128(
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(diff_i11, epsis[k]), _mm_cmplt_epi16(diff_i22, epsis[k])), // i1_-i1>eps && i2_-i2<eps
        _mm_and_si128(_mm_cmplt_epi16(diff_i11, epsis[k]), _mm_cmpgt_epi16(diff_i22, epsis[k])) ),// i1_-i1<eps && i2_-i2>eps
      _mm_or_si128(
        _mm_and_si128(_mm_cmpgt_epi16(mdiff_i11, epsis[k]), _mm_cmplt_epi16(mdiff_i22, epsis[k])),// -(i1_-i1)>eps && -(i2_-i2)<eps
        _mm_and_si128(_mm_cmplt_epi16(mdiff_i11, epsis[k]), _mm_cmpgt_epi16(mdiff_i22, epsis[k])) ) ));// -(i1_-i1)<eps && -(i2_-i2)>eps

      score = _mm_add_epi16( locScore, score );
    }

#ifdef __SSE4_1__
    scores[i] = dtaPen* (score.m128i_i16[0] + score.m128i_i16[1] + score.m128i_i16[2] + score.m128i_i16[3] + 
                         score.m128i_i16[4] + score.m128i_i16[5] + score.m128i_i16[6] + score.m128i_i16[7]);
#else
    // should be done externally and once - as long as there are less than 32768 pixel in a segment
      __m128i tempA = _mm_add_epi16( _mm_shufflelo_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ), _mm_shufflehi_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ) );
      score = _mm_add_epi16( tempA, _mm_unpackhi_epi16( tempA, tempA ) );
      tempA = _mm_add_epi16(  _mm_shufflelo_epi16( score, _MM_SHUFFLE(2,2,2,2) ), score );
      scores[i] += dtaPen * _mm_extract_epi16( tempA, 0 );
#endif
#ifdef _data_truncation_
      scores[i] = min( scores[i], oobThresh );
#endif
    return;
  }
#endif
