/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich

This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Christoph Vogel).

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __genScore__cpp
#define __genScore__cpp

//////////////////////////////////////////////////////////////////////////////
// evaluate data cost //
//////////////////////////////////////////////////////////////////////////////
// uses half neigh only, to speed up evaluation
#define OriginalVersion

#include "DataDefinitionsVC.h"

using namespace std;
using namespace Math;

// census which is slower but more official then the implemented version - in the results virtually no difference at all
//#define __Census__Compare__0__
// worked on all computers, but just in case:
//#define _nosse_

#include "DataDefinitionsVC.h"

#include "genScore.h"
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;
///////////////////////////////////////

template<class Scalar>
void genScore<Scalar>::
  censusBoxError4_short_SSE_inner( const int* segImg, 
  const short *img1, const short *img2, std::vector<__m128i>& intscores, 
  P4i innerbox, const std::vector<int>& vSegs, 
  const Scalar* const Idx , const Scalar* const Idy)
{

  const int steps = 24;
  Scalar scaleThresh = 4./24.;
  const int dd[4] = {-3-3*M, -3-2*M, -3-M, -3};
  static const __m128i ones  = _mm_set1_epi16 (1);
  const int runs = 4;
  Scalar locThresh = thresh * scaleThresh;

  static const __m128i epsisf[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
    _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, 163, 64, 1) };

  static const __m128i iepsisf[4] = {_mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 ), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
    _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -32760, -32760, -32760, -32760, -163, -64, -1) };
  ///
  for (int i = innerbox[0]; i < innerbox[2] ;i++)
    for (int j = innerbox[1]; j < innerbox[3] ;j++)
    {
      int idsj = j + i*M;
      int seg = segImg[idsj];

      // skip if not of interest:
      if (vSegs[seg]<0) continue;

      if ( ( Idx[idsj] < 0.5 || Idx[idsj] > N+0.5 || Idy[idsj] < 0.5 || Idy[idsj] > M+0.5 )) continue;
//      if (Idx[idsj] < 0.5 && Idx[idsj] > N+0.5) continue;

      /// simplistic idea: just do something here
      short  i1s = img1[ idsj ];
      short  i2s = img2[ idsj ];

      __m128i score = _mm_set1_epi16 (0);

      __m128i i1 = _mm_set1_epi16 (i1s);
      __m128i i2 = _mm_set1_epi16 (i2s);

      for (int cid=0; cid < runs; cid++)
      {
        int nid  = idsj + dd[cid];

        // this is where the speedup is !
        __m128i i1_  = _mm_loadu_si128( (__m128i*) (&img1[nid]) );
        __m128i i2_  = _mm_loadu_si128( (__m128i*) (&img2[nid]) );

        __m128i  diff_i11 = _mm_sub_epi16(i1_, i1);
        __m128i  diff_i22 = _mm_sub_epi16(i2_, i2);

        // equality case is neglected: should a<b, a>b not >= in the equation, so hurts with 0 comparison so wither set to 1 in epsis or not like this
        __m128i a = _mm_cmpgt_epi16(diff_i11, epsisf[cid]);
        __m128i b = _mm_cmplt_epi16(diff_i22, epsisf[cid]);

        __m128i c = _mm_cmplt_epi16(diff_i11, iepsisf[cid]);// not the opposite of above == case, why use diff and mdiff though
        __m128i d = _mm_cmpgt_epi16(diff_i22, iepsisf[cid]);

        __m128i locScore =
          _mm_andnot_si128(
          _mm_and_si128(
          _mm_andnot_si128(_mm_and_si128(a,b), _mm_or_si128(a,b)),
          _mm_andnot_si128(_mm_and_si128(c,d), _mm_or_si128(c,d)) ),
          ones);

        score = _mm_add_epi16( locScore, score );
      }
      // all hits anyway
/*
      std::vector<__m128i> intscoresD(nSegments, _mm_set1_epi16(0) );
      std::vector<int>    hits  (nSegments, 0);
      std::vector<int>    miss  (nSegments, 0);
      censusBoxError4_short_SSE( segImg, img1, img2, idsj, seg, intscoresD, hits, miss );//occluded,
      __m128i testScore = intscoresD[seg];

      __m128i dscore = _mm_sub_epi16( testScore, score );

      int what;
      if (dscore.m128i_i64[1] != 0 || dscore.m128i_i64[0] != 0)
        what=0;
*/
      // per pix storage
      _mm_store_si128( &intscores[idsj], score );
//      __m128i tmp2 = _mm_load_si128( &intscores[seg] );
//      _mm_store_si128( &intscores[seg], _mm_add_epi16( score, tmp2) );
    }
}


/// does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled: p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
template<class Scalar>
void genScore<Scalar>::
  computeFullScoresCensus3OMP_Box( short *img_s1, short *img_s2, const Scalar* const Idx , const Scalar* const Idy,
  int nSegments_, int* segImg, const int bboxX1, const int bboxY1, const int bboxX2, const int bboxY2, 
  const std::vector<int>& vSegs, const std::vector<bool>& occlusionsR )
{
  const P4i bigbox( max(0, bboxX1-boxRadius ), max(0, bboxY1-boxRadius ), min(N, bboxX2+boxRadius ) ,min(M, bboxY2+boxRadius ));
  nSegments = nSegments_;

  //const P4i innerbox( max(3, bboxX1-boxRadius ), max(3, bboxY1-boxRadius ), min(N-3, bboxX2+boxRadius ) ,min(M-3, bboxY2+boxRadius ));
  //std::vector<__m128i> intscores_inner(nSegments, _mm_set1_epi16(0) );
  //std::vector<__m128i> intscores_pix(N*M, _mm_set1_epi16(0) );
  //censusBoxError4_short_SSE_inner( segImg, img_s1, img_s2, intscores_pix, innerbox, vSegs, Idx , Idy );


  locScores.assign( nSegments, 0);
  freePixel.assign( nSegments, 0);
  locFreeScores.assign( nSegments, 0);

  std::vector<__m128i> intscores(nSegments, _mm_set1_epi16(0) );

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

        // skip if not of interest:
        if (vSegs[seg]<0) continue;
#ifdef __2dboundControl__
        if ( ((maxDisp>=0) && (( Idx[pix] > (i +1) ) || ( Idx[pix] < (i +1 - maxDisp) ))) ||
          ((maxDisp< 0) && (( Idx[pix] < (i +1) ) || ( Idx[pix] > (i +1 - maxDisp) ))) ) // can now work with maxDisp < 0
        {
          locScores[seg] +=	Scalar(50.*thresh); // invalid,
        }
        else
#endif
        {
          if (Idx[pix] >= 0.5 && Idx[pix] <= N+0.5) // in
          {
            if (occlusionsR[pix]) {fails[seg]++; continue;};
            freePixel[seg]++;

            if ( i < boxRadius || j < boxRadius || j >= M-boxRadius || i >= N-boxRadius )
              censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, scores, hits, miss );//, occluded
            else
#ifndef _nosse_
              censusBoxError4_short_SSE( segImg, img_s1, img_s2, pix, seg, intscores, hits, miss );//occluded,
//             _mm_store_si128( &intscores[seg], _mm_add_epi16( _mm_load_si128( &intscores_pix[pix] ), _mm_load_si128( &intscores[seg] )) );
#else
              censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, occluded, scores, hits, miss );
#endif
          }
          else // out
          {
            if (!occlusionsR[pix]) // those guys contribute to free score - but be aware better count good pixel per segment and assigment
            {
              fails[seg]++;
            }
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

        int seg = segImg[pix];
        // skip if not of interest:
        if (vSegs[seg]<0) continue;

#ifdef __2dboundControl__
        if ( maxMot != N*M && (Idx[pix]-i-1) * (Idx[pix]-i-1) + (Idy[pix]-j-1) * (Idy[pix]-j-1) > maxMot*maxMot)
        {
          locScores[seg] +=	Scalar(50.*thresh); // invalid
          continue;
        }
#endif
        if ( ( Idx[pix] >= 0.5 && Idx[pix]  <= N+0.5 && Idy[pix]  >= 0.5 && Idy[pix]  <= M+0.5 ))
        {
          if (occlusionsR[pix]) {fails[seg]++; continue;};
          freePixel[seg]++;
          //             if ( i < boxRadius || j < boxRadius || j >= M-boxRadius ) // j >= M-boxRadius || i >= N-boxRadius
          if ( i < boxRadius || j < boxRadius || j >= M-boxRadius || i >= N-boxRadius )
            censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, scores, hits, miss );
          else
#ifndef _nosse_
            censusBoxError4_short_SSE( segImg, img_s1, img_s2, pix, seg, intscores, hits, miss );
//             _mm_store_si128( &intscores[seg], _mm_add_epi16( _mm_load_si128( &intscores_pix[pix] ), _mm_load_si128( &intscores[seg] )) );
#else
            censusBoxError4_short( segImg, img_s1, img_s2, pix, seg, occluded, scores, hits, miss );
#endif
        }
        else
        {
          if (!occlusionsR[pix]) 
            fails[seg]++;
        }
      }
  }

#ifdef _use_census3_
  Scalar locThresh = thresh;
#else
#ifndef _use_census5_
  Scalar locThresh = thresh*4./24.;
#else
  Scalar locThresh = thresh*4./12.;
#endif
#endif
  //////////////////////////////////////////
  // adding up the stuff per segement involved:
  for (int i =0;i<nSegments;i++)
  {
    __m128i score = _mm_load_si128( &intscores[ i ] );
    // kill dumb highest 16 bits: never show up anyway so uncomment
    //      score = _mm_and_si128( score, _mm_set_epi16 (0,  0xffff,  0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff) );

#ifndef _sse_linux_
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
#ifdef OriginalVersion
    Scalar scale = 1.;
#else
    Scalar scale = 0.5;
#endif
#endif

    // idea restrict the maximal error per segment?: 1/32. ? sucks or same?
    //      locFreeScores[i] = min( hits[i] * thresh * 1./24., scores[i]*scale);
    locFreeScores[i] = scores[i]*scale;

#ifdef _data_truncation_
    locFreeScores[i] = min( locFreeScores[i], freePixel[i] * cutOffValue);
#endif
    //locScores[i]    += locFreeScores[i] + fails[i] * oobThresh;// fails are permenent as well as too large motion penalty
    Scalar failScore = (autoScore == NULL) ? oobThresh : (xtraPen + autoScore[i]);

    locScores[i]    += locFreeScores[i] + fails[i] * failScore;
  }
  //////////////////////////////////////////
};

template<class Scalar>
void genScore<Scalar>::
  censusBoxError4_short(const int* segImg, short *img1, short *img2, const int idsj, const int seg, 
  std::vector<Scalar>& scores, std::vector<int>& hits, std::vector<int>& miss )
{
  Scalar score (0);

  /// simplistic idea: just do something here
  short  i1 = img1[ idsj ];
  short  i2 = img2[ idsj ];
  short Repsi[3] = {163, 64, 0};

#ifdef _use_census3_
  const int steps = 4;
  const int dd[4] = { -1-M, -M, 1-M, -1};
  const int ww[4] = { 1, 1, 1, 1 };
  const int dy[4] = { -1,  0,  1,  -1};
  const int dx[4] = { -1, -1, -1,   0};
  Scalar scaleThresh = 1;
#else
#ifndef _use_census5_
  const int steps = 24;
  const int dd[24] = { -3-3*M, -2-3*M, -1-3*M, -3*M, 1-3*M, 2-3*M, 3-3*M, -3-2*M, -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, 3-2*M, -3-M, -2-M, -1-M, -M, 1-M, 2-M, 3-M, -3, -2, -1};
  const int ww[24] = {  3,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  2,  2,  3,  3,  2,  1,  1,  1,  2,  3,  3,  2,  1 };
  const int dy[24] = { -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1,  0,  1,  2,  3, -3, -2, -1 };
  const int dx[24] = { -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1,  -1, 0,  0,  0 };
  Scalar scaleThresh = 4./24.;
#else
  const int steps = 12;
  const int dd[24] = { -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, -2-M, -1-M, -M, 1-M, 2-M, -2, -1 };
  const int ww[24] = {  2,  2,  2,  2,  2,  2,  1,  1,  1,  2,  2,  1};
  const int dy[24] = { -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1};
  const int dx[24] = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0};
  Scalar scaleThresh = 4./12.;
#endif
#endif

  int lookId = ( refImg[idsj]*Scalar(255.0) );

#ifdef    OriginalVersion
  Scalar locThresh = thresh * scaleThresh;
#else
  Scalar locThresh = thresh * scaleThresh * 2.; // multiply all with 0.5 outside so here times 2
#endif

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


// could do the following: compute census score for top left quarter and store in intermediate buffer
// final score is sum of '4 quarters'
// would 1st with out checking run over all box pixel
// do top left census
// sum up with checks
// that would save half of census computations 
// 7x7 -> 4x4 window:
// x x x x
// x x x x
// x x x x
// x x x o
// so actualy census for a segment does do
// sum p\in pix sum q\in neigh(pix) |error(p, q, d(p,q))|
// due to symmetry in 2d this 4x4 block score above summarizes 
// all scores for o to x but also vice versi
// s.t each individual pixel to pixel score should be added to segment score of other pixel
// now exploiting this would be costly though (individual pixel checking sucks big time)

// actually a and c can be precomputed and stored per pixel -> 4 elements m128 per pixel 
// these can be loaded aligned as well 
// loop unroll ? does something at all ?
// exploit cache ? by loading new pixel into specific register - NO i load 8 at once anyway
/// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
template<class Scalar>
void genScore<Scalar>::
  censusBoxError4_short_SSE(const int* segImg, const short *img1, const short *img2, 
  const int idsj, const int seg, std::vector<__m128i>& intscores, 
  std::vector<int>& hits, std::vector<int>& miss )
#ifdef OriginalVersion
{
  /// simplistic idea: just do something here
  short  i1s = img1[ idsj ];
  short  i2s = img2[ idsj ];

  __m128i score = _mm_set1_epi16 (0);
  static const __m128i ones  = _mm_set1_epi16 (1);

  __m128i i1 = _mm_set1_epi16 (i1s);
  __m128i i2 = _mm_set1_epi16 (i2s);

  const int steps = 24;
  Scalar scaleThresh = 4./24.;
  const int dd[4] = {-3-3*M, -3-2*M, -3-M, -3};

  static const __m128i epsisf[4] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
    _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, 163, 64, 1) };

  static const __m128i iepsisf[4] = {_mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 ), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
    _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -32760, -32760, -32760, -32760, -163, -64, -1) };

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

#ifndef __Census__Compare__0__
    // equality case is neglected: should a<b, a>b not >= in the equation, so hurts with 0 comparison so wither set to 1 in epsis or not like this
    __m128i a = _mm_cmpgt_epi16(diff_i11, epsisf[cid]);
    __m128i b = _mm_cmplt_epi16(diff_i22, epsisf[cid]);
    __m128i c = _mm_cmplt_epi16(diff_i11, iepsisf[cid]);// not the opposite of above == case, why use diff and mdiff though
    __m128i d = _mm_cmpgt_epi16(diff_i22, iepsisf[cid]);

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
#else

  // here true full neighborrhood - just to be sure somehow
  /// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
  //  inline void censusBoxError4_short_SSE(const int* segImg, const short *img1, const short *img2, 
  //                              const int idsj, const int seg,// const std::vector<bool>& occluded, 
  //                              std::vector<__m128i>& intscores, std::vector<int>& hits, std::vector<int>& miss ) //std::vector<Scalar>& scores
{
  /// simplistic idea: just do something here
  short  i1s = img1[ idsj ];
  short  i2s = img2[ idsj ];

  __m128i score = _mm_set1_epi16 (0);//,0,0,0, 0,0,0,0);

//  __m128i ones  = _mm_set_epi16 (0, 1,1,1,1,1,1,1 );// 7x7 not 8x7 
  static const __m128i ones  = _mm_set1_epi16 (1);

  __m128i i1 = _mm_set1_epi16 (i1s);
  __m128i i2 = _mm_set1_epi16 (i2s);

#ifdef __Census__Compare__0__
  __m128i epsis[4] = {_mm_set_epi16 (32760, 0,0,0,0,0,0,0 ), _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0), 
    _mm_set_epi16 (32760, 0, 64, 0, 0, 0, 64, 0), _mm_set_epi16 (32760, 32760, 32760, 32760, 32760, 0, 64, 0)};
#else
  static const __m128i epsis[7] = {_mm_set_epi16 (32760, 1,1,1,1,1,1,1 ), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
    _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 1, 64, 163, 32760, 163, 64, 1), 
    _mm_set_epi16 (32760, 1, 64, 163, 163, 163, 64, 1), _mm_set_epi16 (32760, 1, 64, 64, 64, 64, 64, 1), 
    _mm_set_epi16 (32760, 1,1,1,1,1,1,1 )};

  static const __m128i iepsis[7] = {_mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 ), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
    _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -1, -64, -163, -32760, -163, -64, -1), 
    _mm_set_epi16 (-32760, -1, -64, -163, -163, -163, -64, -1), _mm_set_epi16 (-32760, -1, -64, -64, -64, -64, -64, -1), 
    _mm_set_epi16 (-32760, -1,-1,-1,-1,-1,-1,-1 )};
#endif

  const int runs = 7;// copy n times
  for (int cid=0; cid < runs; cid++)
  {
    int nid  = idsj + ddss[cid];

    // this is where the speedup is !
    __m128i i1_  = _mm_loadu_si128( (__m128i*) (&img1[nid]) );
    __m128i i2_  = _mm_loadu_si128( (__m128i*) (&img2[nid]) );

    __m128i  diff_i11 = _mm_sub_epi16(i1_, i1);
    __m128i  diff_i22 = _mm_sub_epi16(i2_, i2);

#ifndef __Census__Compare__0__
    // equality case is neglected: should a<b, a>b not >= in the equation, so hurts with 0 comparison so wither set to 1 in epsis or not like this
    __m128i a = _mm_cmpgt_epi16(diff_i11, epsis[cid]);
    __m128i b = _mm_cmplt_epi16(diff_i22, epsis[cid]);
    __m128i c = _mm_cmplt_epi16(diff_i11, iepsis[cid]);// not the opposite of above == case, why use diff and mdiff though
    __m128i d = _mm_cmpgt_epi16(diff_i22, iepsis[cid]);

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
  // all hits anyway
#ifndef _globalRescaling_
  __m128i tmp2 = _mm_load_si128( &intscores[seg] );
  _mm_store_si128( &intscores[seg], _mm_add_epi16( score, tmp2) );

  hits[seg]   += 48;
#else
  Scalar scaleThresh = 4./24.;
  Scalar locThresh = thresh * scaleThresh;
  score = _mm_hadd_epi16(score, score);
  score = _mm_hadd_epi16(score, score);
  score = _mm_hadd_epi16(score, score);
  scores[seg] += locThresh * score.m128i_i16[0];
  hits[seg]   += 48;
#endif
}
#endif


// per pixel function:
/// ids == global ids are used to look up the oob pixels PatchIdGenerator<Scalar>* genPatchIds // std::vector< int >* otherEdgeIds
/// lIds local ids, ids: global ids?, Mstep: height, Nstep width
template<class Scalar>
void genScore<Scalar>::
  computeFullScoresMappedCensusNew_short( short *img_s1 , short *img_s2, Scalar* Idx0, Scalar* Idy0, std::vector< std::pair<int,int> >& ids, std::vector< int >& lIds, int Mstep, int Nstep, const std::vector<bool>& oobPixR )
{
  const int csize = 3; // fixed census 7x7
  int nPixel = Nstep*Mstep; // size of the box

  std::vector<int> rowPtr(csize*2+1);
  for(int i=0;i<csize*2+1;i++)
    rowPtr[i] = (i-csize)*Mstep - csize;
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

  const int steps2 = 24;
  const int dd2[24] = { -2-2*M, -1-2*M, -2*M, 1-2*M, 2-2*M, -2-M, -1-M, -M, 1-M, 2-M, -2, -1, 1, 2, -2+M, -1+M, M, 1+M, 2+M, -2+2*M, -1+2*M, 2*M, 1+2*M, 2+2*M };
  const int ww2[24] = { 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  const int dy2[24] = { -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1, 1, 2, -2, -1,  0,  1,  2, -2, -1,  0,  1,  2 };
  const int dx2[24] = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2};

  short Repsi[3] = {163, 64, 0.};
  int cSize = ((2*csize+1)*(2*csize+1)-1)/2;// 24;//dxyC.size();
  Scalar dtaPen (thresh*2./Scalar(cSize));

  const int *dd, *ww, *dx, *dy;
  int steps=8;

  std::vector<__m128i> epsis(csize*2+1, _mm_set1_epi32(0));
  // epsilons:
  switch(csize)
  {
  case 1: 
    // setting 160 is about 1.25: 1.25*255*128 = 160
    epsis[0] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,163,163,163 );
    epsis[2] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,163,163,163 );
    epsis[1] = _mm_set_epi16 (32760, 32760, 32760,32760,32760,163,32760,163 );
    dd = dd1;
    ww = ww1;
    dx = dx1;
    dy = dy1;
    break;
  default:
  case 2:
    // idea load all 4, last per hand -> distribute, 3 runs only
    epsis[0] = _mm_set_epi16 (32760, 32760, 32760,  64,  64,  64,  64, 64 );
    epsis[1] = _mm_set_epi16 (32760, 32760, 32760,  64,   163,   163,   163, 64);
    epsis[2] = _mm_set_epi16 (32760, 32760, 32760,  64,   163, 32760, 163, 64);
    epsis[3] = _mm_set_epi16 (32760, 32760, 32760,  64,   163,   163,   163, 64);
    epsis[4] = _mm_set_epi16 (32760, 32760, 32760,  64,  64,  64,  64, 64 );
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
    epsis[2] = _mm_set_epi16 (32760, 0, 64, 163, 163, 163, 64, 0);
    epsis[3] = _mm_set_epi16 (32760, 0, 64, 163, 32760, 163, 64, 0);
    epsis[6] = _mm_set_epi16 (32760, 0,0,0,0,0,0,0 );
    epsis[5] = _mm_set_epi16 (32760, 0, 64, 64, 64, 64, 64, 0);
    epsis[4] = _mm_set_epi16 (32760, 0, 64, 163, 163, 163, 64, 0);
    dd = dd3;
    ww = ww3;
    dx = dx3;
    dy = dy3;
    steps= 48;
    break;
  }
  /////////////////////////// 

  // accumulated scores
  locScores.clear();
  locScores.resize(nPixel,0);

  // penalties
  locScores2.clear();
  locScores2.resize(nPixel,0);

  // correspondence inside image
  //    freePixel.clear();
  //    freePixel.resize(nPixel,0);

  std::vector<bool> occluded0;

  Scalar locScore0(0), locScore1(0);
  int fails0(0), fails1(0);

  for (int i = 0; i < nPixel ;i++)
  {
    if (lIds[i] < 0) continue; // local ids

    bool eval0 = true;
    if (maxDisp != N)// case of bounded motion checking
      locScores2[i]  = getPenaltyDisp( Idx0, eval0, occluded0, i, ids[i].first ) ;
    else
      locScores2[i]  = getPenaltyMot( Idx0, Idy0, eval0, occluded0, i, ids[i].first ) ;

    if ( !eval0 ) continue;

    // correspondence inside image boundaries
    // freePixel[i] = 1;

    if (i/Mstep < csize || i/Mstep >= Nstep-csize || i%Mstep < csize || i%Mstep >= Mstep-csize )
      censusBinaryNew_short( rowPtr, epsis, img_s1, img_s2, ids[i].first, locScores, i, dd, ww, dx, dy, steps, Repsi, Nstep, Mstep, cSize);
    else
      censusBinaryNew_shortSSE( rowPtr, epsis, img_s1, img_s2, ids[i].first, locScores, i, Repsi, dtaPen ) ;
  }
}


template<class Scalar>
Scalar genScore<Scalar>::
  getPenaltyDisp( Scalar* Idx, bool& eval, std::vector<bool>& occluded, int i, int gId )
{
  Scalar penalty(0);eval = true;

  int px = gId/M;
  int py = gId%M;
#ifdef __2dboundControl__
  // note that Idx should be 'close' to px at least if far away
  if ( ((maxDisp>=0) && (( Idx[i] > (px +1) ) || ( Idx[i] < (px +1 - maxDisp) ))) ||
    ((maxDisp< 0) && (( Idx[i] < (px +1) ) || ( Idx[i] > (px +1 - maxDisp) ))) ) // can not work with maxDisp < 0
  {
    penalty = Scalar(100.)*thresh; // invalid since in image
    eval = false;
  }
  else // valid disparity
#endif

    if (!(Idx[i] >= 0.5 && Idx[i] <= N+0.5)) // outside
    {
      eval = false; 
      penalty = oobThresh;
      penalty = (autoScore == NULL) ? oobThresh : (xtraPen+autoScore[i]);// non free variables are oob
    }

  return penalty;
}

template<class Scalar>
Scalar genScore<Scalar>::
  getPenaltyMot( Scalar* Idx, Scalar* Idy, bool& eval, std::vector<bool>& occluded, int i, int gId ) 
{
  Scalar penalty(0);eval = true;
#ifdef __2dboundControl__
  int px = gId/M;
  int py = gId%M;

  if ( (Idx[i]-(px +1)) * (Idx[i]-(px +1)) + (Idy[i]-(py+1)) * (Idy[i]-(py+1)) > maxMot*maxMot)
  {
    penalty = Scalar(100.)*thresh; // invalid - twice the time of per segment
    eval = false;
  }
  else
#endif
    if ( !( Idx[i]  >= 0.5 && Idx[i]  <= N+0.5 && Idy[i]  >= 0.5 && Idy[i]  <= M+0.5 ) )
    {
      penalty = oobThresh;
      penalty = (autoScore == NULL) ? oobThresh : (xtraPen+autoScore[i]);// non free variables are oob

      eval = false;
    }

  return penalty;
}

template<class Scalar>
char genScore<Scalar>::
  censusScore_short (short diff_a, short diff_b, short epsi )
{
  if (
    ( (diff_a > epsi) && (diff_b < epsi) ) ||
    ( (diff_a <-epsi) && (diff_b >-epsi) ) ||
    ( (diff_a < epsi) && (diff_b > epsi) ) ||
    ( (diff_a >-epsi) && (diff_b <-epsi) ) )
    return 1;
  return 0;
}


/// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
template<class Scalar>
void genScore<Scalar>::
  censusBinaryNew_short( const std::vector<int>& rowPtr, const std::vector<__m128i>& epsis, 
  const short *img0_a, short *img0_b, const int gId,
  std::vector<Scalar>& scores, const int i, const int *dd, const int *ww, 
  const int *dx, const int *dy, const int steps, const short Repsi[3], const int Mstep, const int Nstep, const int cSize)
{

  if (i<0 || i> Nstep*Mstep)
  {
    int yi = i%Nstep;
    int xi = i/Nstep;
    printf("Warning: i:%d - yi*Mstep + xi, Nstep %d, Mstep %d, Nstep*Mstep%d, gidX%d, gidY%d, xi%d,yi%d\n", i, yi*Mstep + xi, Nstep, Mstep, Nstep*Mstep, gId/N, gId%N, xi, yi);
    return;
  }

  // i : local position
  short  i0_a = img0_a[ i ];
  short  i0_b = img0_b[ i ];

  // dependent on gray value in center:
  //    int nGid;
  int score=0;
  int hits= 0;

  //    int cSize= 24;//dxyC.size();
  Scalar dtaPen (thresh*2./Scalar(cSize));

  int yi = i%Nstep;//after
  int xi = i/Nstep;

  for ( int j=0; j<steps; j++ )
  {
    // somehow that worked above, with this check
    if (yi+dy[j] >= 0 && yi+dy[j] < Nstep && xi+dx[j] >= 0 && xi+dx[j] < Mstep)//after && before
    {
      int nid = i+dx[j]*Nstep+dy[j];//afteer

      if (nid<0 || nid> Nstep*Mstep)
        printf("Warning: nid%d - (yi+dy[j])*Mstep + xi+dx[j]%d, Nstep %d, Mstep %d, Nstep*Mstep%d, gidX%d, gidY%d, xi%d,yi%d, dx%d,dy%d\n", nid, (yi+dy[j])*Mstep + xi+dx[j], Nstep, Mstep, Nstep*Mstep, gId/N, gId%N, xi, yi, dx[j], dy[j]);
      else
        if ( censusScore_short (img0_a[nid]-i0_a, img0_b[nid]-i0_b, Repsi[ww[j]-1]) )
          score++;
      hits++;
    }
  }

#ifdef _data_truncation_
  scores[i] = min( cutOffValue, score * steps * (dtaPen / std::max( Scalar(1.), Scalar(hits) )));
#else
  scores[i] = score * steps * (dtaPen / std::max( Scalar(1.), Scalar(hits) ));
#endif
  return;
}

/// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
template<class Scalar>
void genScore<Scalar>::
  censusBinaryNew_shortSSE( const std::vector<int>& rowPtr, const std::vector<__m128i>& epsis, 
  const short *img0_a, const short *img0_b, const int gId, std::vector<Scalar>& scores, 
  const int i, const short Repsi[3], const Scalar dtaPen )
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

#ifndef _sse_linux_
  scores[i] = dtaPen* (score.m128i_i16[0] + score.m128i_i16[1] + score.m128i_i16[2] + score.m128i_i16[3] + 
    score.m128i_i16[4] + score.m128i_i16[5] + score.m128i_i16[6] + score.m128i_i16[7]);
#else
  // should be done externally and once - as long as there are less than 32768 pixel in a segment
  __m128i tempA = _mm_add_epi16( _mm_shufflelo_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ), _mm_shufflehi_epi16 ( score, _MM_SHUFFLE(0,1,2,3) ) );
  score = _mm_add_epi16( tempA, _mm_unpackhi_epi16( tempA, tempA ) );
  tempA = _mm_add_epi16(  _mm_shufflelo_epi16( score, _MM_SHUFFLE(2,2,2,2) ), score );

#ifdef _data_truncation_
  scores[i] = min( cutOffValue, _mm_extract_epi16( tempA, 0 ) * dtaPen );
#else
  scores[i] = _mm_extract_epi16( tempA, 0 ) * dtaPen;
#endif

#endif
  return;
}

#endif