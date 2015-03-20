////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#ifndef __GEN_SCORE_HH
#define __GEN_SCORE_HH

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;

/*!
Accumulates the data scores. Currently used is the census transform
*/
template<typename Scalar> class genScore
{
public:

  typedef Math::Mat3x3T< Scalar>     M3;
  typedef Math::VectorT< Scalar, 3>  P3;
  typedef Math::VectorT<int, 4>      P4i;

  genScore( int N_, int M_, Scalar thresh_ = 0.4, Scalar oobThresh_ = Scalar(0.8) )
    : N(N_), M(M_), nSegments(0), nNormals(0), maxDisp(N_),  maxMot(M_*N_), thresh(thresh_), 
    oobThresh(oobThresh_), cutOffValue(0.8), autoScore(NULL), xtraPen(0.1)
  {
    ddss[0] = -3-3*M;ddss[1] = -3-2*M;ddss[2] = -3-M;ddss[3] = -3;ddss[4] = -3+M;ddss[5] = -3+2*M;ddss[6] = -3+3*M;
  };

  ~genScore(){};

  void setCutOff(Scalar _cutOffValue) {cutOffValue = _cutOffValue;};

  void setAutoScores( Scalar* _autoScore ){autoScore = _autoScore;};
  Scalar* getAutoScores( ){return autoScore;};

  void setRefImage(Scalar* refImg_) {refImg = refImg_;}; //buildEpsilonLookup();

  // should reflect the minimal depth possible - even better car: upper half should be 10m away or so
  void setMaxDisp(Scalar maxDisp_) {maxDisp = maxDisp_;};
  void setMaxMot (Scalar maxMot_)  {maxMot  = maxMot_;};

  Scalar getOobThresh() {return oobThresh;}; 
  Scalar getDataThresh() {return thresh;};
  Scalar*& getRefImgPtr() {return refImg;}; 
  int getMaxDisp() {return maxDisp;};
  int getMaxMot() {return maxMot;};

  /// that is used in Compressed as well !!!! the box is surrounding all segments in question
  void computeFullScoresCensus3OMP_Box( short *img_s1, short *img_s2, 
    const Scalar* const Idx , const Scalar* const Idy,
    int nSegments_, int* segImg, 
    const int bboxX1, const int bboxY1, const int bboxX2, const int bboxY2, 
    const std::vector<int>& vSegs, const std::vector<bool>& occlusionsR );

//  void buildEpsilonLookup();

  inline void censusBoxError4_short(const int* segImg, short *img1, short *img2, 
    const int idsj, const int seg, std::vector<Scalar>& scores, 
    std::vector<int>& hits, std::vector<int>& miss );

  /// int16 stuff is the smae - but use a per computed mask here, so that all 8 pixel per row are loaded at once, shuffle or blend with 0
  inline void censusBoxError4_short_SSE(const int* segImg, const short *img1, const short *img2, 
    const int idsj, const int seg, std::vector<__m128i>& intscores, 
    std::vector<int>& hits, std::vector<int>& miss );

  // per pixel function:
  /// ids == global ids are used to look up the oob pixels PatchIdGenerator<Scalar>* genPatchIds // std::vector< int >* otherEdgeIds
  /// lIds local ids, ids: global ids?, Mstep: height, Nstep width
  void computeFullScoresMappedCensusNew_short( short *img_s1 , short *img_s2, Scalar* Idx0, Scalar* Idy0, 
    std::vector< std::pair<int,int> >& ids, std::vector< int >& lIds, 
    int Mstep, int Nstep, const std::vector<bool>& oobPixR );

  std::vector<Scalar>& getScores()        {return locScores;};

  std::vector<Scalar>& getPenScores()     {return locScores2;};

  std::vector<Scalar>& getFreeScores()    {return locFreeScores;};

  std::vector<int>&    getFreeVariables() {return freePixel;};

  ////////////////////////////////////////////////////////////

  inline Scalar getPenaltyDisp( Scalar* Idx, bool& eval, std::vector<bool>& occluded, int i, int gId ); 

  inline Scalar getPenaltyMot( Scalar* Idx, Scalar* Idy, bool& eval, std::vector<bool>& occluded, int i, int gId );

  inline char censusScore_short (short diff_a, short diff_b, short epsi = 160);
  ////////////////

  void clearScores()
  {
    locScores.clear();
    locScores2.clear();
  }

  /// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
  inline void censusBinaryNew_short( const std::vector<int>& rowPtr, const std::vector<__m128i>& epsis, 
    const short *img0_a, short *img0_b, const int gId,
    std::vector<Scalar>& scores, const int i, const int *dd, const int *ww, 
    const int *dx, const int *dy, const int steps, const short Repsi[3], const int Mstep, const int Nstep, const int cSize);


  /// needs to store also the amount of hits later to be summed up i need to go in anyway - no eval, 1 vector
  inline void censusBinaryNew_shortSSE( const std::vector<int>& rowPtr, const std::vector<__m128i>& epsis, 
    const short *img0_a, const short *img0_b, const int gId, std::vector<Scalar>& scores, 
    const int i, const short Repsi[3], const Scalar dtaPen );

  ////////////////////////////////////////////////////////////
private:

  void censusBoxError4_short_SSE_inner( const int* segImg, 
    const short *img1, const short *img2, std::vector<__m128i>& intscores, 
    P4i innerbox, const std::vector<int>& vSegs, const Scalar* const Idx , const Scalar* const Idy);

  int N;
  int M;

  int nSegments;
  int nNormals;

  std::vector<Scalar> locScores;
  std::vector<Scalar> locScores2;
  std::vector<Scalar> locFreeScores; 
  std::vector<int>    freePixel;

  char* oMap;

  Scalar maxDisp;
  Scalar maxMot;

  Scalar thresh;
  Scalar oobThresh;

  std::vector<Scalar> lookUpThresh;

  /// to get the brightness at the original pixel in the image
  Scalar* refImg;

  Scalar cutOffValue;

  Scalar* autoScore;
  Scalar xtraPen;

  // init by first init and kept fix so know M
  int ddss[7];
};
///////////////////////////////////

#if defined(INCLUDE_TEMPLATES) && !defined(__genScore__cpp)
#include "genScore.cpp"
#endif

#endif // __GEN_SCORE_HH
