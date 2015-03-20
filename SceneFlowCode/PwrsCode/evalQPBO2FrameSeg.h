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

/// here we allow proposals to be moved to theri neighbours
/// the selection is taken randomly, we pick a segment, distribute the motion
/// to the neeighbours and fix these
/// finally we stop if no segment can be pixked any more
/// we do this N times
/// also possible Growing seeds like in patchMatch
#ifndef __QPBO_SEG_2FRAME__
#define __QPBO_SEG_2FRAME__
////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include "DataDefinitions.h"
#include "AccumData.h"

#ifdef _approximateSmoothing_
#include "EvalEnergyFull2Frame.h"
#else
#include "EvalEnergyReal2Frame.h"
#endif

#ifdef _motionChecking_
#include "SensibleMotionCheck.h"
#endif

#include "EvalEnergyStore.h"
#include "LocalProposalHandler.h"

#include <algorithm>
#include <list>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T
#include "../QPBO-v1.3.src/QPBO.h"

#include <ctime>

using namespace std;
using namespace Math;

#ifdef NDEBUG
#define _GRID_SIZE_ 150
#else
#define _GRID_SIZE_ 60
#endif

#define __Small_Improvement__ 0.501

// also use cross view
//#define _use_5th_view

//////////////////////////////////////////////////////////////////////////////////////

template<typename Scalar>
Scalar getEnergy ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda)//, std::vector<bool> oobS,  std::vector<bool> oobT)
{
  Scalar score(0.);
  for ( int j = 0; j< nSegments; j++ )
    score += dataScores[j];

  evalSmooth.compute_score_Fuse( current, current );
//  evalSmooth.compute_score_combiDepth( 0, current);//, oobT, oobS );

  const std::vector<Scalar>& f00=evalSmooth.getF00();
  for ( int i = 0; i < f00.size(); i++)
    score+= f00[i] * lambda;

  return score;
}

template<typename Scalar>
Scalar getDataEnergy ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores )
{
  Scalar score(0.);
  for ( int j = 0; j< nSegments; j++ )
    score += dataScores[j];

  return score;
}

template<typename Scalar>
Scalar getSmoothEnergy ( std::vector<int>& current, int nSegments, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda)//, std::vector<bool> oobS, std::vector<bool> oobT)
{
  Scalar score(0.);

  evalSmooth.compute_score_Fuse( current, current );
//  evalSmooth.compute_score_combiDepth( 0, current);//, oobT, oobS );
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  for ( int i = 0; i < f00.size(); i++)
    score+= f00[i] * lambda;

  return score;
}


// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void evalQPBO2FramesNew( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;
  typedef genWarp<Scalar>::P3 P3;

  size_t N,M;
  size_t elements;
  char buffer[256];

  Scalar thresh = 0.4;
  Scalar lambda = 0.02;
  Scalar theta  = 100.0;
  Scalar dispMax = 256.;
  Scalar pixJump = 1.0;
  Scalar rotJump(20.);//   setRotJump (rotJump)
  Scalar depthJump(20.);//   setRotJump (rotJump)
  Scalar occThresh(0.05/0.75);
  Scalar motMax = 350; // the maximal motion in the scene

  Scalar* halfEdgeX(NULL);
  Scalar* halfEdgeY(NULL);
  Scalar* halfEdgeXY(NULL);
  Scalar* halfEdgeiXY(NULL);
  int* origSegImg(NULL);

  Scalar* Kr_(NULL);
  int* segImg(NULL);

  Scalar* img1    = (Scalar*)  mxGetPr(prhs[0]);
  Scalar* img2    = (Scalar*)  mxGetPr(prhs[1]);
  Scalar* img3    = (Scalar*)  mxGetPr(prhs[2]);
  Scalar* img4    = (Scalar*)  mxGetPr(prhs[3]);

  Scalar* p2d_    = (Scalar*)  mxGetPr(prhs[4]); // K_l^-1
  Scalar* K_      = (Scalar*)  mxGetPr(prhs[5]); // K_r
  const mxArray* Segments  = (prhs[6]);

  Scalar* normals = (Scalar*) mxGetPr(prhs[ 7]);
  thresh          = (Scalar) *mxGetPr(prhs[ 8]);
  lambda          = (Scalar) *mxGetPr(prhs[ 9]);
  theta           = (Scalar) *mxGetPr(prhs[10]);

  const mxArray* edges    = (prhs[11]);
  const mxArray* weights  = (prhs[12]);

  // optional: always 0,0,0 if not needed
  Scalar* mc_       =  (Scalar*) mxGetPr(prhs[13]);// right cam

  int* oobSegsR  = NULL;
  int* oobSegsLT = NULL;
  int* oobSegsRT = NULL;

  // optional:
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[14]);
  Scalar* centers   =  (Scalar*) mxGetPr (prhs[15]);
  Scalar* MC_       =  (Scalar*) mxGetPr (prhs[16]);//right cam

  if (nrhs > 17)
    rotJump       =  (Scalar) *mxGetPr(prhs[17]);

  if (nrhs > 18)
    depthJump     =  (Scalar) *mxGetPr(prhs[18]);

  if (nrhs > 19)
    dispMax     =  (Scalar) *mxGetPr(prhs[19]);

  //////////////////
  if (nrhs > 20)
    occThresh   =  (Scalar) *mxGetPr(prhs[20]);

  if (nrhs > 21)
    segImg      =  (int*) mxGetPr(prhs[21]); // for census - not flipped

  if (nrhs > 22)
    oobSegsR   =  (int*) mxGetPr(prhs[22]);

  if (nrhs > 23)
    oobSegsLT   =  (int*) mxGetPr(prhs[23]);

  if (nrhs > 24)
    oobSegsRT  =  (int*) mxGetPr(prhs[24]);

  if (nrhs > 25)
    motMax     =  (Scalar) *mxGetPr(prhs[25]);

  if (nrhs > 26)
    halfEdgeX       =  (Scalar*) mxGetPr(prhs[26]);

  if (nrhs > 27)
    halfEdgeY       =  (Scalar*) mxGetPr(prhs[27]);

  if (nrhs > 28)
    halfEdgeXY      =  (Scalar*) mxGetPr(prhs[28]);

  if (nrhs > 29)
    halfEdgeiXY     =  (Scalar*) mxGetPr(prhs[29]);

  if (nrhs > 30)
    origSegImg      =  (int*) mxGetPr(prhs[30]);

  if (nrhs > 31)
    Kr_             = (Scalar*)  mxGetPr(prhs[31]); // K_r
  /// if Kr == Kl
  if (Kr_==NULL) Kr_ = K_;

  //////////////////////////

  Scalar outOfBoundsThresh = occThresh;
  std::srand ( unsigned ( std::time(0) ) );

  printf("rotJump %f \n", rotJump);
  printf("depthJump %f \n", depthJump);
  printf("dispMax %f \n", dispMax);
  printf("motMax %f \n", motMax);
  printf("occThresh %f \n", occThresh);
  printf("outOfBoundsThresh %f \n", outOfBoundsThresh);
  printf("data thresh %f \n", thresh);
  printf("motion smooth theta %f \n", theta);
  printf("general smoothness lambda %f \n", lambda);

  Scalar* output, *output2 = NULL;
  Scalar* output3, *output4 = NULL;

  M           = (size_t) mxGetM(prhs[0]);
  N           = (size_t) mxGetN(prhs[0]);

  int nNormals = max( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );
  int nStep    = min( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );

  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  printf ("nSegments %d, nNormals %d \n", nSegments, nNormals);
  //  mexEvalString("drawnow");

  std::vector<bool> oobOutsideLR(M*N, false);
  if (oobSegsR != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsR[i] ) oobOutsideLR[i]=true;//disp

  std::vector<bool> oobOutsideLLT(M*N, false);
  if (oobSegsLT != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsLT[i] ) oobOutsideLLT[i]=true;//flow

  // maybe use later
  std::vector<bool> oobOutsideLRT(M*N, false);
  if (oobSegsRT != NULL)
    for (int i= 0; i< M*N; i++)
      if ( oobSegsRT[i] ) oobOutsideLRT[i]=true;//left right t+1

  /// these could be merged with left, right t+1 
  std::vector<bool> oobOutsideRRT(M*N, false); //right right t+1
  if (oobSegsRT !=NULL)
  {
    for (int i= 0; i< M*N; i++)
      if ( oobSegsR[i] || oobSegsRT[i] ) oobOutsideRRT[i]=true;
  }
  else 
    if (oobSegsR != NULL)
    {
      for (int i= 0; i< M*N; i++)
        if ( oobSegsR[i] ) oobOutsideRRT[i]=true;
    }

  /// rrt
  std::vector<bool> oobOutsideLTRT(M*N, false);  //left t+1 - right t+1
  if (oobSegsRT != NULL  )
  {
    for (int i= 0; i< M*N; i++)
      if ( oobSegsLT[i] || oobSegsRT[i] ) oobOutsideLTRT[i]=true;
  }
  else // no third view
  {
    if (oobSegsRT != NULL)
      for (int i= 0; i< M*N; i++)
        if ( oobSegsLT[i] ) oobOutsideLTRT[i]=true;
  }

  mwSize dims[2];dims[0] = N;dims[1] = M;
  // Create the output array
  {
    // scores
    plhs[0]      =  mxCreateDoubleMatrix(nSegments, 1, mxREAL); //mxCreateDoubleMatrix(N, M, mxREAL); // not double  but ...
    output  = (Scalar*) mxGetPr(plhs[0]);
  }

  ///////////////////////////////////////////////////////

  EvalEnergyStore<Scalar> enStore;
  enStore.prepareNormals( normals, nNormals, nStep );
  enStore.setRotTra( rotations, nNormals );

  genHomography<Scalar> gHom1 (K_, MC_, mc_, Kr_);//t,t, left right
  genHomography<Scalar> gHom2 (K_); // t,t+1 left
  genHomography<Scalar> gHom3 (K_, MC_, mc_, Kr_);// t,t+1 left right

  gHom1.setNormals( enStore.getNormals() );
  gHom2.setNormals( enStore.getNormals() );
  gHom2.setRotTra(  enStore.getTra(), enStore.getRot() );  
  gHom3.setNormals( enStore.getNormals() );
  gHom3.setRotTra(  enStore.getTra(), enStore.getRot() );

  gHom2.setCenters( centers );
  gHom3.setCenters( centers );

  genWarp<double> gWarp1(N,M, img2, p2d_);// do warps
  genWarp<double> gWarp2(N,M, img3, p2d_);// do warps
  genWarp<double> gWarp3(N,M, img4, p2d_);// do warps
  genWarp<double> gWarp0(N,M, img1);// no warps -> orig image

  gWarp1.setHom( &gHom1 );//rt
  gWarp2.setHom( &gHom2 );//lt1
  gWarp3.setHom( &gHom3 );//rt1

  accumulateWarp<Scalar> aW01(N,M,thresh,outOfBoundsThresh); // disparity
  accumulateWarp<Scalar> aW02(N,M,thresh,outOfBoundsThresh); // flow
  accumulateWarp<Scalar> aW13(N,M,thresh,outOfBoundsThresh); // flow right
  accumulateWarp<Scalar> aW23(N,M,thresh,outOfBoundsThresh); // disparity t+1

#ifdef _use_5th_view
  accumulateWarp<Scalar> aW03(N,M,thresh,outOfBoundsThresh); // cross view lt->rt1
  aW03.setWarp1(&gWarp3);  aW03.setWarp2(&gWarp0);// to enable testing of maxDisp
  aW03.setRefImage(img1);
#endif


  aW01.setMaxDisp(dispMax); // if set a test is done before adding up values
  aW01.setWarp1(&gWarp1);  aW01.setWarp2(&gWarp0);// to enable testing of maxDisp
  aW02.setWarp1(&gWarp2);  aW02.setWarp2(&gWarp0);// takes indices from the first warp
  aW02.setMaxMot(motMax);
  aW13.setWarp1(&gWarp1);  aW13.setWarp2(&gWarp3);
  aW23.setWarp1(&gWarp2);  aW23.setWarp2(&gWarp3);

  accumulate2FrameData<double> aFD(N,M);
  // 4(3) images to be warped
  aFD.addWarp(&gWarp1); 
  // 4 equations from images:
  aFD.addAccumulateWarp(&aW01);  

  aW01.setRefImage(img1);
  aW02.setRefImage(img1);
  aW13.setRefImage(img1);
  aW23.setRefImage(img1);

  aFD.addWarp(&gWarp2);  aFD.addWarp(&gWarp3);  aFD.addWarp(&gWarp0);
  aFD.addAccumulateWarp(&aW02);  aFD.addAccumulateWarp(&aW13);  aFD.addAccumulateWarp(&aW23);

#ifdef _use_5th_view
  aFD.addAccumulateWarp(&aW03);
#endif

#ifdef __QUIET__
  std::clock_t startP(std::clock());
#endif

  ///////////////////////
  /// find a cover for the whole reference view s.t. the expansion area for the proposals is at lest _GRID_SIZE_^2 (150x150) pixel
  int runs = nNormals/nSegments;

  if( runs < 1)
    printf("Too few normals/rotations, less than segments\n");

  LocalProposalHandler<Scalar> LPH( N, M, _GRID_SIZE_, centers, K_, nSegments, segImg );
  LPH.generateProposals(edges);
  std::vector<int> proMap = LPH.getProposalMap();
  int nDatas = runs*LPH.getNGrids();

  printf("GZ: %d, getNGrids: %d  ratio:%4.2f eff:%2.3f ", _GRID_SIZE_, LPH.getNGrids(), Scalar(nSegments)/Scalar(LPH.getNGrids()), LPH.efficiency());

  std::vector< LocalProposalHandler<Scalar>::P4i >& lbb = LPH.getLargeBBox();
  std::vector< std::list< std::pair<int,int> > >&   psp = LPH.getPageSegPair();


  if (runs < 20)
  {
    int sz = proMap.size();
    proMap.resize( runs * proMap.size());
    for (int i=1;i<runs;i++)
      for (int k=0;k < sz;k++)
      proMap[i*sz + k] = proMap[(i-1)*sz + k]+nSegments;
  }
  else
  {
    printf("Alarm: too many proposals: %d normals vs %d segments\n", nNormals, nSegments);
  }
  ///////////////////////

  std::vector<int> currentSolution( nSegments, 0 );
  std::vector<int> oldSolution( nSegments, 0 );

  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      currentSolution[i]=i;

#ifdef __QUIET__
  std::clock_t stopP(std::clock());
  printf("Preprocess took %f\n", double(stopP-startP)/CLOCKS_PER_SEC );
  std::clock_t start(std::clock());
#endif

  aFD.computeAllVecScores_boxed (Segments, nNormals, segImg, psp, lbb, runs*nSegments, proMap, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, oobOutsideLRT);

  std::vector<Scalar>& allData = aFD.getAllScores();
  Scalar bestData = aFD.getBestConfiguration( nDatas, currentSolution );

#ifdef __QUIET__
  std::clock_t stop(std::clock());
  printf("Data took %f\n", double(stop-start)/CLOCKS_PER_SEC );
#endif

  std::vector<Scalar> currentData        (nSegments, 0);
  std::vector<int>    currentPropSolution(nSegments, 0);

  // fix the assignments:
  for ( int i = 0;i < nSegments;i++ )
  {
    currentPropSolution[i] = nSegments*currentSolution[i] + i;
    currentSolution[i]     = proMap[currentPropSolution[i]];
    currentData[i] = allData[ currentPropSolution[i] ];
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  // now smoothness cost setup:

  EvalEnergyFullFrame<Scalar> evalSmooth;
  evalSmooth.setRotJump  (rotJump);
  evalSmooth.setDepthJump(depthJump);

  evalSmooth.setGamma ( 1.0 );
  evalSmooth.setRotWeight( theta );
//  evalSmooth.set2dMotionMatrix ( K_, Kr_, MC_, mc_, pixJump ); // jump at 0.5 pixel -> needs to set jump smaller
#ifdef _approximateSmoothing_
  evalSmooth.set2dMotionMatrix ( K_, Kr_, MC_, mc_, pixJump ); // jump at 0.5 pixel -> needs to set jump smaller
#else
  evalSmooth.set2dMotionMatrix ( K_, MC_, mc_, pixJump ); // jump at 0.5 pixel -> needs to set jump smaller
#endif
  
  enStore.flipNormals();
  evalSmooth.setNormals( enStore.getNormals() );
  evalSmooth.setRotTra( enStore.getRot(), enStore.getTra() );

  assert ( segImg != NULL );

  evalSmooth.setHalfEdges(halfEdgeX, halfEdgeY);
  if ( halfEdgeiXY != NULL)
    evalSmooth.setCrossHalfEdges(halfEdgeXY, halfEdgeiXY);
//  evalSmooth.prepareOwnWeights( nSegments, edges, centers, origSegImg, N, M );
#ifdef _approximateSmoothing_
    evalSmooth.prepareOwnWeights( nSegments, edges, centers, origSegImg, N, M );
#else
    evalSmooth.initMapping( nSegments, edges, centers, origSegImg, N, M );
//  evalSmooth.buildWeights( origSegImg, N, M, nSegments); // not really
#endif

  int numEdges = evalSmooth.edges_num( );

#ifdef _motionChecking_
  MotionCheck<Scalar> mChecker( p2d_, nSegments, nNormals );
  mChecker.setPenalty(8. * thresh); // in dependence on data
  mChecker.setNormals( enStore.getNormals() );
  mChecker.setRotTra( enStore.getRot(), enStore.getTra() );
  mChecker.setWidthHeight(N,M);
  mChecker.setK( K_ );  
#endif
  
  Scalar OLD_ENERGY = std::numeric_limits<Scalar>::max();
  Scalar ENERGY(0);

  /// fixed anyway
  const std::vector<int>& idk=evalSmooth.getIdk();
  const std::vector<int>& idl=evalSmooth.getIdl();
  const std::vector<Scalar>& f11=evalSmooth.getF11();
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  const std::vector<Scalar>& f10=evalSmooth.getF10();
  const std::vector<Scalar>& f01=evalSmooth.getF01();     
  std::vector<int> noSolutions;noSolutions.reserve(nSegments);
  
  int maxValue = 1 << 29;
  int non_sub = 0;
  int non_sol = 0;

  QPBO<int> q(nSegments, numEdges); // max number of nodes & edges

  ENERGY = getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda);
#ifdef __QUIET__
  printf("------------ Start Energy is %f --------------\n", ENERGY);
  printf("Start Energies (full=data+smooth) %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, currentData ), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda) );
#endif
  std::clock_t fullStart(std::clock());
  // fuse all proposals in a row

  std::vector<int> order(nDatas, 0);
  for ( int i = 0;i < nDatas;i++ )
    order[i]=i;

  while( ENERGY+__Small_Improvement__ < OLD_ENERGY )
  {
    // full sweep over all proposal:
    random_shuffle( order.begin(), order.end() );
    for ( int ord = 0; ord < order.size(); ord++ )
    {
      int propId = order[ord];
      std::vector<int>    trialVec    (proMap.begin()+nSegments*propId,  proMap.begin()+nSegments*(propId+1) );
      std::vector<Scalar> dataTrialVec(allData.begin()+nSegments*propId, allData.begin()+nSegments*(propId+1) );

      Scalar worstDataScore(0);
      Scalar bestDataScore(0);

#ifdef _motionChecking_	  
	  mChecker.computeGeometryConsistency( Segments, centers, currentSolution );
      std::vector<Scalar> motionCheckScores1 = mChecker.getScores();
      mChecker.computeGeometryConsistency( Segments, centers, trialVec );
      std::vector<Scalar> motionCheckScores2 = mChecker.getScores();
#endif
	  
      for (int k=0; k< nSegments;k++)
      {
        worstDataScore += max( currentData[k], dataTrialVec[k] );
        bestDataScore  += min( currentData[k], dataTrialVec[k] );
      }

#ifdef __DEBUG_ENERGIES__
      for ( int k = 0;k < nSegments;k++ )
        oldSolution[k]=currentSolution[k];
#endif

      q.AddNode(nSegments);
      Scalar worstSmoothScore = evalSmooth.compute_score_Fuse( currentSolution, trialVec );

      Scalar scale = (Scalar)(maxValue) / (lambda * worstSmoothScore + worstDataScore);

      Scalar scaleSmo = scale * lambda;
      // add unary data terms:
      for ( int j = 0; j< nSegments; j++ )
      {
        q.AddUnaryTerm(j, scale*currentData[j] , scale*dataTrialVec[j] );
      }

#ifdef _motionChecking_
      for ( int j = 0; j< nSegments; j++ )
      {
        q.AddUnaryTerm(j, scale*motionCheckScores1[j], scale*motionCheckScores2[j] );
      }
#endif
	  
      for( int j=0; j < idk.size(); j++ )
      {
        if (f00[j]+f11[j] > f10[j]+f01[j])
          non_sub++;

        if (f11[j] != f00[j] || f00[j] != f01[j] || f00[j] != f10[j] )
          q.AddPairwiseTerm(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);
      }

      q.Solve();
      q.ComputeWeakPersistencies();

#define useProbe
#ifndef useProbe
      for (int j = 0; j < nSegments;j++)
      {
        int x = q.GetLabel(j);
        if (x<0)
          non_sol++;

        if (x==1) {
          currentPropSolution[j] = nSegments*propId + j;
          currentData[j]         = allData [currentPropSolution[j] ];
          currentSolution[j]     = proMap  [currentPropSolution[j] ];
        };
      }

#else
      // try to improve if something is not solved completely:
      int doImprove = 0;
#ifdef _DEBUG
      std::vector<int> noSolution ( nSegments, 0 );
#endif
      noSolutions.clear();
#ifdef _DEBUG
      std::vector<int> new1old0(nSegments, 0);
      std::list<int>   newListed;
#endif
      int allNew = 0;
      for (int j = 0; j < nSegments;j++)
      {
        oldSolution[j] = currentSolution[j];
        int x = q.GetLabel(j);
        if (x<0)
        {
          doImprove = 1;
          non_sol++;
#ifdef _DEBUG
          noSolution[j]++; 
#endif
          noSolutions.push_back(j);
        }
        if (x==1) {
          currentPropSolution[j] = nSegments*propId + j;
          currentData[j]         = allData [currentPropSolution[j] ];
          currentSolution[j]     = proMap  [currentPropSolution[j] ];

          allNew++;
#ifdef _DEBUG
          new1old0[j]=1;
          newListed.push_back(j);
#endif
        };
      }

      int doProbe = 1;int impIts =0;int allImpIts = 0;
      const int maxProbe =3;//3
      if (doProbe && doImprove)
      {
        int nVars = nSegments;
        for (int j = 0; j < nVars;j++)
        {
          int x = q.GetLabel(j);
          if (x>=0)
            q.SetLabel(j,x);
        }

        QPBO<int>::ProbeOptions options;
        int* mapping = (int*) malloc( sizeof(int) * nVars );
        q.Probe( mapping, options );

        int* mapping2 = (int*) malloc( sizeof(int) * nVars );
        for (int i=1;i < maxProbe; i++)
        {
          q.Probe( mapping2, options );
          q.MergeMappings( nVars, mapping, mapping2 );
        }
        free(mapping2);
        const int maxInner = 10;//10;
        const int maxReps = 20; int imp = 0;
        while (imp < maxInner && impIts < maxReps)
        {
          impIts++;
          for(imp=0; imp<maxInner; imp++ )
            if ( q.Improve() )  {break;printf("Improved after probing\n");};

          allImpIts += imp;
          if (imp == maxInner) break;// out of while
        }

        for (int j=0;j < noSolutions.size(); j++)
        {
          int pos = noSolutions[j];
          int x = q.GetLabel( mapping[pos]/2 );
          if( (x>=0) && ((x+ mapping[pos]) % 2) )
          {
            currentPropSolution[pos] = nSegments*propId + pos;
            currentData[pos]         = allData [currentPropSolution[pos] ];
            currentSolution[pos]     = proMap  [currentPropSolution[pos] ];

            allImpIts++;
            impIts++;
          }
        }
        free( mapping );
      }
#endif

      if (ord == order.size()-1)
      {
        OLD_ENERGY = ENERGY;
        ENERGY = getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda);
#ifdef __QUIET__
        printf("Energies (full=data+smooth) comp %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, currentData), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda) );
#endif
		}

      q.Reset();
    }
  }

#ifdef __QUIET__
  std::vector<int> realSolution( nSegments, 0 );
  if ( nNormals >= nSegments )
    for ( int i = 0;i < nSegments;i++ )
      realSolution[i]=i;
  // not working since flipped
  printf("Energies (full=data+smooth) comp %f = (%f,%f) 0:N %f = (%f,%f)\n", getEnergy ( currentSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( currentSolution, nSegments, currentData ), getSmoothEnergy ( currentSolution, nSegments, evalSmooth, theta, lambda), getEnergy ( realSolution, nSegments, currentData, evalSmooth, theta, lambda), getDataEnergy ( realSolution, nSegments, currentData ), getSmoothEnergy ( realSolution, nSegments, evalSmooth, theta, lambda) );

  printf("QPBO took %f\n", double(std::clock()-fullStart)/CLOCKS_PER_SEC);
  printf("ENERGY: %f  NON submodularity detected %d, no solutions: %d\n", ENERGY, non_sub, non_sol);
#endif
  for (int i=0;i < currentSolution.size(); i++)
    output[i] =currentSolution[i];
}
#undef useProbe
#undef __Small_Improvement__
#endif 
