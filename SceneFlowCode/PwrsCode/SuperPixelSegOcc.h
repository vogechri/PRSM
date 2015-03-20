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

#ifndef __SUPERPIXESEG2____
#define __SUPERPIXESEG2____

////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#define _occHandlingOn_

#include "DataDefinitions.h"

#include "PatchSegmentHandler.h"
#include "AccumData.h"

#include "EvalEnergyStore.h"
#include "EvalEnergyFull2FramePixel.h"

#include "OcclusionMappingBufferSegments.h"
#include "BinaryConversion.h"

#include <algorithm>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T
//#include "QPBO.h" 
#include "../QPBO-v1.3.src/QPBO.h"

#include <ctime>

using namespace std;
using namespace Math;

// at this rate the improvements are too small too see an effect and it takes ages to finish
#define __Small_Improvement__ 0.501

/////////////////////////////////////////////////////////////////////////////

// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void SuperPix2Frames( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

  bool doOcclusionStuff = true;

  size_t N,M;

  int patchSize = 30;
  int censusSize = 3;
  Scalar thresh  = 0.1;
  Scalar lambda  = 0.02;
  Scalar theta   = 100.0;
  Scalar dispMax = 125;
  Scalar motMax  = 220; // the maximal motion in the scene
  Scalar outOfBoundsThresh = 0.14;

  Scalar rotJump(15.);//   setRotJump (rotJump)
  Scalar depthJump(sqrt(1300.));//   setRotJump (rotJump)
  Scalar occThresh(0.05);

  char* occMapR (NULL);
  char* occMapLT(NULL);
  char* occMapRT(NULL);

  int* segImg(NULL), *origSegImg(NULL);
  int* occSegImg(NULL);
  //  int* oobSegs = NULL;
  int* oobSegsR  = NULL;
  int* oobSegsLT = NULL;
  int* oobSegsRT = NULL;
  int fullIterations = 3;

  const mxArray* occList  (NULL);
  const mxArray* occListL (NULL);
  const mxArray* occListR (NULL);
  Scalar* Kr_(NULL);

  Scalar* img1    = (Scalar*)  mxGetPr(prhs[0]);
  Scalar* img2    = (Scalar*)  mxGetPr(prhs[1]);
  Scalar* img3    = (Scalar*)  mxGetPr(prhs[2]);
  Scalar* img4    = (Scalar*)  mxGetPr(prhs[3]);

  Scalar* p2d_    = (Scalar*)  mxGetPr(prhs[4]); // K_l^-1
  Scalar* Kl_     = (Scalar*)  mxGetPr(prhs[5]); // K_l
  const mxArray* Segments  = (prhs[6]);

  Scalar* normals = (Scalar*) mxGetPr(prhs[ 7]);
  thresh          = (Scalar) *mxGetPr(prhs[ 8]);
  lambda          = (Scalar) *mxGetPr(prhs[ 9]);
  theta           = (Scalar) *mxGetPr(prhs[10]);

  //  const mxArray* edges    = (prhs[11]);
  //  const mxArray* weights  = (prhs[12]);

  Scalar* halfEdgeX(NULL);
  Scalar* halfEdgeY(NULL);
  Scalar* halfEdgeXY(NULL);
  Scalar* halfEdgeiXY(NULL);

  // optional: always 0,0,0 if not needed
  Scalar* mc_       =  (Scalar*) mxGetPr(prhs[13]);

  // optional:
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[14]);
  Scalar* centers   =  (Scalar*) mxGetPr (prhs[15]);
  Scalar* MC_       =  (Scalar*) mxGetPr (prhs[16]);

  Scalar* homs_     =  (Scalar*) mxGetPr (prhs[17]);
  Scalar* vNoms_    =  (Scalar*) mxGetPr (prhs[18]);

  if (nrhs > 19)
    rotJump         =  (Scalar) *mxGetPr(prhs[19]);

  if (nrhs > 20)
    depthJump       =  (Scalar) *mxGetPr(prhs[20]);

  if (nrhs > 21)
    dispMax         =  (Scalar) *mxGetPr(prhs[21]);

  //////////////////
  if (nrhs > 22)
    outOfBoundsThresh  =  (Scalar) *mxGetPr(prhs[22]);

  if (nrhs > 23)
    origSegImg      =  (int*) mxGetPr(prhs[23]);

  if (nrhs > 24)
    halfEdgeX       =  (Scalar*) mxGetPr(prhs[24]);

  if (nrhs > 25)
    halfEdgeY       =  (Scalar*) mxGetPr(prhs[25]);

  if (nrhs > 26)
    halfEdgeXY      =  (Scalar*) mxGetPr(prhs[26]);

  if (nrhs > 27)
    halfEdgeiXY     =  (Scalar*) mxGetPr(prhs[27]);

  if (nrhs > 28)
    Kr_             = (Scalar*)  mxGetPr(prhs[28]); // K_r
  else
    Kr_ = Kl_;

  if (nrhs > 29)
    patchSize       = (Scalar)  *mxGetPr(prhs[29]); // ps

  if (nrhs > 30)
    motMax     =  (Scalar) *mxGetPr(prhs[30]);

  if (nrhs > 31)
    occThresh   =  (Scalar) *mxGetPr(prhs[31]);

  if (nrhs > 32)
    fullIterations =  (int) *mxGetPr(prhs[32]);

  if (nrhs > 33)
    oobSegsR   =  (int*) mxGetPr(prhs[33]);

  if (nrhs > 34)
    oobSegsLT   =  (int*) mxGetPr(prhs[34]);

  if (nrhs > 35)
    oobSegsRT  =  (int*) mxGetPr(prhs[35]);

  //////////////////////////
  printf("fullIterations %d \n", fullIterations);
  printf("patchSize %d \n", patchSize);
  printf("censusSize %d \n", censusSize);
  printf("rotJump %f \n", rotJump);
  printf("depthJump %f \n", depthJump);
  printf("dispMax %f \n", dispMax);
  printf("motMax %f \n", motMax);
  printf("occThresh %f \n", occThresh);
  printf("outOfBoundsThresh %f \n", outOfBoundsThresh);
  printf("data thresh %f \n", thresh);
  printf("motion smoothness theta %f \n", theta);
  printf("smooth lambda %f \n", lambda);

  Scalar* output, *output2 = NULL;
  Scalar* output3, *output4 = NULL;

  M           = (size_t) mxGetM(prhs[0]);
  N           = (size_t) mxGetN(prhs[0]);

  if (origSegImg != NULL)
  {
    segImg = (int*) malloc(N*M*sizeof(int));
    memcpy(segImg, origSegImg, N*M*sizeof(int));
  }

  int nNormals = max( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );
  int nStep    = min( mxGetM(prhs[ 7]), mxGetN(prhs[ 7]) );

  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  mwSize dims[2];dims[0] = N;dims[1] = M;
  // Create the output array
  {
    // segment map later
    plhs[0]      =  mxCreateDoubleMatrix(N, M, mxREAL);
    output  = (Scalar*) mxGetPr(plhs[0]);

    memset(output, 0, sizeof(Scalar) * N*M);
  }

  printf ("nSegments %d, nNormals %d \n", nSegments, nNormals);
  //  mexEvalString("drawnow");
  ///////////////////////////////////////////////////////

  ///////// generate oob maps per view:
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
    ////////////////////////
    ///////////////////////////////////////////////////////

    std::srand ( unsigned ( std::time(0) ) );

    PatchIdGenerator<Scalar> genPatchIds(N,M, patchSize, centers, origSegImg, Kl_, nSegments, segImg, img1, true, censusSize);

    int curId, pStartX, pStartY, pEndX, pEndY;

    EvalEnergyStore<Scalar> enStore;
    enStore.prepareNormals( normals, nNormals, nStep );
    //  enStore.prepareRotTra( rotations, nNormals );
    enStore.setRotTra( rotations, nNormals );

    genHomography<Scalar> gHom1 (Kl_, MC_, mc_, Kr_);//t,t, left right
    genHomography<Scalar> gHom2 (Kl_); // t,t+1 left
    genHomography<Scalar> gHom3 (Kl_, MC_, mc_, Kr_);// t,t+1 left right


    gHom1.setNormals( enStore.getNormals() );

    gHom2.setNormals( enStore.getNormals() );
    gHom2.setRotTra(  enStore.getTra(), enStore.getRot() );

    gHom3.setNormals( enStore.getNormals() );
    gHom3.setRotTra(  enStore.getTra(), enStore.getRot() );

    gHom2.setCenters( centers );
    gHom3.setCenters( centers );

    genWarp<Scalar> gWarp0(N,M, img1);// no warps -> orig image
    gWarp0.setDummy();
    genWarp<Scalar> gWarp1(N,M, img2, p2d_);// do warps, I2:r,t
    gWarp1.setHom( &gHom1 );//rt

#ifndef __do_depth_only__
    genWarp<Scalar> gWarp2(N,M, img3, p2d_);// do warps  I3:l,t+1
    genWarp<Scalar> gWarp3(N,M, img4, p2d_);// do warps  I4:r,t+1
    gWarp2.setHom( &gHom2 );//lt1
    gWarp3.setHom( &gHom3 );//rt1
#endif

    accumulateWarp<Scalar> aW01(N,M,thresh,outOfBoundsThresh, censusSize);
    aW01.setMaxDisp(dispMax); // if set a test is done before adding up values
    aW01.setWarp1(&gWarp1);  aW01.setWarp2(&gWarp0);// to enable testing of maxDisp
    aW01.setRefImage(img1);

#ifndef __do_depth_only__
    accumulateWarp<Scalar> aW02(N,M,thresh,outOfBoundsThresh, censusSize);
    accumulateWarp<Scalar> aW13(N,M,thresh,outOfBoundsThresh, censusSize);
    accumulateWarp<Scalar> aW23(N,M,thresh,outOfBoundsThresh, censusSize);
    aW02.setWarp1(&gWarp2);  aW02.setWarp2(&gWarp0);// takes indices from the first warp
    aW02.setMaxMot(motMax);   // if set the maximal motion is restricted
    aW13.setWarp1(&gWarp1);  aW13.setWarp2(&gWarp3);
    aW23.setWarp1(&gWarp2);  aW23.setWarp2(&gWarp3);
    aW02.setRefImage(img1);
    aW13.setRefImage(img1);
    aW23.setRefImage(img1);
#endif

    accumulate2FrameData<Scalar> aFD(N,M);
    // 4(3) images to be warped
    aFD.addWarp(&gWarp1);
    // 4 equations from images:
    aFD.addAccumulateWarp(&aW01);  

#ifndef __do_depth_only__
    aFD.addWarp(&gWarp2);  aFD.addWarp(&gWarp3);  
    aFD.addAccumulateWarp(&aW02);  aFD.addAccumulateWarp(&aW13);  aFD.addAccumulateWarp(&aW23);
#endif
    aFD.addWarp(&gWarp0);//must be last

    typedef genWarp<Scalar>::P3 P3;

    ////////////////////////////////////

    EvalEnergyFullFramePixel<Scalar> evalSmooth(N,M);
    evalSmooth.setRotJump  (rotJump);
    evalSmooth.setDepthJump(depthJump);

    enStore.flipNormals();
    evalSmooth.setNormals( enStore.getNormals() );
    evalSmooth.setRotTra( enStore.getRot(), enStore.getTra() );
    evalSmooth.setGamma(1.0);

    evalSmooth.setRotWeight( theta );
    evalSmooth.set2dMotionMatrix ( Kl_, Kr_, MC_, mc_, 1 ); // jump at 0.5 pixel -> needs to set jump smaller

    evalSmooth.setK ( Kl_ );  
    evalSmooth.setCenters( centers );
    evalSmooth.prepare( N*M, nSegments  );
    evalSmooth.setSegImg( N, M, segImg, nSegments );
    evalSmooth.setHalfEdges(halfEdgeX, halfEdgeY);

    evalSmooth.setCrossHalfEdges(halfEdgeXY, halfEdgeiXY);
    evalSmooth.setHom( homs_, nSegments );

    OcclusionMapBuffer<Scalar> testOPix((2*patchSize+1)*(2*patchSize+1),  doOcclusionStuff);
    OcclusionMapBuffer<Scalar> testOPix2((2*patchSize+1)*(2*patchSize+1), doOcclusionStuff);
    OcclusionMapBuffer<Scalar> testOPix3((2*patchSize+1)*(2*patchSize+1), doOcclusionStuff);

    // that the buffer for the reference image - to generate the idx, idy arrays
    OcclusionMapBuffer<Scalar> dummyBuffer((2*patchSize+1)*(2*patchSize+1), N,M);
    dummyBuffer.init();

    testOPix.SetOccPenalty(occThresh);
    testOPix.setSegImg( N,M, segImg, (int) (mxGetM(Segments) * mxGetN(Segments)) );

    testOPix.setP2d( p2d_ );
    testOPix.setHomvNom( homs_, vNoms_ );
    testOPix.init();

    testOPix2.setP2d( p2d_ );
    testOPix2.setSegImg( N,M, segImg, (int) (mxGetM(Segments) * mxGetN(Segments)) );
    testOPix2.setHomvNom( &homs_[9*nSegments], &vNoms_[3*nSegments] );
    testOPix2.init();

    testOPix3.setP2d( p2d_ );
    testOPix3.setSegImg( N,M, segImg, (int) (mxGetM(Segments) * mxGetN(Segments)) );
    testOPix3.setHomvNom( &homs_[2*9*nSegments], &vNoms_[2*3*nSegments] );
    testOPix3.init();

    /// reihenfolge oben: 1,2,3,0
    aFD.addOcclusionMappingBuffer( &testOPix  ); // has to be added to represent the image, so gWarpi for all i, the first one therefore has to have identity idx, idy
#ifndef __do_depth_only__
    aFD.addOcclusionMappingBuffer( &testOPix2 ); // has to be added to represent the image, so gWarpi for all i, the first one therefore has to have identity idx, idy
    aFD.addOcclusionMappingBuffer( &testOPix3 ); // has to be added to represent the image, so gWarpi for all i, the first one therefore has to have identity idx, idy
#endif
    aFD.addOcclusionMappingBuffer( &dummyBuffer ); // has to be added to represent the image, so gWarpi for all i, the first one therefore has to have identity idx, idy

    gWarp0.request2ndWarp();
    gWarp1.request2ndWarp();
#ifndef __do_depth_only__
    gWarp2.request2ndWarp();
    gWarp3.request2ndWarp();
#endif

    genPatchIds.initMappingsNew( -1, false );
    std::vector<int> allGlobIdsTrial_; // ids of all the pixel in the patch

    /// trial means here the ids of the free pixel
    allGlobIdsTrial_ = genPatchIds.getAllGlobIds();
    std::vector<int> allLocIdsTrial_; // ids for the edges
    allLocIdsTrial_    = genPatchIds.getLocIds();
    aFD.computeScorePairWindowCensusNew_short( &allGlobIdsTrial_, &allLocIdsTrial_, M, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT, true );// for each pairing a different score

    unaryDataScoresBuffer<Scalar> dataBuffer0( N,M, aFD.getPartialScores( 0 ), censusSize );
#ifndef __do_depth_only__
    unaryDataScoresBuffer<Scalar> dataBuffer1( N,M, aFD.getPartialScores( 1 ), censusSize );
    unaryDataScoresBuffer<Scalar> dataBuffer2( N,M, aFD.getPartialScores( 2 ), censusSize );
    unaryDataScoresBuffer<Scalar> dataBuffer3( N,M, aFD.getPartialScores( 3 ), censusSize );
#endif

    std::vector<Scalar> scores0;
    std::vector<Scalar> scores1;
    std::vector<Scalar> scores2;
    std::vector<Scalar> scores3;
    std::vector<char> edgeDataScores0;
    std::vector<char> edgeDataScores1;
    std::vector<char> edgeDataScores2;
    std::vector<char> edgeDataScores3;
    std::vector<bool> fixedIds;
    //--
    HigherOrderConverter<Scalar> hoc(nSegments);
    //--
    int allNodes (0), allEdges(0), allUnlabeled(0), allRuns(0);
    int numEdges = 4*8*(2*patchSize+1)*(2*patchSize+1);
    int numNodes = 4*(2*patchSize+1)*(2*patchSize+1);
    QPBO<int> q(numNodes, numEdges ); // max number of nodes & edges

    Scalar scale = 1;
    int maxValue = 1 << 27;

    Scalar ENERGY = 1000000;
    Scalar OLD_ENERGY = ENERGY;
    srand ( time(NULL) );
    clock_t start = std::clock();

    int all_non_sol = 0;
    int all_impIts  = 0;
    int non_sub = 0;

    for (int nIterations = 0; nIterations < fullIterations; nIterations++)
    {
      // could possibly change :)
      int nProps = genPatchIds.getNProposals();
      std::vector<int> order(nProps);
      for ( int i = 0; i < nProps; i++ )
        order[i]=i;

      // full sweep:
      random_shuffle( order.begin(), order.end() );

      std::vector<Scalar> locScores0;
      std::vector<Scalar> locScores1;
      std::vector<Scalar> locScores2;
      std::vector<Scalar> locScores3;


#ifndef _DEBUG
      //    for ( int i = 0; i < 10; i++ )
      for ( int i = 0; i < order.size(); i++ )
#else
      for ( int i = 0; i < order.size(); i++ ) //min(int(order.size()), 10); i++ )
#endif
      {
        int trial = order[i];
        int newTrial = trial;
#ifdef __QUIET__
        if ( i%50 == 49 ) printf("it: %d \n", i);
#endif
        if (i< order.size()-1)
          newTrial = order[i+1];

        genPatchIds.initMappingsNew( trial, false );
        int nOrigVars = genPatchIds.getNVars();
        if (nOrigVars < 1) continue;

        std::vector<int> allGlobIdsTrial; // ids of all the pixel in the patch
        std::vector<int> globIdsTrial; // ids of the free pixel in the patch
        std::vector<int> locIdsTrial;  // local ids of the free pixel in the patch
        std::vector<int> patchEdgeIds; // ids for the edges 
        std::vector<int> otherPatchEdgeIds; // ids for the edges

        /// trial means here the ids of the free pixel
        allGlobIdsTrial = genPatchIds.getAllGlobIds();
        globIdsTrial    = genPatchIds.getGlobIds();
        locIdsTrial     = genPatchIds.getLocIds();//.getLocIds();
        patchEdgeIds    = genPatchIds.getEdgeIds();
        otherPatchEdgeIds = genPatchIds.getOtherEdgeIds();

        fixedIds     = genPatchIds.getFixedIds();

        curId = genPatchIds.getCurrentPatch( pStartX, pEndX, pStartY, pEndY );

        testOPix.initNewSeg   (curId, pStartX, pEndX, pStartY, pEndY);
#ifndef __do_depth_only__
        testOPix2.initNewSeg  (curId, pStartX, pEndX, pStartY, pEndY);
        testOPix3.initNewSeg  (curId, pStartX, pEndX, pStartY, pEndY);
#endif
        dummyBuffer.initNewSeg(curId, pStartX, pEndX, pStartY, pEndY);
        // now we also have the coordinates of the warps to compute the data term

        testOPix.setLocalIds( locIdsTrial );
#ifndef __do_depth_only__
        testOPix2.setLocalIds( locIdsTrial );
        testOPix3.setLocalIds( locIdsTrial );
#endif

        testOPix.generateOcclusionLists( curId, pStartX, pEndX, pStartY, pEndY );
        testOPix2.generateOcclusionLists( curId, pStartX, pEndX, pStartY, pEndY );
        testOPix3.generateOcclusionLists( curId, pStartX, pEndX, pStartY, pEndY );

        aFD.computeScorePairWindowCensusNew_short( &allGlobIdsTrial, &locIdsTrial, pEndY-pStartY, oobOutsideLR, oobOutsideLLT, oobOutsideRRT, oobOutsideLTRT );// for each pairing a different score

        // scores sre allPartialscores which is the penalty scores!
        locScores0.resize(allGlobIdsTrial.size()*2); std::copy( aFD.getPartialScores( 0 ).begin(), aFD.getPartialScores( 0 ).end(), locScores0.begin());

#ifndef __do_depth_only__
        locScores1.resize(allGlobIdsTrial.size()*2); std::copy( aFD.getPartialScores( 1 ).begin(), aFD.getPartialScores( 1 ).end(), locScores1.begin());
        locScores2.resize(allGlobIdsTrial.size()*2); std::copy( aFD.getPartialScores( 2 ).begin(), aFD.getPartialScores( 2 ).end(), locScores2.begin());
        locScores3.resize(allGlobIdsTrial.size()*2); std::copy( aFD.getPartialScores( 3 ).begin(), aFD.getPartialScores( 3 ).end(), locScores3.begin());
#endif

        dataBuffer0.constructDataScore( scores0, locScores0, globIdsTrial, genPatchIds);
        dataBuffer1.constructDataScore( scores1, locScores1, globIdsTrial, genPatchIds);
        dataBuffer2.constructDataScore( scores2, locScores2, globIdsTrial, genPatchIds);
        dataBuffer3.constructDataScore( scores3, locScores3, globIdsTrial, genPatchIds);

        // scores0 ,scores1, etc. now contain the scores for the current, current penalties, trial, trial penalties in one vector, in consecutive order for each view

        //      Scalar worstSmoothScore = 10000;
        // takes half the time !!!
        //    Scalar worstSmoothScore = evalSmooth.compute_score_combiDepth( curId, theta, pStartX, pEndX, pStartY, pEndY, locIdsTrial, globIdsTrial.size());
        Scalar worstSmoothScore = evalSmooth.compute_score_combiDepth( curId, pStartX, pEndX, pStartY, pEndY, locIdsTrial, globIdsTrial.size());

        const std::vector<int>& idk=evalSmooth.getIdk();
        const std::vector<int>& idl=evalSmooth.getIdl();

        const std::vector<Scalar>& f11=evalSmooth.getF11();
        const std::vector<Scalar>& f00=evalSmooth.getF00();
        const std::vector<Scalar>& f10=evalSmooth.getF10();
        const std::vector<Scalar>& f01=evalSmooth.getF01();

        // to be added to the data term
        const std::vector<Scalar>& f0=evalSmooth.getF0();
        const std::vector<Scalar>& f1=evalSmooth.getF1();
        {
          Scalar occPenalty = occThresh; // best thing is to do nothing !
          /////////////////////////////////////

          // global dataa : dataBuffer0.dataScore
          //        printf(" convert occ ");mexEvalString("drawnow");
          // free .. is always 1 - or 0 if penalty is not 0! convert this internally
          BinaryConverter<Scalar> bc_0( occPenalty, testOPix.getOccList0(), testOPix.getOccList1(), scores0 );//list,list, local penalties, local penalties
          BinaryConverter<Scalar> bc_1( occPenalty, testOPix2.getOccList0(), testOPix2.getOccList1(), scores1 );
          BinaryConverter<Scalar> bc_2( occPenalty, testOPix.getOccList0(), testOPix.getOccList1(), testOPix3.getOccList0(), testOPix3.getOccList1(), scores2 );
          BinaryConverter<Scalar> bc_3( occPenalty, testOPix2.getOccList0(), testOPix2.getOccList1(), testOPix3.getOccList0(), testOPix3.getOccList1(), scores3 );
          //        printf (" bc init ");mexEvalString("drawnow"); 
          //scores0 and locScores0 first true data scores, than penalties, 
          bc_0.setPenaltiesPerSegment( scores0 );
          bc_1.setPenaltiesPerSegment( scores1 );
          bc_2.setPenaltiesPerSegment( scores2 );
          bc_3.setPenaltiesPerSegment( scores3 );

          bc_0.addOcclusionListOuterPixels( testOPix.getOuterList() );
          bc_1.addOcclusionListOuterPixels( testOPix2.getOuterList() );
          bc_2.addOcclusionListOuterPixels( testOPix.getOuterList(), testOPix3.getOuterList() );
          bc_3.addOcclusionListOuterPixels( testOPix2.getOuterList(), testOPix3.getOuterList() );
          //        printf (" bc addocclist ");mexEvalString("drawnow"); 
          ////////////////////////
          bc_0.generateHigherOrderTerms();
          bc_1.generateHigherOrderTerms();
          bc_2.generateHigherOrderTerms();
          bc_3.generateHigherOrderTerms();
          //        printf (" bc hot ");mexEvalString("drawnow"); 
          bc_0.generateHigherOrderTerms_outer( dataBuffer0.dataScore );
          bc_1.generateHigherOrderTerms_outer( dataBuffer1.dataScore );
          bc_2.generateHigherOrderTerms_outer( dataBuffer2.dataScore );
          bc_3.generateHigherOrderTerms_outer( dataBuffer3.dataScore );
          //        printf (" bc hot_out ");mexEvalString("drawnow"); 
          // it is often a good idea to merge submods and nonsubmods in advance before building terms from this
          // adding score together. e.g. bc_o and bc_2, since here is an overlap !
          hoc.setNSegments( nOrigVars );

          //        printf(" mergeSubs ");mexEvalString("drawnow");
          //        printf (" #bc0.nsm %d, #bc2.nsm %d, #bc0.sm %d, #bc2.sm %d ", bc_0.get_nSubmods().size(), bc_2.get_nSubmods().size(), bc_0.get_Submods().size(), bc_2.get_Submods().size());mexEvalString("drawnow");
          hoc.addMergeNonSubs_Subs( bc_0.get_nSubmods(), bc_2.get_nSubmods(), bc_0.get_Submods(), bc_2.get_Submods() );
          //        printf (" #bc1.nsm %d, #bc3.nsm %d, #bc1.sm %d, #bc3.sm %d ", bc_1.get_nSubmods().size(), bc_3.get_nSubmods().size(), bc_1.get_Submods().size(), bc_3.get_Submods().size());mexEvalString("drawnow");
          hoc.addMergeNonSubs_Subs( bc_1.get_nSubmods(), bc_3.get_nSubmods(), bc_1.get_Submods(), bc_3.get_Submods() );
          //        printf (" hoc addM ");mexEvalString("drawnow"); 
          hoc.addUnaries (bc_0.getUnaries());
          hoc.addUnaries (bc_1.getUnaries());
          hoc.addUnaries (bc_2.getUnaries());
          hoc.addUnaries (bc_3.getUnaries());
          //        printf (" hoc conv ");mexEvalString("drawnow"); 
          hoc.convert();
          //        printf (" hoc conv done ");mexEvalString("drawnow"); 

          // takes ages in practice '0 &&' should be deleted 
          if ( 0 && (hoc.get_nVars() > numNodes || hoc.get_nEdges() > numEdges) )
          {
#ifdef __QUIET__
            printf(" Skipping :%d vars, %d edges\n", hoc.get_nVars(), hoc.get_nEdges() );mexEvalString("drawnow");
#endif
            i++;
            q.Reset();
            hoc.clear();
            continue;
          }
          else
          {
#ifdef __QUIET__
            printf(" node:%d, edge:%d", hoc.get_nVars(), hoc.get_nEdges());mexEvalString("drawnow");
#endif
          }
          q.AddNode( hoc.get_nVars() );
        }
        ////////////

        Scalar worstDataScore = 0;
        for ( int j = 0; j< nOrigVars; j++ )
        {                
          Scalar dt0 = scores0[j]+scores1[j]+scores2[j]+scores3[j] + scores0[j+nOrigVars]+scores1[j+nOrigVars]+scores2[j+nOrigVars]+scores3[j+nOrigVars];
          Scalar dt1 = scores0[j+nOrigVars*2]+scores1[j+nOrigVars*2]+scores2[j+nOrigVars*2]+scores3[j+nOrigVars*2] + scores0[j+nOrigVars*3]+scores1[j+nOrigVars*3]+scores2[j+nOrigVars*3]+scores3[j+nOrigVars*3];

          worstDataScore += max(dt0,dt1);
        }


        Scalar scale = (Scalar)(maxValue) / (lambda * worstSmoothScore + worstDataScore);
        Scalar scaleSmo = scale * lambda;

        for ( int j = 0; j< nOrigVars; j++ )
        {
          Scalar dt0 = scores0[j]+scores1[j]+scores2[j]+scores3[j] + scores0[j+nOrigVars]+scores1[j+nOrigVars]+scores2[j+nOrigVars]+scores3[j+nOrigVars];
          Scalar dt1 = scores0[j+nOrigVars*2]+scores1[j+nOrigVars*2]+scores2[j+nOrigVars*2]+scores3[j+nOrigVars*2] + scores0[j+nOrigVars*3]+scores1[j+nOrigVars*3]+scores2[j+nOrigVars*3]+scores3[j+nOrigVars*3];

          q.AddUnaryTerm(j, scale*(dt0)+scaleSmo*f0[j], scale*(dt1)+scaleSmo*f1[j]);
        }

        for( int j=0; j < idk.size(); j++ )
        {
          if (f00[j]+f11[j] > f10[j]+f01[j])
          {
            //          printf("NON submodularity detected (diff:%f) (00=%f,11=%f, 10=%f,01=%f), (curr:(%d,%d), trial:%d)\n", f00[j]+f11[j] - f10[j]-f01[j], f00[j],f11[j], f10[j],f01[j] ,currentSolution[idk[j]], currentSolution[idl[j]], trial);
            non_sub++;
          }
          q.AddPairwiseTerm(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]); // add term (x+1)*(y+2)
        }

        // change QPBO for occlusion handling
        int nVars = hoc.get_nVars();//occHandler.getNVars();
		allNodes += nVars;
	    allRuns++;
        if (nVars > nOrigVars || hoc.get_nEdges() > 0 || hoc.get_b0().size() + hoc.get_b1().size() > 0)
        {
          //        printf(" occ done"); mexEvalString("drawnow");
          std::vector< std::pair<int, Scalar> >& occ_b0 = hoc.get_b0();
          std::vector< std::pair<int, Scalar> >& occ_b1 = hoc.get_b1();

          std::vector< edgePenalty<Scalar> >& occ_b00 = hoc.get_b00();
          std::vector< edgePenalty<Scalar> >& occ_b01 = hoc.get_b01();
          std::vector< edgePenalty<Scalar> >& occ_b10 = hoc.get_b10();
          std::vector< edgePenalty<Scalar> >& occ_b11 = hoc.get_b11();

		  allEdges += occ_b00.size() + occ_b01.size() + occ_b10.size() + occ_b11.size() + f01.size();
		  
          for (int j=0; j< occ_b00.size();j++)
		  {
            if (occ_b00[j].f11 > 0)
              non_sub++;
            q.AddPairwiseTerm(occ_b00[j].seg_i, occ_b00[j].seg_j, scale*occ_b00[j].f11,  0,0,0);
		  }
          for (int j=0; j< occ_b01.size();j++)
		  { 
		    if (occ_b01[j].f11 < 0)
              non_sub++;
			q.AddPairwiseTerm(occ_b01[j].seg_i, occ_b01[j].seg_j, 0, scale*occ_b01[j].f11, 0,0);
		  }
          for (int j=0; j< occ_b10.size();j++)
		  {
		    if (occ_b10[j].f11 < 0)
              non_sub++;
		    q.AddPairwiseTerm(occ_b10[j].seg_i, occ_b10[j].seg_j, 0,0, scale*occ_b10[j].f11, 0);
		  }
          for (int j=0; j< occ_b11.size();j++)
		  {
            if (occ_b11[j].f11 > 0)
              non_sub++;		  
            q.AddPairwiseTerm(occ_b11[j].seg_i, occ_b11[j].seg_j, 0,0,0, scale*occ_b11[j].f11 );
		  }
		  
          for (int j=0; j< occ_b0.size();j++)
            q.AddUnaryTerm(occ_b0[j].first, scale*(occ_b0[j].second), 0);
          for (int j=0; j< occ_b1.size();j++)
            q.AddUnaryTerm(occ_b1[j].first, 0, scale*(occ_b1[j].second));

          // otherwise all edges should be unique/called only once!
          //        printf(" premerge edges"); mexEvalString("drawnow");
          q.MergeParallelEdges();
          //        printf(" edges merge "); mexEvalString("drawnow");
        }
        hoc.clear();
        /////////////////////////////

        // this is the situation before expansion
        std::vector<int> currentSolution(nVars,0);
        q.Solve();
        //      printf(" solved "); mexEvalString("drawnow");
        q.ComputeWeakPersistencies();
        //      printf(" pers "); mexEvalString("drawnow");

        int non_sol=0;int doImprove = 0;
        std::vector<int> noSolutions;noSolutions.reserve(nVars);

        // find noSol only in pixel variables - others not important
        for (int j = 0; j < nVars;j++)
        {
          int x = q.GetLabel(j);
          if (x<0)
          {
            doImprove = 1;
            non_sol++;
            noSolutions.push_back(j);
          }
          if (x==1) {currentSolution[j] = 1;};
        }
        bool successPI = false; // eval energy after -> how much is the gain ??

		allUnlabeled += non_sol;
		
        int impIts = 0;
        if (noSolutions.size() < 300)
        {
		// new for testing
#define useProbe
          //#undef useProbe
#ifdef useProbe
          int doProbe = 1;int allImpIts = 0;
          const int maxProbe =3;//3
          if (doProbe && doImprove)
          {
            for (int j = 0; j < nVars;j++)
            {
              int x = q.GetLabel(j);
              if (x>=0)
                q.SetLabel(j,x);
            }

            QPBO<int>::ProbeOptions options;
            int* mapping = (int*) malloc( sizeof(int) * nVars );
            q.Probe( mapping, options );
            q.ComputeWeakPersistencies();

            int* mapping2 = (int*) malloc( sizeof(int) * nVars );
            for (int ii=1;ii < maxProbe; ii++)
            {
#ifdef __QUIET__
              printf(" %d ", ii );      mexEvalString("drawnow");
#endif
              q.Probe( mapping2, options );
              q.ComputeWeakPersistencies();
              q.MergeMappings( nVars, mapping, mapping2 );
            }
            free(mapping2);
            const int maxInner = 10;//10
            const int maxReps = 20; int imp = 0;
            while (imp < maxInner && impIts < maxReps)
            {
              impIts++;
#ifdef __QUIET__
              for(imp=0; imp<maxInner; imp++ )
                if ( q.Improve() )  {printf("Improved after probing ");break;};
#endif
              allImpIts += imp;
              if (imp == maxInner) break;// out of while
            }

            for (int j=0;j < noSolutions.size(); j++)
            {
              int pos = noSolutions[j];
              int x = q.GetLabel( mapping[pos]/2 );
              if( (x>=0) && ((x+ mapping[pos]) % 2) )
              {
                successPI = true;
                currentSolution[pos] = 1;
                allImpIts++;
                impIts++;
              }
            }
            free( mapping );
          }
#endif
        }
        else
        {
		// new for testing
#define useImprove
          //#undef useImprove
#ifdef useImprove
          if (doImprove && noSolutions.size() < 600)
          {
#ifdef __QUIET__
            printf(" Improve : %d undef ", noSolutions.size() );mexEvalString("drawnow");
#endif
            for (int j = 0; j < nVars;j++)
            {
              int x = q.GetLabel(j);
              if (x>=0)
                q.SetLabel(j,x);
            }

            const int maxReps = 20;const int maxInner = 10;int allImpIts = 0;
            int imp=0;//int impIts = 0;
            if (doImprove)
              while (imp < maxInner && impIts < maxReps)
              {
                impIts++;
                for(imp=0; imp<maxInner; imp++ )
                  if ( q.Improve() )  {break;printf("Improved after Improving\n");};

                allImpIts += imp;
                if (imp == maxInner) break;// out of while
              }

              for (int j = 0; j < nOrigVars; j++)
              {
                int x = q.GetLabel(j);
                if (x==1) {currentSolution[j] = 1;successPI=true;};
              }
          }
#endif
        }
        q.Reset();

        currentSolution.resize( nOrigVars );

        testOPix.finishNewSeg( curId, globIdsTrial, currentSolution);
        testOPix2.finishNewSeg( curId, globIdsTrial, currentSolution);
        testOPix3.finishNewSeg( curId, globIdsTrial, currentSolution);

        dataBuffer0.updateDataScore( locScores0, globIdsTrial, genPatchIds, currentSolution );
#ifndef __do_depth_only__
        dataBuffer1.updateDataScore( locScores1, globIdsTrial, genPatchIds, currentSolution );
        dataBuffer2.updateDataScore( locScores2, globIdsTrial, genPatchIds, currentSolution );
        dataBuffer3.updateDataScore( locScores3, globIdsTrial, genPatchIds, currentSolution );
#endif

        for (int j=0;j<nOrigVars;j++)
          if (currentSolution[j]) 
            segImg[globIdsTrial[j]] = curId;

        all_non_sol += non_sol;
        all_impIts  += impIts;
      }
#ifdef __QUIET__
      printf("Round:%d non_sol: %d, improved by probe %d\n", nIterations, all_non_sol, all_impIts);
      //    mexEvalString("drawnow");
#endif
      // find patches which do not exist any more, find those which 'hit' their borders and place new seed there not a new center, a seed
      if (nIterations+1 < fullIterations)
        genPatchIds.refreshMappings();

    }
#ifdef __QUIET__
    printf("QPBO took %f\n", double(std::clock()-start)/CLOCKS_PER_SEC);
#endif
    ////
    if (segImg != NULL)
    {
      for ( int i = 0; i < N*M; i++ )
        output[ i ] = (Scalar) (segImg[i]);
      free( segImg );
    }

    if (nlhs > 1)
    {
      plhs[1] =  mxCreateDoubleMatrix(N, M, mxREAL);
      output2  = (Scalar*) mxGetPr(plhs[1]);
      testOPix.getCurrentOcclusions( output2 );
    }
    if (nlhs > 2)
    {
      plhs[2] =  mxCreateDoubleMatrix(N, M, mxREAL);
      output3  = (Scalar*) mxGetPr(plhs[2]);
      testOPix2.getCurrentOcclusions( output3 );
    }
    if (nlhs > 3)
    {
      plhs[3] =  mxCreateDoubleMatrix(N, M, mxREAL);
      output4  = (Scalar*) mxGetPr(plhs[3]);
      testOPix3.getCurrentOcclusions( output4 );
    }

	if (nlhs > 4)
	{
		mwSize dims[2];
		dims[0] = 1;
		dims[1] = 8;
		plhs[4] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		Scalar* outputI = (Scalar*) mxGetPr(plhs[4]);
		outputI[0] = all_impIts; outputI[1] = non_sub; outputI[2] = allUnlabeled; 
		outputI[3] = allNodes; outputI[4] = allEdges; outputI[5] = allRuns;
		outputI[6] = allUnlabeled;
		outputI[7] = all_non_sol;
  }
	
	
    mexPrintf("non_sol: %d, improved by probe %d - regular non_sub:%d\n", all_non_sol, all_impIts, non_sub);
}
#undef __do_depth_only__
#undef __Small_Improvement__
#undef useImprove
#undef useProbe
#endif// __SUPERPIXESEG____
