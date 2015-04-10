/// here we allow proposals to be moved to theri neighbours
/// the selection is taken randomly, we pick a segment, distribute the motion
/// to the neeighbours and fix these
/// finally we stop if no segment can be pixked any more
/// we do this N times
/// also possible Growing seeds like in patchMatch
#ifndef __VCSF_SEG_GC__
#define __VCSF_SEG_GC__
///////////////////////////////////////////////////////////////
//// Compute the mapping P from segments to moving planes  ////
///////////////////////////////////////////////////////////////

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "DataDefinitionsVC.h"

#ifdef _qpboVersion_
#include "QPBOCutting.h" // use QPBO to solve
#else
#include "GraphCutting.h"  // use LSA Aux to solve
#endif

#include <numeric>
#include <limits>

#ifdef WIN32
#include <windows.h>
#endif

#ifdef max 
#undef max
#endif
#ifdef min
#undef min
#endif

#include "EvalEnergyFull2Frame.h"
#include "EvalEnergyStore.h"
#ifdef _testLocalReplacement_
  #include "ReplacementProposals.h"
#endif
#ifdef _MotionChecking_
  #include "SensibleMotionCheck.h"
#endif
#include "DataSegments.h"

#include <algorithm>
#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

#include <ctime>

using namespace std;
using namespace Math;

// no data stuff and iterations fast debugging
//#define _shortCut_forTesting_
//#define __Energy_Check__ 
//#define __DEBUG_ENERGIES__
// sum over all views
#define _SumViewSmoothing_

template<typename Scalar>
void addUnariesBinaries( GraphCutContainer<Scalar>& q, const typename std::vector< Unary<Scalar> >& unaries, const typename std::vector< Binary<Scalar> >& binaries, const Scalar scale, int& non_sub)
{
  for (int k=0; k< unaries.size();k++)
    q.AddUnaryTerm(unaries[k].var, scale*( unaries[k].node[0] ) ,scale*( unaries[k].node[1] ) );
  for( int j=0; j < binaries.size(); j++ )
  {
    typename Datasegments<Scalar>::P4 ff = scale*binaries[j].edge;
    if (ff[0]+ff[3] > ff[1]+ff[2])
      non_sub++;
    if (ff[1] > 0 || ff[2] > 0  || ff[0] > 0 || ff[3] > 0)
      q.AddPairwiseTermLSAAux(binaries[j].segI, binaries[j].segJ, ff);
    //      q.AddPairwiseTerm(binaries[j].segI, binaries[j].segJ, scale*ff[0], scale*ff[1], scale*ff[2], scale*ff[3]);
  }
};

template<typename Scalar>
Scalar addUnariesBinariesWorstData( GraphCutContainer<Scalar>& q, const typename std::vector< Unary<Scalar> >& unaries, const typename std::vector< Binary<Scalar> >& binaries )
{
  Scalar worstDataScore (0);
  for (int k=0; k< unaries.size();k++)
    worstDataScore += max( unaries[k].node[0], unaries[k].node[1] );
  for (int k=0; k< binaries.size();k++)
    worstDataScore += max(max( max( binaries[k].edge[0], binaries[k].edge[1] ),binaries[k].edge[2]), binaries[k].edge[3] );
  return worstDataScore;
};

template<typename Scalar>
void addSmoothnessConstraints ( GraphCutContainer<Scalar>& q, const std::vector<Scalar>& f00, 
  const std::vector<Scalar>& f01, const std::vector<Scalar>& f10, const std::vector<Scalar>& f11, 
  const std::vector<int>& idk, const std::vector<int>& idl, Scalar scaleSmo, int& non_sub )
{
  for( int j=0; j < idk.size(); j++ )
  {
    if ( idk[j] >=0 && idl[j] >=0 )
    {
      if (f00[j]+f11[j] > f10[j]+f01[j])
        non_sub++;

	  typename Datasegments<Scalar>::P4 edge(scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);

      if (f10[j] > 0 || f01[j] > 0  || f00[j] > 0 || f11[j] > 0)
        q.AddPairwiseTermLSAAux(idk[j], idl[j], edge );
//        q.AddPairwiseTerm(idk[j], idl[j], scaleSmo*f00[j], scaleSmo*f01[j], scaleSmo*f10[j], scaleSmo*f11[j]);
    }
    else // onesided bandit
    {
      if ( idk[j] >= 0 )
        q.AddUnaryTerm(idk[j], scaleSmo*f00[j], scaleSmo*f10[j] );
      if ( idl[j] >= 0 )
        q.AddUnaryTerm(idl[j], scaleSmo*f00[j], scaleSmo*f01[j] );
    }
  }
}

int getNSegments( int*segImg, int N, int M )
{
  int maxSeg=0;
  for(int i=0;i<N*M;i++)
    maxSeg = std::max( maxSeg, segImg[i] );
  return maxSeg+1;
}

int getDsgId(int cam1, int time1, int nCams, int nTimeSteps) 
{
  return cam1 + time1*nCams;
}

// std::min(cam1,cam2) map to be unique - so ?? 
int getQuadId(int cam1, int time1, int cam2, int time2, int nCams, int nTimeSteps) 
{
  int c1,c2,t1,t2;
  if ( cam1 < cam2 || (cam1 == cam2 && time1< time2) )
  {c1 = cam1;c2=cam2;t1=time1;t2=time2;}
  else
  {c1 = cam2;c2=cam1;t1=time2;t2=time1;}
//  {c1 = cam1;c2=cam2;t1=time1;t2=time2;}// WTF

  return c1*nCams*nTimeSteps*nTimeSteps + c2*nTimeSteps*nTimeSteps + t1*nTimeSteps + t2;
}

Datasegments<double>::P4i getQuadIdInv(int id, int nCams, int nTimeSteps) 
{
  int c1,c2,t1,t2;

  t2  = id % nTimeSteps;
  id -= t2;id /= nTimeSteps;
  t1  = id % nTimeSteps;
  id -= t1;id /= nTimeSteps;
  c2  = id % nCams;
  id -= c2;id /= nCams;
  c1  = id % nCams;
  id -= c1;
  assert(id==0);
  return Datasegments<double>::P4i(c1,t1,c2,t2);
}

template<typename Scalar>
Scalar getEnergy ( std::vector<int>& current, int nSegments, const std::vector<Scalar>& dataScores, EvalEnergyFullFrame<Scalar>& evalSmooth, Scalar theta, Scalar lambda)//, std::vector<bool> oobS,  std::vector<bool> oobT)
{
  Scalar score(0.);
  for ( int j = 0; j< nSegments; j++ )
    score += dataScores[j];

  std::vector<int> allOn ( nSegments, 0);
  evalSmooth.compute_score_Fuse_local( current, current, allOn );
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
  std::vector<int> allOn ( nSegments, 0);
  std::vector<Scalar> _F00;std::vector<int> ids;
  evalSmooth.compute_score_Fuse_local( current, current, allOn, _F00, _F00, _F00, _F00, ids, ids, false );
  const std::vector<Scalar>& f00=evalSmooth.getF00();
  for ( int i = 0; i < f00.size(); i++)
    score+= f00[i] * lambda;

  return score;
}


// input: image dx, dy, normals, rotations, translations, p2d, K, mct Seg.Ids
void run_VCSF( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar; // double is givn - does float do something reasonable ?

  size_t N,M;

  if (nrhs < 6)
    mexErrMsgTxt("Must have at least 6 input arguments");

  Scalar ENERGY(0);
  int maxNodesUsed(0), maxEdgesUsed(0), nRuns(0), avNodesUsed(0), avEdgesUsed(0), non_sol(0), non_sub(0);

  Scalar gridSize = 16;
#ifdef _DEBUG
  Scalar gridSizeX = 6;//152/16; //104: 6
  Scalar gridSizeY = 4;//104/16; // 72: 4
#else
  Scalar gridSizeX = 9;//152/16; //104: 6
  Scalar gridSizeY = 6;//104/16; // 72: 4
#endif

  int nTimeSteps  = 2;
  Scalar thresh = 0.1;
  Scalar lambda = 0.02;
  Scalar theta  = 100.0;
  Scalar dispMax = 32.;
  Scalar pixJump = 1.0;
  Scalar rotJump(20.);   //   setRotJump (rotJump)
  Scalar depthJump(20.); //   setDepthJump (depthJump)
  Scalar occThresh(0.8);
  Scalar outOfBoundsThresh(0.8);
  Scalar motMax = 220; // the maximal motion in the scene
  Scalar gamma = 1.0;
  Scalar tempPotts = 0.25;//0.25;// less here in the first place, later more
  Scalar ot = 0.1;// originally 0.1
  Scalar doAuto =0.;

  typedef Datasegments<Scalar>::P4i camTimePair;
  std::vector<camTimePair> dataCams;

  int *prop2PlaneMap(NULL);
  Scalar* expCenters (NULL);

  Scalar* halfEdgeX(NULL);
  Scalar* halfEdgeY(NULL);
  Scalar* halfEdgeXY(NULL);
  Scalar* halfEdgeiXY(NULL);

  std::vector<Scalar*> pict;
  std::vector<Scalar*> centerst;
  std::vector<int*> segImg_t;
  std::vector<int*> solution_t;
  std::vector<const mxArray*> edgest;

  pict.push_back ( (Scalar*)  mxGetPr(prhs[0]) );
  pict.push_back ( (Scalar*)  mxGetPr(prhs[1]) );

  mwSize nDims = mxGetNumberOfDimensions(prhs[0]);
  assert(nDims == 3);
  const mwSize* imgDims = mxGetDimensions(prhs[0]);
  const int nCams   = imgDims[2];
  assert( mxGetNumberOfElements(prhs[0]) == mxGetNumberOfElements(prhs[1]) );
  M = imgDims[0];
  N = imgDims[1];

  Scalar* parameters = (Scalar*) mxGetPr(prhs[ 4]);// all parameters:
  int nParameters    = (size_t)  mxGetNumberOfElements(prhs[4]);

  if ( nParameters > 0) thresh            = parameters[0];
  if ( nParameters > 1) lambda            = parameters[1];
  if ( nParameters > 2) theta             = parameters[2];
  if ( nParameters > 3) rotJump           = parameters[3];
  if ( nParameters > 4) depthJump         = parameters[4];
  if ( nParameters > 5) dispMax           = parameters[5];
  if ( nParameters > 6) motMax            = parameters[6];
  if ( nParameters > 7) outOfBoundsThresh = parameters[7];//unsued -> autovector: per segment
  if ( nParameters > 8) occThresh         = parameters[8];//unsued -> autovector: per segment
  if ( nParameters > 9) gridSize          = parameters[9];
  if ( nParameters > 10) doAuto           = parameters[10];//doAuto=1;
  if ( nParameters > 11) tempPotts        = parameters[11];
  if ( nParameters > 12) ot               = parameters[12];
  if ( nParameters > 13) gridSizeX        = parameters[13];
  if ( nParameters > 14) gridSizeY        = parameters[14];


  segImg_t.push_back( (int*) mxGetPr(prhs[2])); 
  segImg_t.push_back( (int*) mxGetPr(prhs[3]));

  /// need to get # of 9segments here first
  solution_t.push_back ( (int*) mxGetPr(prhs[5]) );
  solution_t.push_back ( (int*) mxGetPr(prhs[6]) );

  std::vector< vector< vector<int> > > currentSolutions;

  std::vector< Scalar* > K_( nCams );
  assert( nCams*9 == mxGetNumberOfElements(prhs[7]) );
  for(int i=0;i<nCams;i++)
    K_[i] = &(((Scalar*)  mxGetPr(prhs[7]))[9*i]);

  // hm tough one, only other cams have this, so one per cam other then canonical
  // optional:

  std::vector< Scalar* > MC_( nCams-1 );
  std::vector< Scalar* > mc_( nCams-1 );
  assert( (nCams-1)*9 == mxGetNumberOfElements(prhs[8]) );
  assert( (nCams-1)*3 == mxGetNumberOfElements(prhs[9]) );
  for(int i=0;i<nCams-1;i++)
  {
    MC_[i] = &(((Scalar*)  mxGetPr(prhs[8]))[9*i]);
    mc_[i] = &(((Scalar*)  mxGetPr(prhs[9]))[3*i]);
  }

  Scalar* normals   =  (Scalar*) mxGetPr (prhs[10]);
  Scalar* rotations =  (Scalar*) mxGetPr (prhs[11]);
  centerst.push_back(  (Scalar*) mxGetPr (prhs[12]) );
  centerst.push_back(  (Scalar*) mxGetPr (prhs[13]) );

  // e.g. : which prop should be expanded at which expCenter, e.g. plane # 5 at center #77 -> prop2PlaneMap[77]=5
  if (nrhs > 14)
    prop2PlaneMap    =  (int*) mxGetPr(prhs[14]);
  // e.g. centers of segments of old solution (rougher solution)
  if (nrhs > 15)
    expCenters       =  (Scalar*) mxGetPr(prhs[15]);

  // multiple edges, assume it is a cell of edges-cells, so
  edgest.push_back ( (prhs[16]) );
  edgest.push_back ( (prhs[17]) );
  assert( mxGetNumberOfElements(prhs[16]) == nCams );
  assert( mxGetNumberOfElements(prhs[17]) == nCams );

  // here all autoscores are put in one, also for ? t0, nSegments t1, t-1 t2 ..
  Scalar* autoScore = NULL;
  autoScore = (Scalar*) mxGetPr(prhs[18]); // turn off if not needed - or set to all 1's or constant stuff

  if (nrhs > 19)
    halfEdgeX       =  (Scalar*) mxGetPr(prhs[19]);

  if (nrhs > 20)
    halfEdgeY       =  (Scalar*) mxGetPr(prhs[20]);

  if (nrhs > 21)
    halfEdgeXY      =  (Scalar*) mxGetPr(prhs[21]);

  if (nrhs > 22)
    halfEdgeiXY     =  (Scalar*) mxGetPr(prhs[22]);

  int* oobSegsR  = NULL;
  int* oobSegsLT = NULL;
  int* oobSegsRT = NULL;

  if (nrhs > 23)
    oobSegsR   =  (int*) mxGetPr(prhs[23]);

  if (nrhs > 24)
    oobSegsLT   =  (int*) mxGetPr(prhs[24]);

  if (nrhs > 25)
    oobSegsRT  =  (int*) mxGetPr(prhs[25]);// hijacked it is RRT now

  const int pastTimeNr = 26;// 26 regularly 126 to use only 2 time steps for speedup
  // past time step:
  Scalar *rotations_glob(NULL);
  if (nrhs > pastTimeNr)
    pict.push_back ( (Scalar*)  mxGetPr(prhs[pastTimeNr]) );
  if (nrhs > pastTimeNr+1)
    segImg_t.push_back( (int*) mxGetPr(prhs[pastTimeNr+1]) ); 
  if (nrhs > pastTimeNr+2)
    solution_t.push_back( (int*) mxGetPr(prhs[pastTimeNr+2]) );
  if (nrhs > pastTimeNr+3)
    centerst.push_back( (Scalar*) mxGetPr (prhs[pastTimeNr+3]) );
  if (nrhs > pastTimeNr+4)
    edgest.push_back ( prhs[pastTimeNr+4] );
  if (nrhs > pastTimeNr+5)
    rotations_glob      =  (Scalar*) mxGetPr(prhs[pastTimeNr+5]);

  if (nrhs > pastTimeNr+5)
    nTimeSteps =3;

  Scalar* rotations_glob_t1t2(NULL);
  if (nrhs > pastTimeNr+6)
    pict.push_back ( (Scalar*)  mxGetPr(prhs[pastTimeNr+6]) );
  if (nrhs > pastTimeNr+7)
    segImg_t.push_back( (int*) mxGetPr(prhs[pastTimeNr+7]) ); 
  if (nrhs > pastTimeNr+8)
    solution_t.push_back( (int*) mxGetPr(prhs[pastTimeNr+8]) );
  if (nrhs > pastTimeNr+9)
    centerst.push_back( (Scalar*) mxGetPr (prhs[pastTimeNr+9]) );
  if (nrhs > pastTimeNr+10)
    edgest.push_back ( prhs[pastTimeNr+10] );
  if (nrhs > pastTimeNr+11)
    rotations_glob_t1t2  =  (Scalar*) mxGetPr(prhs[pastTimeNr+11]);

  if (nrhs > pastTimeNr+11)
    nTimeSteps =4;

  std::vector< vector<Scalar*> > centers( nTimeSteps );
  std::vector< vector<int> >   nSegments( nTimeSteps );// per timestep, nCam # of Segments
  std::vector< vector<Scalar*> >    imgs( nTimeSteps );// per timestep nCam image
  std::vector< vector<int*> >     segImg( nTimeSteps );// per timestep nCam image

  int maxSegs=0;
  for( int j=0;j<nTimeSteps;j++ )
    for( int i=0;i<nCams;i++ )
    {
      nSegments[j].push_back( getNSegments( &segImg_t[j][N*M*i], M,N ) );
      segImg[j].push_back ( &segImg_t[j][N*M*i] );
      imgs[j].push_back (   &pict[j][N*M*i] );
      maxSegs = std::max (nSegments[j][i], maxSegs );
    }

    // 3d coordiantes to rotate around - if all 0's ..
    // max of nSegments guess:
    std::vector<Scalar> origins(3*maxSegs, 0); // a hack for origin centric rotations

    for( int k=0;k<nTimeSteps;k++ )
    {
      int temp=0;
      currentSolutions.push_back( std::vector< vector<int> > (nCams) );
      for(int i=0;i<nCams;i++ )
      {
        currentSolutions[k][i].resize( nSegments[k][i] );

        for(int j=0;j<nSegments[0][i];j++)
          (currentSolutions[k][i])[j] = (solution_t[k][temp + j]);
        temp += nSegments[k][i];
      }
    }

    for( int j=0;j<nTimeSteps;j++ )
    {
      int temp = 0;
      for(int i=0;i<nCams;i++ )
      {
        centers[j].push_back ( &(centerst[j][temp*3]) );
        temp += nSegments[j][i];
      }
    }

    // if not always take 1st?
    std::vector< vector< const mxArray*> > edges(nTimeSteps);
    for( int j=0;j<nTimeSteps;j++ )
      for(int i=0;i<nCams;i++)
        edges[j].push_back ( mxGetCell(edgest[j], i) );

      std::vector< vector<Scalar*> > autoScores( nTimeSteps );// start 2 timesteps - can be more

      int temp=0;
      for( int j=0;j<nTimeSteps;j++ )
      {
        autoScores[j].resize(nCams, NULL);
        if (autoScore != NULL )
          for(int i=0;i<nCams;i++)
          {
            if ( temp+nSegments[j][i] <= mxGetNumberOfElements(prhs[18]) )
              (autoScores[j])[i] = &( autoScore[temp] );
            temp += nSegments[j][i];
          }
      }

      std::vector< vector<Scalar*> > heX( nTimeSteps );// start 2 timesteps - can be more
      std::vector< vector<Scalar*> > heY( nTimeSteps );// start 2 timesteps - can be more
      std::vector< vector<Scalar*> > heXY( nTimeSteps );// start 2 timesteps - can be more
      std::vector< vector<Scalar*> > heYX( nTimeSteps );// start 2 timesteps - can be more

      for (int i=0;i< nTimeSteps; i++)
      {
        heX[i].resize(nCams, NULL);
        heY[i].resize(nCams, NULL);
        heXY[i].resize(nCams, NULL);
        heYX[i].resize(nCams, NULL);
      }

      assert ( (nTimeSteps*nCams)*(N-1)*M     <= mxGetNumberOfElements(prhs[19]) );
      assert ( (nTimeSteps*nCams)*N*(M-1)     <= mxGetNumberOfElements(prhs[20]) );
      assert ( (nTimeSteps*nCams)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[21]) );
      assert ( (nTimeSteps*nCams)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[22]) );

      for (int j=0;j< nTimeSteps; j++)
      {
        for(int i=0;i<nCams;i++)
        {
          if ( (i+1)*(N-1)*M <= mxGetNumberOfElements(prhs[19]) )
            (heX[j])[i] = &( halfEdgeX[(i+j*nCams)*(N-1)*M] );
          if ( (i+1)*(M-1)*N <= mxGetNumberOfElements(prhs[20]) )
            (heY[j])[i] = &( halfEdgeY[(i+j*nCams)*(M-1)*N] );
          if ( (i+1)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[21]) )
            (heXY[j])[i] = &( halfEdgeXY[(i+j*nCams)*(N-1)*(M-1)] );
          if ( (i+1)*(N-1)*(M-1) <= mxGetNumberOfElements(prhs[22]) )
            (heYX[j])[i] = &( halfEdgeiXY[(i+j*nCams)*(N-1)*(M-1)] );
        }
      }

      // a grid if nothing else specified:
      // input: layout for one! timestep; copy this layout for all time steps
      // and add connections from t to t+1
      for (int j=1;j<nCams;j++)
        dataCams.push_back( camTimePair(0,0, j,0) );
      for (int j=0;j<nCams;j++)
        dataCams.push_back( camTimePair(j,0, j,1) );
      if (nTimeSteps > 2)
        for (int j=0;j<nCams;j++)
          dataCams.push_back( camTimePair(j,0, j,2) );
      if (nTimeSteps > 3)
        for (int j=0;j<nCams;j++)
          dataCams.push_back( camTimePair(j,1, j,3) );
      for(int i=1;i<nTimeSteps;i++)
        for (int j=1;j<nCams;j++)
          dataCams.push_back( camTimePair(0,i, j,i) );

      int nNormals = max( mxGetM(prhs[10]), mxGetN(prhs[10]) );
      int nStep    = min( mxGetM(prhs[10]), mxGetN(prhs[10]) );
      // define the number of expansions here:
#ifdef _DEBUG
      int nProposalExpansions = min(200, (int) max( mxGetM(prhs[14]), mxGetN(prhs[14]) ));// Fast
#else
      int nProposalExpansions = max( mxGetM(prhs[14]), mxGetN(prhs[14]) );
#endif
      //////////////////////////

//      const int decreaser = 1;
//      int expPatches = 8/decreaser; // expand patches by (warp bbox, expand it, test all segments in the box in 2nd view)
//      int patchX = 104/decreaser;// !! if too many segments storing the stuff is too expensive w.r.t. memory
//      int patchY = 72/decreaser;

  	  int expPatches (Scalar(gridSize)/2.);
   	  int patchX = gridSize * (gridSizeX+0.5);
	    int patchY = gridSize * (gridSizeY+0.5);

      std::srand ( unsigned ( std::time(0) ) );
      //std::srand ( 555 ); // for debugging
      //////////////////////////////////////////////
      printf("rotJump %f \n", rotJump);
      printf("depthJump %f \n", depthJump);
      printf("dispMax %f \n", dispMax);
      printf("motMax %f \n", motMax);
      printf("occThresh %f \n", occThresh);
      printf("outOfBoundsThresh %f \n", outOfBoundsThresh);
      printf("data thresh %f \n", thresh);
      printf("motion smooth theta %f \n", theta);
      printf("general smoothness lambda %f \n", lambda);
      printf("gridSize: %.1f, gridX %.1f  gridY %.1f -> Patches %dx%d\n", gridSize, gridSizeX, gridSizeY, patchX, patchY );
      printf("nSegments %d, nSegmentsR0 %d, nNormals %d nProps: %5d\n", nSegments[0][0], nSegments[0][1], nNormals, nProposalExpansions);
      printf("doAuto %.1f, ot:%.2f, tempPotts:%.2f\n", doAuto, ot, tempPotts );
      //mexEvalString("drawnow");

      /// from iccv13 model this is NOT used in VC thus all input is supposed to be false
      std::vector<bool> oobOutsideAll(M*N, false);
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


        std::vector< std::vector<bool>*> oobVex(nCams*nCams*nTimeSteps*nTimeSteps, &oobOutsideAll);
        std::vector< int > dsgId   (nCams*nCams*nTimeSteps*nTimeSteps, -1);
        std::vector< Datasegments<Scalar>::P4i > dsgIdInv;
        std::vector< bool > dataOn (nCams*nCams*nTimeSteps*nTimeSteps, false );

        oobVex[ getQuadId( 0, 0, 1, 0, nCams, nTimeSteps) ] = &oobOutsideLR;
        oobVex[ getQuadId( 1, 0, 1, 1, nCams, nTimeSteps) ] = &oobOutsideLRT;
        oobVex[ getQuadId( 0, 0, 0, 1, nCams, nTimeSteps) ] = &oobOutsideLLT;

        for( int i=0;i<dataCams.size();i++)
        {
          int wtf = getQuadId( dataCams[i][0], dataCams[i][1], dataCams[i][2], dataCams[i][3], nCams, nTimeSteps);
          dataOn [ getQuadId( dataCams[i][0], dataCams[i][1], dataCams[i][2], dataCams[i][3], nCams, nTimeSteps) ] = true;
        }
        // must do: first initHom .. then warps based on this: 
        // later expand -> vars, then expandregion based on this

        mwSize dims[2];dims[0] = N;dims[1] = M;

        ///////////////////////////////////////////////////////
        EvalEnergyStore<Scalar> enStore;
        enStore.prepareNormals( normals, nNormals, nStep );
        enStore.setRotTra( rotations, nNormals );
#ifdef _testLocalReplacement_
        // ?? can i do a repeat whole thing i.e. fit new props with recomputation of whole data, etc ? 
        // homos not needed since copied - number ? 
        /// new LocalReplacement Strategy simplified here
        LocalReplacement<Scalar> lrp(enStore, segImg[0][0], nSegments[0][0], (int*) (&(currentSolutions[0][0][0])), centers[0][0], edges[0][0]);//, K_[0]);
        prop2PlaneMap       = lrp.getPropMap();
        expCenters          = lrp.getExpCenters();
        nProposalExpansions = lrp.getNProposals();
        nNormals            = (enStore.getNormals())->size();
 
#ifdef _DEBUG
       patchX =1;patchY =1;
       nProposalExpansions = min(200, nProposalExpansions);
#endif

#endif
       /*
        // these must be done to address consistency issues between propid, id stored in the segmentation and solution for the prop2plane map
        std::vector<int> prop2Plane( nProposalExpansions, 0 );
        for ( int i=0; i < nProposalExpansions; i++ )
          if ( prop2PlaneMap != NULL)
            prop2Plane[i] = prop2PlaneMap[i];
          else
            prop2Plane[i] = i;


        // for each proposal, center its (new) box after joining similar props together
        std::vector<dataElem<Scalar>::P2i> patchXY( nProposalExpansions, dataElem<Scalar>::P2i(patchX,patchY) );
        std::vector<Scalar>  propCenters;
*/
        ////////////////////////////////////////////
        // all homos to all involved other views !

        // indeed like this: star config a homo from canonical to all other views:
        std::vector< genHomoG<Scalar> > homos;
        for (int j=0; j< nTimeSteps; j++)
        {
          for(int i=0;i<nCams;i++)
          {
            if (i==0 && j==0) {homos.push_back( genHomoG<Scalar>() ); continue;}

            if (j==2 && rotations_glob != NULL)
              enStore.setGlobRotTra( rotations_glob );
            if (j==3 && rotations_glob_t1t2  != NULL)
              enStore.setGlobRotTra( rotations_glob_t1t2 );

            int id = j*nCams+i;
            if (i==0)  // cases: j!=0 i==0(timestep!=0, camera canonical)
              homos.push_back( genHomoG<Scalar> ( K_[0] ) );
            else       // cases: j!=0 i!=0(timestep!=0, camera not canonical)
              homos.push_back( genHomoG<Scalar> ( K_[0], MC_[i-1], mc_[i-1], K_[i]) );

            homos[id].setNormals( enStore.getNormals() );
            if (j==0) continue;
            homos[id].setRotTra(  enStore.getTra(), enStore.getRot() );
            homos[id].setCenters( &(origins[0]) );// could put all 0's here -> Rt w.r.t. origin

            if (j==2)
              homos[id].do_inverse_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
            if (j==3)
              homos[id].do_double_motion(true, enStore.getGlobRot(), enStore.getGlobTra());
          }
        }
        // ? start with those receiving the homos ? then 
        // use a helper vector: a> warp==#'s b> use for data term, also which needs info of which
        std::vector < Datasegments<Scalar>* > dsgX;//(nCams*nTimeSteps);
        // store ids of initialisers and data dsg's !

        std::vector < int > dsg_init;
        std::vector < int > dsg_data;

        for (int j=0; j< nTimeSteps; j++)
        {
          for(int i=0;i<nCams;i++)
          {
            if (i==0 && j==0) continue;

            if ( dataOn [ getQuadId( 0, 0, i, j, nCams, nTimeSteps) ] )
              dsg_data.push_back(dsgX.size());
            dsg_init.push_back(dsgX.size());
            dsgId[ getQuadId( 0,0,i,j, nCams, nTimeSteps) ] = dsgX.size();

            // pushback copies which is hard if pointer are involved - sigh
            dsgX.push_back( new Datasegments<Scalar> (N, M, thresh, outOfBoundsThresh, imgs[0][0], imgs[j][i], segImg[0][0], nSegments[0][0], segImg[j][i], nSegments[j][i], (*oobVex[ getQuadId( 0, 0, i, j, nCams, nTimeSteps)] ),  oobOutsideAll) );

            dsgX.back()->setGridSize( gridSize );
            dsgX.back()->setTimeCam( Datasegments<Scalar>::P4i(0,0,i,j) );
            dsgX.back()->setCentersKs( K_[0], K_[i], centers[0][0],  centers[j][i] );
            dsgX.back()->setOccPen ( occThresh );

            if (j==0)// time step 0
            {
               if (i >0 && mc_[i-1][0] > 0)
                 dsgX.back()->setMaxDisp( -dispMax );
               else
                 dsgX.back()->setMaxDisp(  dispMax );
            }
            else
            {
              if (i==0) // cam 0
                dsgX.back()->setMaxMot( fabs((double)j)*motMax );// never called anyway
              else
                dsgX.back()->setMaxMot( fabs((double)j)*motMax+dispMax );// never called anyway
            }

            if ( autoScores[0][0] != NULL && autoScores[j][i] != NULL )
              dsgX.back()->set_AutoScores ( autoScores[0][0], autoScores[j][i]);
            dsgX.back()->init();
          }
        }

        for (int k = 0; k< dataCams.size(); k++ )
        {
          int c1 = dataCams[k][0];// 1st is time though ?
          int c2 = dataCams[k][2];
          int t1 = dataCams[k][1];
          int t2 = dataCams[k][3];

          if ( dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] >=0)
            continue;
          dsgId[ getQuadId( c1,t1,c2,t2, nCams, nTimeSteps) ] = dsgX.size();

          dsg_data.push_back(dsgX.size());  
          dsgX.push_back( new Datasegments<Scalar> (N, M, thresh, outOfBoundsThresh, imgs[t1][c1], imgs[t2][c2], segImg[t1][c1], nSegments[t1][c1], segImg[t2][c2], nSegments[t2][c2], (*oobVex[ getQuadId( c1, t1, c2, t2, nCams, nTimeSteps)] ),  oobOutsideAll) );

          dsgX.back()->setGridSize( gridSize );
          dsgX.back()->setCentersKs( K_[c1], K_[c2], centers[t1][c1],  centers[t2][c2] );
          dsgX.back()->setTimeCam( Datasegments<Scalar>::P4i(c1,t1,c2,t2) );
          dsgX.back()->setOccPen ( occThresh );

          if (t1==t2)
          {
            Scalar tx_c1 = c1>0 ? mc_[c1-1][0] : 0;
            Scalar tx_c2 = c2>0 ? mc_[c2-1][0] : 0;
            if ( tx_c2 < tx_c1 )
              dsgX.back()->setMaxDisp(  dispMax ); // TODO: check
            else
              dsgX.back()->setMaxDisp( -dispMax ); // 
          }
          else
          {
            if (c1==c2)
              dsgX.back()->setMaxMot( motMax );// 
            else
              dsgX.back()->setMaxMot( motMax+dispMax );// never called anyway by construction of data term
          }

          if ( autoScores[t1][c1] != NULL &&  autoScores[t2][c2] != NULL )
            dsgX.back()->set_AutoScores ( autoScores[t1][c1], autoScores[t2][c2]);
          dsgX.back()->init();
        }

//////////////////////////////////////////////// repeat here ? 
#ifdef _testLocalReplacement_
        const int allLocRepIts = __locRepsReplace__;int locRepIts = allLocRepIts;
        while(locRepIts>0)
        {
          locRepIts--;
#endif

        // these must be done to address consistency issues between propid, id stored in the segmentation and solution for the prop2plane map
        std::vector<int> prop2Plane( nProposalExpansions, 0 );
        for ( int i=0; i < nProposalExpansions; i++ )
          if ( prop2PlaneMap != NULL)
            prop2Plane[i] = prop2PlaneMap[i];
          else
            prop2Plane[i] = i;


        // for each proposal, center its (new) box after joining similar props together
        std::vector<dataElem<Scalar>::P2i> patchXY( nProposalExpansions, dataElem<Scalar>::P2i(patchX,patchY) );
        std::vector<Scalar>  propCenters;


        for (int i = 0;i < nNormals; i++)
        {
          for (int k = 0;k < dsg_init.size(); k++)
            dsgX[dsg_init[k]]->initHoms(    i, &homos[dsg_init[k]+1],   enStore.getNormals() );
          for (int k = 0;k < dsg_data.size(); k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
            int id1= dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
            int id2= dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
            dsgX[dsg_data[k]]->initHoms( i, dsgX[id1]->getviHoms(), dsgX[id2]->getvHoms(), dsgX[id1]->getpiNoms(), dsgX[id2]->getpiNoms() );
          }
        }

        // fast combination of segment expansion areas speeds up, same energy ! -- could be even better
        // here setExpCenters: speed up, combining identical proposals to a 'super-proposal' in terms of expansion area
        dsgX[0]->reorganizeProposals ( prop2Plane, nProposalExpansions, expCenters, propCenters, patchX, patchY, patchXY );
        expCenters = (Scalar*) (&(propCenters[0]));

        if (expCenters != NULL)
        {
          for (int k=0;k<dsgX.size();k++)
            dsgX[k]->setExpCenters (expCenters, nProposalExpansions);
#ifdef _writeEnergies_
          printf("set extra centers \n");
#endif
        }
        // proposals can be a segment/proposal id pair now - fit at rough resolution
        // however which normal and which position do not correlate
        // proposal number i could also be normal j
        // i store the id of the normal/rot in the segimg
        // so here the parameter is not i but the normal's id
        // stored should be then the prop-id at the point, not the number 
        // here parallel ?

#ifdef _NO_OPENMP
        int nProcessors = 1;
#else
        //parallel over i
        // request gs and gsw 1/2 for n processors release later proc-# in function call(s) as last parameter
        int nProcessors = min (_max_num_threads_used_, omp_get_max_threads());
#ifdef _DEBUG
        nProcessors  = 1;
#endif
        omp_set_num_threads( nProcessors  );
        printf("Threads used: %d\n", nProcessors );
#endif
        // allow for parallel expansion, keeping track internally what areas are already covered
        for (int k=0;k<dsgX.size();k++)
          dsgX[k]->requestNElements( nProcessors, nProposalExpansions );

        int maxVariables =0;
        std::vector<int> expPropVars(nProposalExpansions, 0);
        std::clock_t startP(std::clock());
#ifndef  _shortCut_forTesting_
#pragma omp parallel for
        for (int i = 0;i < nProposalExpansions; i++)
        {
          int varNr(0);
#ifdef _NO_OPENMP
          const int pid = 0;
#else
          const int pid = omp_get_thread_num();
#endif

          for (int k = 0;k < dsg_init.size(); k++)
            varNr = dsgX[ dsg_init[k] ]->initNewSeg( i, prop2Plane[i], patchXY[i][0], patchXY[i][1], &homos[dsg_init[k]+1], expPatches, enStore.getNormals(), varNr ,0, pid );

          expPropVars[i] = varNr + dsgX[ dsg_init[0] ]->getElemVarsFirst( pid );
          maxVariables = max ( expPropVars[i], maxVariables );

          for (int k = 0;k < dsg_data.size(); k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            if (ct_ct[0] == 0 && ct_ct[1] ==0 ) continue;
            int id1= dsgId[ getQuadId( 0,0, ct_ct[0],ct_ct[1], nCams, nTimeSteps) ];
            int id2= dsgId[ getQuadId( 0,0, ct_ct[2],ct_ct[3], nCams, nTimeSteps) ];
            dsgX[dsg_data[k]]->initNewSegFromView( i, prop2Plane[i], dsgX[id1]->getBox2(pid), dsgX[id2]->getBox2(pid), expPatches, dsgX[id1]->getElement(pid), dsgX[id2]->getElement(pid), pid );
          }
          for (int k = 0;k < dsg_data.size(); k++)
            dsgX[dsg_data[k]]->createWarps( pid, i );
        }

        for (int k = 0;k < dsgX.size(); k++)
          dsgX[ k ]->releaseNElements( );
#endif

        std::clock_t stopP(std::clock());
        printf("Preprocess Data took %.2f\n", double(stopP-startP)/CLOCKS_PER_SEC );
        //  mexEvalString("drawnow");

        std::vector< std::vector<int>* > currentSolutionsNr;

        // for smoothing:
        for(int i=0;i< nCams*nTimeSteps;i++)
          currentSolutionsNr.push_back( &currentSolutions[i/nCams][i%nCams] );

        std::vector<int> order(nProposalExpansions, 0);
        for ( int i = 0;i < nProposalExpansions;i++ )
          order[i]=i;

        //  dsg.sortInitNewPixel ( order );//sorted by depth at center pixel


#ifdef _testLocalReplacement_
        if (locRepIts == allLocRepIts-1) // only once -- appears ok
#endif
#pragma omp parallel for
          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            dsgX[ dsg_data[k] ]->buildFromCurrentSolution(  currentSolutions[ct_ct[1]][ct_ct[0]], currentSolutions[ct_ct[3]][ct_ct[2]] );
          }

        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        // now smoothness cost setup:

        enStore.flipNormals();// stupid internal assumption - but to hard to change everywhere
        std::vector< EvalEnergyFullFrame<Scalar> > evalSmooth(nCams*nTimeSteps);

        for(int k=0;k<nTimeSteps;k++)
          for(int j=0;j<nCams;j++)
          {
            evalSmooth[k*nCams+j].setRotJump  (rotJump);
            evalSmooth[k*nCams+j].setDepthJump(depthJump);
            evalSmooth[k*nCams+j].setGamma ( gamma );
            evalSmooth[k*nCams+j].setRotWeight( theta );
            // jump at 0.5 pixel -> needs to set jump smaller
            if(k==0)
              evalSmooth[k*nCams+j].set2dMotionMatrix ( K_[0], K_[1], MC_[std::max(0,j-1)], mc_[std::max(0,j-1)], pixJump ); 
            else
              evalSmooth[k*nCams+j].set2dMotionMatrix_inv ( K_[0], K_[1], MC_[std::max(0,j-1)], mc_[std::max(0,j-1)], pixJump );
            evalSmooth[k*nCams+j].setRotTra( enStore.getRot(), enStore.getTra() );

            if( k==0 && j==0 )
              evalSmooth[k*nCams+j].setNormals( enStore.getNormals() );
            else
            {
              int id = dsgId[ getQuadId( 0,0,j,k, nCams, nTimeSteps) ];
              evalSmooth[k*nCams+j].setNormals( dsgX[id]->getvNormals() );
            }

            evalSmooth[k*nCams+j].setHalfEdges(heX[k][j], heY[k][j]);
            if ( halfEdgeiXY != NULL)
              evalSmooth[k*nCams+j].setCrossHalfEdges(heXY[k][j], heYX[k][j]);
            if ( halfEdgeiXY != NULL)
              evalSmooth[k*nCams+j].setCrossHalfEdges(heXY[k][j], heYX[k][j]);

            evalSmooth[k*nCams+j].prepareOwnWeights( nSegments[k][j], edges[k][j], centers[k][j], segImg[k][j], N, M );
          }

          for(int i=0;i<dsg_data.size(); i++)
          {
            dsgX[ dsg_data[i] ]->setOccThresh( ot );
            dsgX[ dsg_data[i] ]->setPotts(tempPotts);
          }

//          View-consistency off :		  
//          for(int i=0;i<dsg_data.size(); i++)
//            dsgX[ dsg_data[i] ]->setSimple(true);

#ifdef _MotionChecking_
          MotionCheck<Scalar> mChecker( NULL, nSegments[0][0], nNormals );
//          mChecker.setEpsilon(1); // no help to 1?
          mChecker.setPenalty(8. * thresh); // in dependence on data
          mChecker.setNormals( enStore.getNormals() );
          mChecker.setRotTra( enStore.getRot(), enStore.getTra() );
          mChecker.setWidthHeight(N,M);
          mChecker.setK( K_[0] );
#endif

          Scalar OLD_ENERGY = std::numeric_limits<Scalar>::max();
//          Scalar ENERGY(0);

          int numEdges =0;
          for (int e=0;e<evalSmooth.size(); e++ )
            numEdges += evalSmooth[e].edges_num( );

          int maxValue = 1 << 29;
//          int non_sub = 0;
//          int non_sol = 0;

          std::vector<int> noSolutions;noSolutions.reserve( nSegments[0][0]+nSegments[0][1]+nSegments[1][0]+nSegments[1][1] );

          std::vector< GraphCutContainer<Scalar> > gs ( 4 );
          for (int i=0;i<gs.size();i++)
            gs[i].create( maxVariables, numEdges+4*maxVariables );

#ifdef _testLocalReplacement_
        if (locRepIts == allLocRepIts-1) // only once -- appears ok
        {
#endif
          std::vector<Scalar> dataEnergies(dsg_data.size(),0);
          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            Scalar endE = dsgX[ dsg_data[k] ]->getEnergy();
            dataEnergies[k] = endE;
#ifdef _writeEnergies_
            printf("Start Energy:%.2f time:%d cam:%d <-> time:%d cam:%d\n", endE, ct_ct[1], ct_ct[0], ct_ct[3], ct_ct[2]);
#endif
          }
          Scalar sumDataEMulti = std::accumulate(dataEnergies.begin(),dataEnergies.end(),0.);

          std::vector<Scalar> smoothEnergies(nCams*nTimeSteps,0);
#ifdef _SumViewSmoothing_
          for (int k=0; k<evalSmooth.size();k++)
#else
          for (int k=0; k<1;k++)
#endif
            smoothEnergies[k] = getSmoothEnergy ( currentSolutions[k/nCams][k%nCams], nSegments[k/nCams][k%nCams], evalSmooth[k], theta, lambda);
          Scalar sumSmoothEMulti = std::accumulate(smoothEnergies.begin(),smoothEnergies.end(),0.);

#ifdef _writeEnergies_
          for(int kk= 0; kk< smoothEnergies.size(); kk++ )
            printf("cam/time (%d,%d) smooth Score: %.2f\n", kk/nCams, kk%nCams, smoothEnergies[kk]);
#endif
          ENERGY = sumDataEMulti + sumSmoothEMulti;
#ifdef _writeEnergies_
          printf("------------ Start Energy is %.2f --------------\n", ENERGY);
          printf("Start Energies (full=data+smooth) %.2f = ([%.2f],[%.2f,%.2f,%.2f,%.2f])\n", ENERGY, 
            sumDataEMulti, smoothEnergies[0], smoothEnergies[1], smoothEnergies[2], smoothEnergies[3] );
#endif
#ifdef _testLocalReplacement_
        }
#endif

          std::clock_t fullStart(std::clock());
          // fuse all proposals in a row

          std::clock_t startP3(std::clock());
          // segments where mCheck fails
//          int badCount(0), goodCount(0);
//          int maxNodesUsed(0), maxEdgesUsed(0), nRuns(0), avNodesUsed(0), avEdgesUsed(0);
          Scalar ENERGYOLD = ENERGY;

#ifdef _shortCut_forTesting_
          while( false )
#else
#ifndef _testLocalReplacement_
          while( ENERGY+__Small_Improvement__ < OLD_ENERGY )
#else
          for(int LocRepRuns=2;LocRepRuns>0;LocRepRuns--)
#endif
#endif
          {
            // full sweep over all proposal:
            random_shuffle( order.begin(), order.end() );

            // which prposals are yet to be processed
            std::vector<bool> processed (order.size(), false);
            int sstop(0);//finished?

#pragma omp parallel
            {
              std::vector< vector< Binary<Scalar> > > allBinaries( dsg_data.size() );
              std::vector< vector< Unary<Scalar> >  > allUnaries(  dsg_data.size() );

              while (!sstop)
              {
                std::vector< std::vector<int> > trialVecs( nCams*nTimeSteps );
                std::vector< std::vector<int>* > seg2vars( nCams*nTimeSteps, NULL );

                //#ifndef _DEBUG
#ifdef _NO_OPENMP
                const int pid = 0;
#else
                const int pid = omp_get_thread_num();
#endif
                int propId(0);
                int nVars(0), ord(0);
                bool okToRun (false), endIt (true);

#pragma omp critical
                {
#ifndef _DEBUG
                  for(; ord < processed.size(); ord++)
#else
                  for(; ord < std::min( 50, int(processed.size())); ord++)
#endif
                  {
                    if ( !processed[ord] )
                    {
                      endIt = false;
                      propId = order[ord];// propId - the id of the moving plane - this is in general not related to the center !

                      okToRun=true;
                      for( int k=0;k<dsg_init.size();k++)
                        okToRun = okToRun && dsgX[dsg_init[k]]->testNewProp( propId );
                      if (okToRun) break;// okToRun == true now
                    }
                  }

                  // end it now
                  if (endIt) // all processed
                  {
                    sstop = 1;
#pragma omp flush(sstop)
                  }

                  if (!okToRun)//if ( pid >= numThreads)
                  {
                    // should occur seldom - otherwise performance sucks!
                    //              printf("TOTAL nonsense error %d,%d \n", pid, nProcessors);
                    ord = order.size();
#ifdef WIN32
                    Sleep( 20 );
#endif
                    nVars =0;
                  }

                  if (ord < order.size() )
                  {
                    processed[ord] =true;
                    nVars = expPropVars[ propId ];
                  }

                  // do something now in case
                  if (nVars > 0 )
                  {
                    propId = order[ord];// states which 'pseudo' segment is to be expanded (with its solution)
                    // printf("perform run pid: %d, prop %d\n", pid, propId);//mexEvalString("drawnow");
                    // moving plane: propId. located at segment expcenter[propid]

                    gs[pid].setup(nVars);

                    for( int k=0;k<dsg_data.size();k++)
                      dsgX[dsg_data[k]]->expandRegion( propId );

                    // fill unaries, binaries by copy !
                    // end critical section

                    for(int k=0; k<dsg_data.size(); k++ )
                    {
                      allUnaries  [k] = dsgX[ dsg_data[k] ]->getUnaries();
                      allBinaries [k] = dsgX[ dsg_data[k] ]->getBinaries();
                    }

                    for (int k=0; k<dsgX.size();k++)
                      dsgX[k]->storeReleaseProp( pid, propId);
                  } // nVars>0

                }//critical section

                if (nVars > 0 )
                {
                  // seg2vars have to follow the order of smoothness
                  // order of dsg_init is time, cam, so time 0, cam 1,2,3 then time 1, cam 1,2,3
                  // so camid + timeId * nCams - keep that order here
                  for (int k=0;k<dsg_data.size();k++)
                  {
                    Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
                    if (seg2vars[ct_ct[0] + ct_ct[1] * nCams] == NULL)
                    {
                      trialVecs[ ct_ct[0] + ct_ct[1] * nCams ] = std::vector<int> ( nSegments[ct_ct[1]][ct_ct[0]], prop2Plane[propId] );
                      seg2vars[ct_ct[0] + ct_ct[1] * nCams ] = ( &(dsgX[ dsg_data[k] ]->getSegmentsInvolvedFW( propId )) );
                    }
                    if (seg2vars[ct_ct[2] + ct_ct[3] * nCams] == NULL)
                    {
                      trialVecs[ ct_ct[2] + ct_ct[3] * nCams ] = std::vector<int> ( nSegments[ct_ct[3]][ct_ct[2]], prop2Plane[propId] );
                      seg2vars[ct_ct[2] + ct_ct[3] * nCams ] = ( &(dsgX[ dsg_data[k] ]->getSegmentsInvolvedBW( propId )) );
                    }
                  }

                  Scalar worstDataScore(0);

                  for(int k=0; k<allUnaries.size(); k++ )
                    worstDataScore += addUnariesBinariesWorstData( gs[pid], allUnaries[k], allBinaries[k] );

#ifdef __DEBUG_ENERGIES__
                  for ( int k = 0;k < nSegments[0][0];k++ )
                    oldSolution[k]=currentSolution[k];
#endif
                  /////////////////////
                  std::vector<Scalar> worstSmoothScores(evalSmooth.size(), 0);
                  Scalar scale(0),  scaleSmo(0);

                  // stored in f00, etc. with critical super slow and wrong ?!?!?!
                  // not parallel at ll: correct
                  Scalar worstSmoothScore(0);

                  std::vector<Scalar> f00,f01,f10,f11;
                  std::vector<int>    idk,idj;
                  f00.reserve( 6*nVars );
                  f01.reserve( 6*nVars );
                  f10.reserve( 6*nVars );
                  f11.reserve( 6*nVars );
                  idk.reserve( 6*nVars );
                  idj.reserve( 6*nVars );

                  for(int kk = 0;kk< evalSmooth.size(); kk++)
                    worstSmoothScore += evalSmooth[kk].compute_score_Fuse_local( *(currentSolutionsNr[kk]), trialVecs[kk], *(seg2vars[kk]), f00,f01,f10,f11,idk,idj );

                  if (f00.size() > 6*nVars )
                    printf("Warning: f00.size() %d > 6*nVars:%d\n", f00.size(), nVars*6);

                  scale = (Scalar)(maxValue) / (lambda * (worstSmoothScore) + worstDataScore);
                  scaleSmo = scale * lambda;

                  addSmoothnessConstraints ( gs[pid], f00, f01, f10, f11, idk, idj, scaleSmo, non_sub );

                  ///////////////
#ifdef _MotionChecking_
                  //      std::clock_t startM(std::clock());
                  std::vector<Scalar> motionCheckScores1, motionCheckScores2;
                  mChecker.computeGeometryConsistencyLocal( segImg[0][0], nSegments[0][0],*(currentSolutionsNr[0]), *(seg2vars[0]), MotionCheck<Scalar>::P3(mc_[0]), dispMax, motionCheckScores1 );
                  mChecker.computeGeometryConsistencyLocal( segImg[0][0], nSegments[0][0],           trialVecs[0],  *(seg2vars[0]), MotionCheck<Scalar>::P3(mc_[0]), dispMax, motionCheckScores2 );

                  //      std::clock_t stopM(std::clock());
                  //    printf("MotionChecking took %f\n", double(stopM-startM)/CLOCKS_PER_SEC);
                  ///////////////////////

                  for (int k=0; k< (seg2vars[0])->size(); k++)
                    if ( (*(seg2vars[0]))[k] > -1)
                      gs[pid].AddUnaryTerm((*(seg2vars[0]))[k], scale*motionCheckScores1[k], scale*motionCheckScores2[k] );
#endif

                  for(int k=0; k<allUnaries.size(); k++ )
                    addUnariesBinaries( gs[pid], allUnaries[k], allBinaries[k], scale, non_sub );
                  ////////////////////////
                  int nEdges = f00.size();
                  for(int k=0;k< allBinaries.size();k++)
                    nEdges += allBinaries[k].size();

                  maxNodesUsed = max( maxNodesUsed, nVars);
                  maxEdgesUsed = max( maxEdgesUsed, nEdges);
                  nRuns ++;
                  avNodesUsed += nVars;
                  avEdgesUsed += nEdges;

                }//nvars >0

                if(nVars>0)
                {
                  gs[pid].finalize();
                  gs[pid].solve();

#ifdef _qpboVersion_				  
                  // NEW: improve moves:
                  gs[pid].resolve(false, false);
#endif
                  std::vector<int> solution01( nVars, 0 );
                  for (int j = 0; j < nVars;j++)
                  {
                    int x = gs[pid].GetLabel(j);
                    if (x<0)
                      non_sol++;

                    if (x==1)
                      solution01[j] = 1;
                  }

                  gs[pid].reset();

                  // currentSolutions do not overlap, so parallel should be ok

                  for (int k=0; k<dsg_data.size();k++)
                  {
                    Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
                    dsgX[dsg_data[k]]->updateSolutions( propId, currentSolutions[ct_ct[1]][ct_ct[0]], currentSolutions[ct_ct[3]][ct_ct[2]], solution01 );
                  }

                  for (int k=0; k<dsgX.size();k++)
                    dsgX[k]->storeReleaseProp( pid, -1);

#ifdef __Energy_Check__
                  {
                    std::vector<Scalar> dataEnergies(dsg_data.size(),0);
                    for (int k=0; k<dsg_data.size();k++)
                      dataEnergies[k] = dsgX[ dsg_data[k] ]->getEnergy();
                    Scalar dataEnergy = std::accumulate(dataEnergies.begin(),dataEnergies.end(),0.);

                    std::vector<Scalar> smoothEnergies(nCams*nTimeSteps,0);
#ifdef _SumViewSmoothing_
                    for (int k=0; k<evalSmooth.size();k++)
#else
                    for (int k=0; k<1;k++)
#endif
                      smoothEnergies[k] = getSmoothEnergy ( currentSolutions[k/nCams][k%nCams], nSegments[k/nCams][k%nCams], evalSmooth[k], theta, lambda);
                    Scalar smoothEnergy = std::accumulate(smoothEnergies.begin(),smoothEnergies.end(),0.);
                    Scalar ENERGY2 = dataEnergy + smoothEnergy;
                    if (ENERGY2 - __Small_Improvement__ > ENERGYOLD)
                    {
                      printf("Energies old:%f  (full=data+smooth) %f = (%f,%f,%f,%f),(%f,%f,%f,%f))\n", ENERGYOLD, ENERGY2, 
                        dataEnergies[0], dataEnergies[1], dataEnergies[2], dataEnergies[3], smoothEnergies[0], smoothEnergies[1], smoothEnergies[2], smoothEnergies[3]);
//                      mexEvalString("drawnow");
                    }
                    ENERGYOLD = ENERGY2;
                  }
#endif
                } // nVars >0
              }// while not stop
              // once only:
#ifdef _DEBUG
              OLD_ENERGY = ENERGY;// debug mode 1 iteration only -- otherwise takes ages
#endif
            }// parallel

            // energy update: - rather the last iteration NOT THE LAST ord TODO
            //      if (ord == order.size()-1)
            {
              Scalar dataEnergy = 0;

			  // update the penalties for occlusion/oob by current data term
              if ( doAuto > 0.)
                for (int k=0; k<dsg_data.size();k++)
                  dsgX[ dsg_data[k] ]->updateAutoScores();

              std::vector<Scalar> dataEnergies(dsg_data.size(),0);
              for (int k=0; k<dsg_data.size();k++)
                dataEnergy += dsgX[ dsg_data[k] ]->getEnergy();

              std::vector<Scalar> smoothEnergies(nCams*nTimeSteps,0);
#ifdef _SumViewSmoothing_
              for (int k=0; k<evalSmooth.size();k++)
#else
              for (int k=0; k<1;k++)
#endif
                smoothEnergies[k] = getSmoothEnergy ( currentSolutions[k/nCams][k%nCams], nSegments[k/nCams][k%nCams], evalSmooth[k], theta, lambda);
              Scalar smoothEnergy = std::accumulate(smoothEnergies.begin(),smoothEnergies.end(),0.);

#ifdef _writeEnergies_
#ifndef _testLocalReplacement_
              for(int kk= 0; kk< smoothEnergies.size(); kk++ )
                printf("cam/time (%d,%d) smooth Score: %.2f\n", kk/nCams, kk%nCams, smoothEnergies[kk]);
#endif
#endif

              OLD_ENERGY = ENERGY;
              ENERGY = dataEnergy + smoothEnergy;

#ifdef _writeEnergies_
#ifdef _SumViewSmoothing_
              printf("Energies (full=data+smooth) %.2f = (%.2f,(%.2f,%.2f,%.2f,%.2f))\n", ENERGY, dataEnergy, smoothEnergies[0], smoothEnergies[1], smoothEnergies[2], smoothEnergies[3]);
#else
              printf("Energies (full=data+smooth) %f = (%f,%f)\n", ENERGY, dataEnergy, smoothEnergies[0] );
#endif
#endif
            }

          }// outer while energy does get better loop

          std::clock_t stopP3(std::clock());
          printf("QPBO took %.2f\n", double(stopP3-startP3)/CLOCKS_PER_SEC );

#ifdef _testLocalReplacement_
          lrp.clear();lrp.generateMore();
          if (locRepIts>0)
            enStore.flipNormals();// stupid internal assumption - but to hard to change everywhere
          prop2PlaneMap       = lrp.getPropMap();
          expCenters          = lrp.getExpCenters();
          nProposalExpansions = lrp.getNProposals();
          nNormals            = (enStore.getNormals())->size();
//        while(locRepIts>3)
          }
#endif

#ifdef _writeEnergies_
          //  printf("ENERGY: %f  NON submodularity detected %d, no solutions: %d\n", ENERGY, non_sub, non_sol);
          printf("ENERGY: %.2f  NON submodularity detected %d (%.2f pct), no solutions: %d (%.2f pct); \n", ENERGY, non_sub, double(100*non_sub)/std::max(1,avEdgesUsed), non_sol, double(100*non_sol)/std::max(1,avNodesUsed));
          printf("av/max Nodes :%d/%d, av/max Edges: %d/%d \n", avNodesUsed/std::max(1,nRuns), maxNodesUsed, avEdgesUsed/std::max(1,nRuns), maxEdgesUsed);
#endif

          if (nlhs > 0)
          {
            int temp=0;
            for( int k=0;k<nTimeSteps;k++ )
              for(int i=0;i<nCams;i++ )
                temp += nSegments[k][i];

            plhs[0] =  mxCreateDoubleMatrix(temp, 1, mxREAL);
            Scalar* output  = (Scalar*) mxGetPr(plhs[0]);
            temp=0;
            for( int k=0;k<nTimeSteps;k++ )
              for(int j=0;j<nCams;j++ )
                for(int i=0;i<nSegments[k][j];i++,temp++ )
                  output[temp] = currentSolutions[k][j][i];
          }

          if (nlhs > 1)
          {
            // new:
            std::vector<Scalar> depthMap2(N*M);
            mwSize dims[3];
            dims[0] = M;
            dims[1] = N;
            dims[2] = nCams*nTimeSteps;//dsg_data.size();
            plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            Scalar* outputI = (Scalar*) mxGetPr(plhs[1]);
            memset( outputI, 0, sizeof(double)*dims[2]*dims[1]*dims[0] );

            for (int k=0;k<dsg_data.size();k++)
            {
              Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
              int id1 = ct_ct[0] + ct_ct[1] * nCams;
              int id2 = ct_ct[2] + ct_ct[3] * nCams;

              if (outputI[id1*N*M] == 0)
              {
                dsgX[ dsg_data[k] ]->buildDepthMap( depthMap2, true );
                for (int i=0;i < N*M; i++)
                  outputI[i + id1 * N*M] = depthMap2[i];
              }
              if (outputI[id2*N*M] == 0)
              {
                dsgX[ dsg_data[k] ]->buildDepthMap2( depthMap2, true );
                for (int i=0;i < N*M; i++)
                  outputI[i + id2*N*M] = depthMap2[i];
              }
            }
          }
          if (nlhs > 2)
          {
            // new
            std::vector<int> dataEnergyMap(N*M);
            std::vector<int> dataEnergyMap2(N*M);
            mwSize dims[3];
            dims[0] = M;
            dims[1] = N;
            dims[2] = dsg_data.size()*2;
            plhs[2] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
            int* outputI = (int*) mxGetPr(plhs[2]);

            for (int k=0;k<dsg_data.size();k++)
            {
              Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
              dsgX[ dsg_data[k] ]->getEnergyMap( dataEnergyMap, dataEnergyMap2 );

              for (int i=0;i < N*M; i++)
                outputI[i+    2*k*N*M] = dataEnergyMap [ (segImg[ct_ct[1]][ct_ct[0]])[i] ];
              for (int i=0;i < N*M; i++)
                outputI[i+N*M+2*k*N*M] = dataEnergyMap2[ (segImg[ct_ct[3]][ct_ct[2]])[i] ];
            }
          }

          if (nlhs > 3)
          {
            std::vector<Scalar> dataEnergyMap(N*M);
            std::vector<Scalar> dataEnergyMap2(N*M);
            mwSize dims[3];
            dims[0] = M;
            dims[1] = N;
            dims[2] = dsg_data.size()*2;
            plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            Scalar* outputI = (Scalar*) mxGetPr(plhs[3]);

            for (int k=0;k<dsg_data.size();k++)
            {
              Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
              dsgX[ dsg_data[k] ]->getDataEnergyMap( dataEnergyMap, dataEnergyMap2 );

              for (int i=0;i < N*M; i++)
                outputI[i+    2*k*N*M] = dataEnergyMap [ (segImg[ct_ct[1]][ct_ct[0]])[i] ];
              for (int i=0;i < N*M; i++)
                outputI[i+N*M+2*k*N*M] = dataEnergyMap2[ (segImg[ct_ct[3]][ct_ct[2]])[i] ];
            }
          }

          if (nlhs > 4)
          {
            mwSize dims[2];
            dims[0] = 1;
            dims[1] = 6;
            plhs[4] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
            Scalar* outputI = (Scalar*) mxGetPr(plhs[4]);
            outputI[0] = ENERGY; outputI[1] = non_sub; outputI[2] = non_sol; 
            outputI[3] = avNodesUsed; outputI[4] = avEdgesUsed; outputI[5] = nRuns;
          }

#ifdef _writeEnergies_
          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            Scalar endE = dsgX[ dsg_data[k] ]->getEnergy();
            printf("Data Energy:%.2f time:%d cam:%d <-> time:%d cam:%d\n", endE, ct_ct[1], ct_ct[0], ct_ct[3], ct_ct[2]);
          }
#endif

          enStore.flipNormals();

#ifdef _writeEnergies_
#pragma omp parallel for
          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            dsgX[ dsg_data[k] ]->buildFromCurrentSolution(  currentSolutions[ct_ct[1]][ct_ct[0]], currentSolutions[ct_ct[3]][ct_ct[2]] );
          }

          for (int k=0; k<dsg_data.size();k++)
          {
            Datasegments<Scalar>::P4i ct_ct = dsgX[dsg_data[k]]->getTimeCam();
            Scalar endE = dsgX[ dsg_data[k] ]->getEnergy();
            printf("Data Energy:%.2f time:%d cam:%d <-> time:%d cam:%d\n", endE, ct_ct[1], ct_ct[0], ct_ct[3], ct_ct[2]);
          }
#endif
          for( int i=0;i<dsgX.size();i++)
            delete dsgX[i];

#ifdef _testLocalReplacement_
  if (nlhs > 5)
  {
    std::vector< EvalEnergyStore<Scalar>::P3>* norStored = enStore.manipulateNormals();
    plhs[5]           =  mxCreateDoubleMatrix(3, (*norStored).size(), mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[5]);
    for (int i=0;i < (*norStored).size(); i++)
    {
      EvalEnergyStore<Scalar>::P3 No = (*norStored)[i];
      outputDI[3*i]   = No[0];
      outputDI[3*i+1] = No[1];
      outputDI[3*i+2] = No[2];
    }
  }
  if (nlhs > 6)
  {
    mwSize dims[3] = {4,4,1};
    std::vector< EvalEnergyStore<Scalar>::P3>* traStored = enStore.manipulateTra( );
    std::vector< EvalEnergyStore<Scalar>::M3>* rotStored = enStore.manipulateRot( );
    dims[2] = (*traStored).size();
    plhs[6]           =  mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[6]);
    for (int i=0;i < (*traStored).size(); i++)
    {
      EvalEnergyStore<Scalar>::M3 R = (*rotStored)[i];
      EvalEnergyStore<Scalar>::P3 T = (*traStored)[i];
      outputDI[16*i  ]  = R(0,0);
      outputDI[16*i+1 ] = R(1,0);
      outputDI[16*i+2 ] = R(2,0);
      outputDI[16*i+3 ] = Scalar(0);
      outputDI[16*i+4 ] = R(0,1);
      outputDI[16*i+5 ] = R(1,1);
      outputDI[16*i+6 ] = R(2,1);
      outputDI[16*i+7 ] = Scalar(0);
      outputDI[16*i+8 ] = R(0,2);
      outputDI[16*i+9 ] = R(1,2);
      outputDI[16*i+10] = R(2,2);
      outputDI[16*i+11] = Scalar(0);
      outputDI[16*i+12] = T[0];
      outputDI[16*i+13] = T[1];
      outputDI[16*i+14] = T[2];
      outputDI[16*i+15] = Scalar(1);
    }
  }

  if (nlhs > 7)
  {
    // could also return centers here 
    std::vector<EvalEnergyStore<Scalar>::P3>& mvpCenters =  lrp.getMvpCenters();
    plhs[7]           =  mxCreateDoubleMatrix(3, mvpCenters.size(), mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[7]);
    for (int i=0;i < nProposalExpansions; i++)
    {
      outputDI[3*i  ] = mvpCenters[i][0];
      outputDI[3*i+1] = mvpCenters[i][1];
      outputDI[3*i+2] = mvpCenters[i][2];
    }
  }

#endif

}
#endif // __VCSF_SEG_GC__
