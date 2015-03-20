/*!
Here we project moving planes to a different view (in time or a different camera or both).
We distinguish between simple projection of the moving planes (all of them),
or projection into the other image's segmentation.
In the latter case each segment in the other view then holds a (changed) moving plane
from the original view.
Holes (Segments without a moving plane) are filled by employing a priority Q.
We project moving planes not their labels.
The functionality is thus different from ProjectSegmentation.h.
*/
/// here we allow proposals to be moved to theri neighbours
/// the selection is taken randomly, we pick a segment, distribute the motion
/// to the neeighbours and fix these
/// finally we stop if no segment can be picked any more
/// we do this N times
/// also possible Growing seeds like in patchMatch
#ifndef __LINEAR_MOTION_PROPOSALS__
#define __LINEAR_MOTION_PROPOSALS__
////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include "EvalEnergyStore.h"
#include "linearProposals.h"

#include <algorithm>
#include <list>

#include "mex.h" // include after VectorT and Mat3x3T

#include <ctime>

using namespace std;
using namespace Math;

typedef double Scalar;

/// can either return the projections only or the projections fitted into the new segmentation
//#define _returnProjectionOnly_

// propagate the 3d moving planes to the next timestep
// 2 seg images are used, at t0 and t1
// for all 3d moving planes, project center to other image == proposal
// compute for all centers of t1 segmented image
// closest proposal and store in Q consider only neighs of segment the proposal is projected into
// use priority Q to solve the rest, then compute N, by M3 inversion and R|T remain constant ?!
void generateLinearMotionProposals( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  typedef unsigned int iType;
  typedef double Scalar;

  size_t N,M;

  if (nrhs < 7)
    mexErrMsgTxt("Must have at least 7 input arguments");

  int* segImg_t1(NULL);
  int* segImg(NULL);

  const mxArray* edges    = (prhs[0]);
  Scalar* K_      = (Scalar*)  mxGetPr(prhs[1]); // K_r

  Scalar* normals    = (Scalar*) mxGetPr(prhs[ 2]);

  Scalar* rotations  =  (Scalar*) mxGetPr (prhs[3]);
  Scalar* centers    =  (Scalar*) mxGetPr (prhs[4]);
  Scalar* centers_t1 =  (Scalar*) mxGetPr (prhs[5]);

  segImg_t1   =  (int*) mxGetPr(prhs[6]);

  int withProjection=0;
  if (nrhs > 7)
    withProjection = (int) (*mxGetPr(prhs[ 7]));
  //////////////////////////

  M           = (size_t) mxGetM(prhs[6]);
  N           = (size_t) mxGetN(prhs[6]);

  int nNormals     = (int) max( mxGetM(prhs[ 2]), mxGetN(prhs[ 2]) );
  int nStep        = (int) min( mxGetM(prhs[ 2]), mxGetN(prhs[ 2]) );
  int nSegments_t1 = (int) max( mxGetM(edges), mxGetN(edges) );
  int nSegments_t0 = (int) max( mxGetM(prhs[ 4]), mxGetN(prhs[ 4]) );

  printf ("nSegments t0-%d t1-%d, nNormals %d \n", nSegments_t0, nSegments_t1, nNormals);

  ///////////////////////////////////////////////////////
  std::vector<int> currentSolution( nSegments_t0, 0 );

  if ( nNormals >= nSegments_t0 )
    for ( int i = 0;i < nSegments_t0;i++ )
      currentSolution[i]=i;

  ///////////////////////////////////////////////////////

  EvalEnergyStore<Scalar> enStore;
  enStore.prepareNormals( normals, nNormals, nStep );
  enStore.setRotTra( rotations, nNormals );

  enStore.flipNormals();// ?

  std::clock_t startP(std::clock());

  linearMotionProposals<Scalar> LMH( nSegments_t0, nNormals, (int) N, (int) M );

//  LMH.useCentred(false);
#ifdef _useNoCenters_
  LMH.useCentred(false);
#endif

  LMH.setK(K_);
  LMH.setNormals( enStore.getNormals() );
  LMH.setRotTra( enStore.getRot(), enStore.getTra() );
  if (withProjection)
    LMH.noMotionMovingPlanes( nSegments_t0, centers, currentSolution );
  else
    LMH.computeNew3dmovingplanes( nSegments_t0, centers, currentSolution );



#ifdef  _returnProjectionOnly_
  LMH.adjustNormals();
  std::vector<linearMotionProposals<Scalar>::P3> newNormals      = LMH.getProjN();
  std::vector<linearMotionProposals<Scalar>::P3> newTranslations = LMH.getProjTra();
  std::vector<linearMotionProposals<Scalar>::M3> newRotations    = LMH.getProjRot();
  std::vector<linearMotionProposals<Scalar>::P2> newCenters      = LMH.getNewCenters();
#else
  LMH.generateProposals(edges, segImg_t1, nSegments_t1, centers_t1);

  std::vector<linearMotionProposals<Scalar>::P3> newNormals      = LMH.getNewN();
  std::vector<linearMotionProposals<Scalar>::P3> newTranslations = LMH.getNewTra();
  std::vector<linearMotionProposals<Scalar>::M3> newRotations    = LMH.getNewRot();
#endif
  std::clock_t stopP(std::clock());
  printf("Preprocess took %f\n", double(stopP-startP)/CLOCKS_PER_SEC );

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  if (nlhs > 0)
  {
    plhs[0]           =  mxCreateDoubleMatrix(3, (int) newNormals.size(), mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[0]);
    for (int i=0;i < newNormals.size(); i++)
    {
      linearMotionProposals<Scalar>::P3 No = newNormals[i];
      outputDI[3*i]   = No[0];
      outputDI[3*i+1] = No[1];
      outputDI[3*i+2] = No[2];
    }
  }
  if (nlhs > 1)
  {
    mwSize dims[3] = {4,4,1};
    dims[2] = (mwSize) newTranslations.size();
    plhs[1]           =  mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[1]);
    for (int i=0;i < newTranslations.size(); i++)
    {
      linearMotionProposals<Scalar>::M3 R = newRotations[i];
      linearMotionProposals<Scalar>::P3 T = newTranslations[i];
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

#ifdef  _returnProjectionOnly_
  if (nlhs > 2)
  {
    plhs[2]           =  mxCreateDoubleMatrix(2,newNormals.size(), mxREAL);
    Scalar *outputDI  = (Scalar*) mxGetPr(plhs[2]);
    for (int i=0;i < newCenters.size(); i++)
    {
      linearMotionProposals<Scalar>::P2 No = newCenters[i];
      outputDI[2*i]   = No[0];
      outputDI[2*i+1] = No[1];
    }
  }
#endif

}
#endif // __LINEAR_MOTION_PROPOSALS__
