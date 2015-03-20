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

#ifndef __DataGeneral__cpp
#define __DataGeneral__cpp

#include "DataSegments.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;
/////////////////////////////////////////////////////////

/// first project all centers and sort them - actually no projection needed
template<class Scalar>
void
  Datasegments<Scalar>::
  buildList ( M3& KMat, Scalar* centers, int nSegments, typename std::vector< PatchCenter<Scalar> >& projPatchCenters)
{
  projPatchCenters.clear();
  projPatchCenters.resize( nSegments );

  // first project current solutions:
  for (int i=0;i<nSegments;i++)
  {
    P3 p0     = KMat * P3( &centers[3*i] );// center is in camera coords
    p0 /= p0[2];
    //      projPatchCenters[i]   = PatchCenter<Scalar>( int(p0[0]-0.5), int (p0[1]-0.5), P3( &centers[3*i] ), i );
    projPatchCenters[i]   = PatchCenter<Scalar>( int(floor(p0[0]-0.5)), int (floor(p0[1]-0.5)), p0, i );
  }
  std::sort( projPatchCenters.begin(), projPatchCenters.end() );
  // need also a map where a certain segment is now, to start searching
}

template<class Scalar>
void
  Datasegments<Scalar>::
  setKs ( Scalar* Kl_, Scalar* Kl2_  )
{
  Kl    =  M3(Kl_);
  Kl2   =  M3(Kl2_);
  iKl   = Kl;iKl.invert();
  iKlt  = iKl;iKlt=iKlt.transpose();
  iKl2  = Kl2;iKl2.invert();
  iKlt2 = iKl2;iKlt2=iKlt2.transpose();
  //    Klt   = Kl.transpose();
  //    Klt2  = Kl2.transpose();
  iiKlt2= iKlt2;iiKlt2.invert();
}

template<class Scalar>
void
  Datasegments<Scalar>::
  setAreas()
{
  assert ( segImg2 != NULL );
  assert ( segImg  != NULL );

  areas1.assign(nSegments, 0);
  areas2.assign(nSegments2,0);

  for (int i=0; i < w*h; i++)
  {
    areas1[ segImg  [i] ] ++;
    areas2[ segImg2 [i] ] ++;
  }
}

template<class Scalar>
int Datasegments<Scalar>::
whichPenalty ( M3 Hom, P3 Nom, P3 iNom, P3 centre, M3& iHom, P3 centre2, int prop1, int prop2 )
{
  P3 p = Hom * centre;p /= p[2];//projected point
  Scalar depth0      = 1./(p|Nom); 
  Scalar depthThresh = max ( maxDepthThresh, fabs(depth0  * __OCC_THRESH__) );
  Scalar depthDiff   = depth0 - 1./(p|iNom);// let depths by positive -> depthDiff>0 implies an occlusion case, sign change: < 0 occlusion case

  if ( depthDiff > depthThresh ) // projection occludes other point but the other is visible: impossible
    return 2; // imp pen cyan
  else
    if ( depthDiff < -depthThresh ) // projection is occluded by other point: occlusion model
      return 1; // occ, green hsv 5 bins
    else
    {
      if (prop2 != prop1)
      {
        return 3;// match but different label, dunkelblau oder lila
      }
      return 0; // data, yellow same label (auch ausserhalb)
    }
}

/// debug function
template<class Scalar>
void
  Datasegments<Scalar>::
  checkSolution(const std::vector<int>& _currentSolution, const std::vector<int> &_currentSolution2 )
{
  for (int i=0;i < currentSolution1.size(); i++)
  {
    if (currentSolution1[i] != _currentSolution[i])
      printf("seg i: %d, %d", i, currentSolution1[i], _currentSolution[i]);
  }
  for (int i=0;i < currentSolution2.size(); i++)
  {
    if (currentSolution2[i] != _currentSolution2[i])
      printf("seg j: %d, %d", i, currentSolution2[i], _currentSolution2[i]);
  }
}

template<class Scalar>
void
  Datasegments<Scalar>::
  setCentersKs( Scalar* Kl_, Scalar* Kl2_, Scalar* _centers1, Scalar* _centers2, int nCenters )
{ 
  setKs ( Kl_, Kl2_  ); // those need to be there
  if (_centers1 ==NULL || _centers2 ==NULL ) return;
  centers = _centers1; centers2 = _centers2;
  Centers.clear();// patch centers in view1
  iCenters.clear();// patch centers in view2

  if (nCenters ==0)
    nCenters = nSegments;

  Centers.resize(nCenters);
  iCenters.resize(nCenters);

  for(int i =0;i<nCenters; i++ )
  {
    P3 p0 = Kl * P3( &centers[3*i] );
    p0 /= p0[2];
    Centers[i] = p0;
  }
  for(int i =0;i<nCenters; i++ )
  {
    P3 p0 = Kl2 * P3( &centers2[3*i] );
    p0 /= p0[2];
    iCenters[i] = p0;
  }
}

/// create debug output for visualization
template<class Scalar>
void
  Datasegments<Scalar>::
  buildDepthMap2( std::vector<Scalar>& depthMap, bool remap )
{
  depthMap.assign(w*h,0);
  for (int i =0; i< w*h;i++)
  {
    // this is ok for per pixel, here the data elements are identical to the segmentIds
    int propId = segImg2[ i ];// i need the mapping for this you idiot - the solution AAH
    if (remap)
    {
      propId = currentSolution2[propId];//this is actually a normal/rot 
      //         propId = prop2Data[propId];       // now its one of dataElements and so the iNom/Hom/ .. stored there
      depthMap[i] = 1./ (inoms_pix[propId] | P3FromGlobal( i ));
    }
    else
      // from segment to normal, but dataElements here is stored per segment (if per segi optimization)
      // per pixel optimization, then fine.
      depthMap[i] = 1./ (dataElements[propId].iNom | P3FromGlobal( i ));
  }
}

/// create debug output for visualization
template<class Scalar>
void
  Datasegments<Scalar>::
  buildDepthMap( std::vector<Scalar>& depthMap, bool remap )
{
  depthMap.assign(w*h,0);
  for (int i =0; i< w*h;i++)
  {
    // this is ok for per pixel, here the data elements are identical to the segmentIds
    int propId = segImg[ i ];// i need the mapping for this you idiot - the solution AAH
    if (remap)
    {
      propId = currentSolution1[propId];//this is actually a normal/rot 
      depthMap[i] = 1./ (noms_pix[propId] | P3FromGlobal( i ));
    }
    else
      // from segment to normal, but dataElements here is stored per segment (if per segi optimization)
      // per pixel optimization, then fine.
      depthMap[i] = 1./ (dataElements[propId].Nom | P3FromGlobal( i ));
  }
}

/// visualization
template<class Scalar>
void
  Datasegments<Scalar>::getAutoScores( std::vector<Scalar> &energyMap, std::vector<Scalar> &energyMap2 )
{
  energyMap.resize(nSegments,-1);
  energyMap2.resize(nSegments2,-1);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < nSegments; i++ )
  {
    int segId = i;

    int prop0   = currentSolution1[ segId ];
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// projects onto

    energyMap[i] = autoScoresFW ? areas1[i]*adaptiveAutoScoresFW[segId] : 1.;
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < nSegments2; i++ )
  {
    int segId = i;

    int prop0   = currentSolution2[ segId ];
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment

    energyMap2[i] = autoScoresBW ? areas2[i]*adaptiveAutoScoresBW[segId] : 1.;
  }/////////////////
#endif
}

#endif