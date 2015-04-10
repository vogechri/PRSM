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

#ifndef __DataPerSegment__cpp
#define __DataPerSegment__cpp

#include "DataSegments.h"
#include "DataDefinitionsVC.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;
///////////////////////////////////////

template<class Scalar>
bool Datasegments<Scalar>::
  testNewProp( int proposal )
{
  P4i boxFW = dataElements[proposal].boxFW;
  P4i boxBW = dataElements[proposal].boxBW;

  int puffer = 0;

  // 1. compare with elements under revision:
  for (int i=0;i<storedProp.size();i++)
    if (storedProp[i] != -1)
    {
      int prid = storedProp[i];
      if ( (boxFW[0] < dataElements[prid].boxFW[2]+puffer && puffer+boxFW[2] > dataElements[prid].boxFW[0]) && 
        (boxFW[1] < dataElements[prid].boxFW[3]+puffer && puffer+boxFW[3] > dataElements[prid].boxFW[1]) )
        return false;
      if ( (boxBW[0] < dataElements[prid].boxBW[2]+puffer && puffer+boxBW[2] > dataElements[prid].boxBW[0]) && 
        (boxBW[1] < dataElements[prid].boxBW[3]+puffer && puffer+boxBW[3] > dataElements[prid].boxBW[1]) )
        return false;
    }
    return true;
}

template<class Scalar>
void Datasegments<Scalar>::
  set_AutoScores(Scalar* _autoScores1, Scalar* _autoScores2)
{
  autoScoresFW = _autoScores1;autoScoresBW = _autoScores2;
  gs1.setAutoScores( autoScoresFW );
  gs2.setAutoScores( autoScoresBW );
  updateAutoScores();
}

#ifndef _checkValiadityAuto_
template<class Scalar>
void Datasegments<Scalar>::
  updateAutoScores( )
{
  const int minSize = 10;

  if (autoScoresFW != NULL)
  {
    if ( adaptiveAutoScoresFW.size() == 0 )
      adaptiveAutoScoresFW.insert(adaptiveAutoScoresFW.end(), &autoScoresFW[0], &autoScoresFW[nSegments]);
    else
    {
      for (int i=0;i< adaptiveAutoScoresFW.size(); i++)
      {
        if (currentSolutionE.freeFW[i] > minSize) // and correspondence valid, ie. no occlusion
          adaptiveAutoScoresFW[i] = std::min( adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i]/(currentSolutionE.freeFW[i]) );
      }
      // dataElements are given per proposal ! not per segment
      for (int i=0;i< dataElements.size(); i++) 
        dataElements[i].updateOobPenaltyFW ( adaptiveAutoScoresFW, areas1 );
    }
    //        std::transform (adaptiveAutoScoresFW.begin(), adaptiveAutoScoresFW.end(), currentSolutionE.dataFW.begin(), adaptiveAutoScoresFW.begin(), std::min<Scalar>);
  }

  if (autoScoresBW != NULL)
  {
    if ( adaptiveAutoScoresBW.size() == 0 )
      adaptiveAutoScoresBW.insert(adaptiveAutoScoresBW.end(), &autoScoresBW[0], &autoScoresBW[nSegments2]);
    else
    {
      for (int i=0;i< adaptiveAutoScoresBW.size(); i++)
      {
        if (currentSolutionE.freeBW[i] > minSize) // and correspondence valid
          adaptiveAutoScoresBW[i] = std::min( adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i]/(currentSolutionE.freeBW[i]) );
      }
      for (int i=0;i< dataElements.size(); i++)
        dataElements[i].updateOobPenaltyBW ( adaptiveAutoScoresBW, areas2 );
    }
  }
}
#else
template<class Scalar>
void Datasegments<Scalar>::
  updateAutoScores( )
{
  const int minSize = 10;

  if (autoScoresFW != NULL)
  {
    if ( adaptiveAutoScoresFW.size() == 0 )
      adaptiveAutoScoresFW.insert(adaptiveAutoScoresFW.end(), &autoScoresFW[0], &autoScoresFW[nSegments]);
    else
    {
      for (int i=0;i< adaptiveAutoScoresFW.size(); i++)
      {
        if (currentSolutionE.freeFW[i] > minSize) // and correspondence valid, ie. no occlusion
        {
          int prop0   = currentSolution1[ i ];
          int seg2_0  = currentSolutionE.seg2segFW[ i ].second;// projects onto
          if (seg2_0 >=0) // maps on other segment:
          {
            int prop00 = currentSolution2[seg2_0];// here only if two valid datas
            if( prop00 == prop0 || whichPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ i ], ihoms_pix[prop00], iCenters[ seg2_0 ], prop00, prop0 ) == 3) // 0 see above
              adaptiveAutoScoresFW[i] = std::min( adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i]/ Scalar(currentSolutionE.freeFW[i]) );
//            else
//              printf("Check: %.2f, %.2f \n", adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i]/Scalar(currentSolutionE.freeFW[i]) );
          }
          else // still data 
          {
//            printf("Check: %.2f, %.2f \n", adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i]/Scalar(currentSolutionE.freeFW[i]) );
            adaptiveAutoScoresFW[i] = std::min( adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i]/ Scalar(currentSolutionE.freeFW[i]) );
          }
        }
      }
      // dataElements are given per proposal ! not per segment
      for (int i=0;i< dataElements.size(); i++) 
      {
        dataElements[i].updateOobPenaltyFW ( adaptiveAutoScoresFW, areas1 );
#ifdef __depthControl__
        for(int j=0; j<(dataElements[i].seg2segFW).size(); j++ )
        {
          if (dataElements[i].freeFW[j] <=0)
          {
            int sid = (dataElements[i].seg2segFW)[j].first;
            Scalar depth0      = 1./(Centers[ sid ]|noms_pix[dataElements[i].proposal]); 
            if ( depth0 > Scalar(-_minDepth_) )
            {
              Scalar impPen = areas1[sid] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresFW ? adaptiveAutoScoresFW[sid] : 1.)); // special case segment behind cam -- thus invisible
              dataElements[i].dataFW[j] = impPen;
            }
          }
        }
#endif
     }
    }
    //        std::transform (adaptiveAutoScoresFW.begin(), adaptiveAutoScoresFW.end(), currentSolutionE.dataFW.begin(), adaptiveAutoScoresFW.begin(), std::min<Scalar>);
  }

  if (autoScoresBW != NULL)
  {
    if ( adaptiveAutoScoresBW.size() == 0 )
      adaptiveAutoScoresBW.insert(adaptiveAutoScoresBW.end(), &autoScoresBW[0], &autoScoresBW[nSegments2]);
    else
    {
      for (int i=0;i< adaptiveAutoScoresBW.size(); i++)
      {
        if (currentSolutionE.freeBW[i] > minSize) // and correspondence valid
        {
          int prop0   = currentSolution2[ i ];
          int seg2_0  = currentSolutionE.seg2segBW[ i ].second;// projects onto

          if (seg2_0 >=0) // maps on other segment:
          {
            int prop00 = currentSolution1[seg2_0];
                                   
            if( prop00 == prop0 || whichPenalty ( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ i ], homs_pix[prop00],  Centers[ seg2_0 ], prop00, prop0 ) == 3) // 0 see above
              adaptiveAutoScoresBW[i] = std::min( adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i]/Scalar(currentSolutionE.freeBW[i]) );
//            else
//              printf("Check: %.2f, %.2f \n", adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i]/Scalar(currentSolutionE.freeBW[i]));
          }
          else // still data?
          {
//            printf("Check: %.2f, %.2f \n", adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i]/Scalar(currentSolutionE.freeBW[i]));
            adaptiveAutoScoresBW[i] = std::min( adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i]/Scalar(currentSolutionE.freeBW[i]) );
          }
        }
      }
      for (int i=0;i< dataElements.size(); i++) 
      {
        dataElements[i].updateOobPenaltyBW ( adaptiveAutoScoresBW, areas2 );
#ifdef __depthControl__
        for(int j=0; j<(dataElements[i].seg2segBW).size(); j++ )
        {
          if (dataElements[i].freeBW[j] <=0)
          {
            int sid = (dataElements[i].seg2segBW)[j].first;
            Scalar depth0      = 1./(iCenters[ sid ]|inoms_pix[dataElements[i].proposal]); 
            if ( depth0 > Scalar(-_minDepth_) )
            {
              Scalar impPen = areas2[sid] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresFW ? adaptiveAutoScoresFW[sid] : 1.)); // special case segment behind cam -- thus invisible
              dataElements[i].dataBW[j] = impPen;
            }
          }
        }
#endif
      }
    }
  }
}
#endif
/// init homographies and normals per view
template<class Scalar>
void
  Datasegments<Scalar>::
  initHoms( int normalId, genHomoG<Scalar>* gHom, const std::vector<P3>* normals )
{
  ////// first project the bbox to get a box in the second image:
  M3 Hom    = gHom->getHomC0( normalId );
  M3 iHom   = Hom;iHom.invert();
  M3 HomCC  = Hom * Kl;
  M3 iHomCC = iHom * Kl2;

  P3 vn_0   = gHom->getViewNormalC0( normalId );

  if (vNoms_cam.size() <= normalId )
  {
    vNoms_cam.resize( normalId+1, P3(0,0,0) );
    homs_pix.resize ( normalId+1 , M3(0,0,0, 0,0,0, 0,0,0) );
    ihoms_pix.resize( normalId+1 , M3(0,0,0, 0,0,0, 0,0,0) );
    noms_pix.resize ( normalId+1 , P3(0,0,0) );
    inoms_pix.resize( normalId+1 , P3(0,0,0) );
  }
  vNoms_cam[normalId] = -vn_0;
  homs_pix [normalId] = Hom;
  ihoms_pix[normalId] = iHom;
  noms_pix [normalId] = iKlt  * (*normals)[normalId];
  inoms_pix[normalId] = iKlt2 * vn_0;
}

/// init homographies and normals  of pairs not involving the canonical view - using the output form these pairs
template<class Scalar>
void
  Datasegments<Scalar>::
  initHoms( int normalId, const std::vector<M3>& ihoms1, const std::vector<M3>& homs2, const std::vector<P3>& inoms1, const std::vector<P3>& inoms2 )
{
  ////// first project the bbox to get a box in the second image:
  M3 Hom    = homs2[ normalId ] * ihoms1[ normalId ];
  M3 iHom   = Hom;iHom.invert(); 

  if (vNoms_cam.size() <= normalId )
  {
    vNoms_cam.resize( normalId+1, P3(0,0,0) );
    homs_pix.resize ( normalId+1 , M3(0,0,0, 0,0,0, 0,0,0) );
    ihoms_pix.resize( normalId+1 , M3(0,0,0, 0,0,0, 0,0,0) );
    noms_pix.resize ( normalId+1 , P3(0,0,0) );
    inoms_pix.resize( normalId+1 , P3(0,0,0) );
  }
  vNoms_cam[normalId] =-(iiKlt2 * inoms2[normalId]);
  homs_pix [normalId] = Hom;
  ihoms_pix[normalId] = iHom;
  noms_pix [normalId] = inoms1[normalId];
  inoms_pix[normalId] = inoms2[normalId];
}

template<class Scalar>
void
  Datasegments<Scalar>::
  reorganizeProposals (std::vector<int>& prop2Plane, int& nProposalExpansions, 
  Scalar* expCenters, std::vector<Scalar>& propCenters, 
  const int boxX, const int boxY, 
  std::vector<typename dataElem<Scalar>::P2i>& patchXY, bool doNormalCheck )
{
  // step 1 pick a new proposal, find double occurencies, compute new bbox and center, remove double occurencies
  std::vector<bool> doubleOccurence (nProposalExpansions, false);
  propCenters.clear();
  std::vector<int> new_prop2Plane;
  patchXY.clear();
  int new_nProposalExpansions(0);

  // handles:
  // note 16x16 boxes, so dist 16 + eps == 4 boxes; etc.
  // 2*gridSize*gridSize + 1: // gridSize^2*2 and a bit covers 9 boxes - could also go for 4 by suming up centers included
  //if ((p0-p1).sqrnorm() > gridSize*gridSize + 10) // ideally covers 4 boxes - see above with 9
  const int maxCovering = gridSize*gridSize + 10;// qpbo this laptop: best choice
  //  const int maxCovering = 2*gridSize*gridSize + 10;

  for ( int i=0; i< nProposalExpansions; i++)
  {
    if (doubleOccurence[i]) continue;

    std::vector<int> twins;
    int propId = prop2Plane[i];
    twins.push_back ( i );

    P3 p0     = Kl * P3( &expCenters[3 * i] );p0/=p0[2];// center is in camera coords

    for ( int j=i+1; j< nProposalExpansions; j++)
    {
      if ( propId != prop2Plane[j] || doubleOccurence[j] )
        continue;

      if (!doNormalCheck)
      {
        P3 p1     = Kl * P3( &expCenters[3 * j] );p1/=p1[2];// center is in camera coords

        if ((p0-p1).sqrnorm() > maxCovering) // ideally covers 4 boxes - could also go for 4 by suming up centers included
          continue;
        p0 = (p0 * Scalar(twins.size()) / Scalar(twins.size()+1)) + (p1 / Scalar(twins.size()+1));
      }
      twins.push_back( j );
      doubleOccurence[j] = true;
    }
    // now join all twins into one model

    P4i box1(0,0,0,0);
    for (int k=0;k< twins.size(); k++)
    {
      int expansionArea =0;

      int twinPropID = prop2Plane[ twins[k] ];

      P3 p0;
      p0     = Kl * P3( &expCenters[3 * twins[k]] );// center is in camera coords
      p0 /= p0[2];

      int startX = int ( floor(p0[0] -0.5) ) - boxX;
      int endX   = int ( floor(p0[0] -0.5) ) + boxX;
      int startY = int ( floor(p0[1] -0.5) ) - boxY;
      int endY   = int ( floor(p0[1] -0.5) ) + boxY;
      // treatment else as usual

      // printf("k:%d twins[k]: twinPropID :%d %d newCenter: %0.2f %0.2f %0.2f,  box: %d,%d,%d,%d \n", k, twins[k], twinPropID, p0[0], p0[1], p0[2], startX, startY, endX, endY );

      //startX =    int(startX/gridSize)*gridSize;
      //endX   = (1+int(endX/gridSize))*gridSize;
      //startY =    int(startY/gridSize)*gridSize;
      //endY   = (1+int(endY/gridSize))*gridSize;

      if (k==0)
        box1 = P4i (startX, startY, endX, endY);
      else
      {
        box1.minimize( P4i (startX, startY, box1[2], box1[3]) );
        box1.maximize( P4i (box1[0], box1[1], endX, endY ) );
      }
    }

    if (doNormalCheck && noms_pix.size() > propId)
    {
      P3 normal = noms_pix [propId];
      P3 normalDir = normal;
      if ( fabs (normalDir[0]) + fabs (normalDir[1]) > 0.001 )
      {
        normalDir[2] = 0;normalDir.normalize();//!! not both == 0 else FAIL

        Scalar maxDepth = -1;
        P3 p0     = Kl * P3( &expCenters[3 * i] );// center is in camera coords
        p0 /= p0[2];

        Scalar ownDepth = 1./(normal|p0);

        Scalar scale = (Scalar (1) - maxDepth * (normal|p0) ) / (maxDepth * (normalDir|normal) );
        P3 displace = normalDir * scale;
        P3 maxPoint = p0 + displace;
        P3 minPoint = p0 - displace;

        Scalar dmax = 1./(maxPoint|normal);
        Scalar dmin = 1./(minPoint|normal);

        Scalar planeMax = maxPoint|normalDir;

        //          if ( displace.norm() < max(boxX, boxY) )
        {
          int dummyStep = 1;
          if (displace.norm() < 100)
            dummyStep =2;

          if (dmax < dmin)
            dummyStep =3;

          M2 mat( normalDir[0], normalDir[1], 0., 1. );
          if ( mat.invert() )
          {
            P2 fail1 = mat * P2(planeMax, box1[1]);
            P2 fail2 = mat * P2(planeMax, box1[3]);
            if ( normalDir[0] < 0 )// more right the lower
            {

              if (box1[2] > int (max ( fail1[0], fail2[0] )))
              {
#ifdef _writeDebugOut_
                printf("----- CUTTING ---- \n");
                printf("box 2 min: %d %d\n", box1[2], int (max ( fail1[0], fail2[0] )));

                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[3], 1)|normal));
#endif
                box1[2] = min( box1[2], int (max ( fail1[0], fail2[0] )) );
#ifdef _writeDebugOut_
                printf("Depth: %.2f\n", 1./(P3(fail1[0], fail1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[3], 1)|normal));
#endif
              }
            }
            else
            {
              if (box1[0] < int (min ( fail1[0], fail2[0] )))
              {
#ifdef _writeDebugOut_
                printf("--- CUTTING --- \n");
                printf("box 0 max: %d %d\n", box1[0], int (min ( fail1[0], fail2[0] )));

                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[3], 1)|normal));
#endif
                box1[0] = max( box1[0], int (min ( fail1[0], fail2[0] )) );
#ifdef _writeDebugOut_
                printf("Depth: %.2f\n", 1./(P3(fail1[0], fail1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[3], 1)|normal));
#endif
              }
            }
          }
          mat = M2( normalDir[0], normalDir[1], 1., 0. );
          if ( mat.invert() )
          {
            P2 fail1 = mat * P2(planeMax, box1[0]);
            P2 fail2 = mat * P2(planeMax, box1[2]);
            if ( normalDir[1] < 0 )
            {
              if (box1[3] > int (max ( fail1[1], fail2[1] )))
              {
#ifdef _writeDebugOut_
                printf("--- CUTTING --- \n");
                printf("box 3 min: %d %d\n", box1[3], int (max ( fail1[1], fail2[1] )));

                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[3], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[3], 1)|normal));
#endif
                box1[3] = min( box1[3], int (max ( fail1[1], fail2[1] ) ));
#ifdef _writeDebugOut_
                printf("Depth: %.2f\n", 1./(P3(fail1[0], fail1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[3], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[3], 1)|normal));
#endif
              }
            }
            else
            {

              if (box1[1] < int (min ( fail1[1], fail2[1] )))
              {
#ifdef _writeDebugOut_
                printf("--- CUTTING --- \n");
                printf("box 1 max: %d %d\n", box1[1], int (min ( fail1[1], fail2[1] )));

                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[1], 1)|normal));
#endif
                box1[1] = max( box1[1], int (min ( fail1[1], fail2[1] ) ));
#ifdef _writeDebugOut_
                printf("Depth: %.2f\n", 1./(P3(fail1[0], fail1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[0], box1[1], 1)|normal));
                printf("Depth: %.2f\n", 1./(P3(box1[2], box1[1], 1)|normal));
#endif
              }
            }
          }
        }
      }
    }

    box1.maximize( P4i(box1[0], box1[1], box1[0], box1[1] ));
    box1.minimize( P4i(box1[2], box1[3], box1[2], box1[3] ));

    P3 newCenter = iKl * P3( Scalar(box1[2]+box1[0]) /2., Scalar(box1[3]+box1[1]) /2., 1.);newCenter /= newCenter[2];
    // safety like this unconstraint nothing can go worse than before :-)
    //    int box_in_X = max( boxX, int(floor (Scalar(box1[2]-box1[0]) /2. + 0.5 )));
    //    int box_in_Y = max( boxY, int(floor (Scalar(box1[3]-box1[1]) /2. + 0.5 )));

    // here with cutting 
    // or .. ? no oob parts are not treated correctly but are cut for no reason
    int box_in_X = max( int(Scalar(boxX)/4000.), int(floor (Scalar(box1[2]-box1[0]) /2. + 0.5 )));
    int box_in_Y = max( int(Scalar(boxY)/4000.), int(floor (Scalar(box1[3]-box1[1]) /2. + 0.5 )));

    new_prop2Plane.push_back( propId );
    patchXY.push_back (typename dataElem<Scalar>::P2i (box_in_X, box_in_Y) );
    propCenters.push_back( newCenter[0] );
    propCenters.push_back( newCenter[1] );
    propCenters.push_back( newCenter[2] );

    //    printf("newCenter: %0.2f %0.2f %0.2f, (box_in_X, box_in_Y): (%d %d) box1: %d,%d,%d,%d \n", newCenter[0], newCenter[1], newCenter[2], box_in_X, box_in_Y, box1[0],box1[1], box1[2], box1[3] );
  }

  prop2Plane = new_prop2Plane;
  nProposalExpansions = prop2Plane.size();

#ifdef _writeDebugOut_
  printf("Reduced nProposals: %d,  \n", nProposalExpansions);
#endif
}


template<class Scalar>
void
  Datasegments<Scalar>::
  initNewSegFromView( int propId, int proposal, const P4i& _box1, const P4i& _box2, int expansionArea, const dataElem<Scalar>& leader, const dataElem<Scalar>& second, int pid )
{
  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(nSegments);

  // map all elements form the leader to other view copying viewnormals, 

  ////// first project the bbox to get a box into the second image!:
  // this is pushing back to ref image then to RT image !
  //    M3 Hom  = gHom1->getHom( proposal, segId );Hom.invert();Hom = gHom2->getHom( proposal, segId ) * Hom;
  M3 Hom  = second.Hom * leader.iHom;
  M3 iHom = Hom;iHom.invert();

  M3 HomCC  = Hom * Kl;
  M3 iHomCC = iHom * Kl2;

  P3 normal  = leader.iNom;
  P3 vn_0    = second.iNom;

  P4i box1A = _box1;
  P4i box2A = _box2;

  /////////////
  // all those within the bounds are interesting for data term and variables:

  if (pid <0 || newElemV.size() <=0 )
  {
    newElem.clear();
    newElem.seg2varFW.resize(nSegments , -1);
    newElem.seg2varBW.resize(nSegments2, -1);
    newElem.varsFW =0;
    newElem.varsBW =0;

    newElem.firstIdFW = leader.firstIdBW;
    newElem.firstIdBW = second.firstIdBW;
    newElem.boxFW = leader.boxBW;
    newElem.boxBW = second.boxBW;
    box1 = leader.boxBW;
    box2 = second.boxBW;
  }
  else
  {
    newElemV[pid].clear();
    newElemV[pid].seg2varFW.resize(nSegments , -1);
    newElemV[pid].seg2varBW.resize(nSegments2, -1);
    newElemV[pid].varsFW =0;
    newElemV[pid].varsBW =0;

    newElemV[pid].firstIdFW = leader.firstIdBW;
    newElemV[pid].firstIdBW = second.firstIdBW;
    newElemV[pid].boxFW = leader.boxBW;
    newElemV[pid].boxBW = second.boxBW;
  }

  for( int i =0; i< leader.seg2segBW.size(); i++ )
  {
    int my_segId = leader.seg2segBW[i].first;
    if (pid <0 || newElemV.size() <=0 )
    {
      newElem.seg2varFW[my_segId] = leader.seg2varBW[my_segId];

      // NEW _ TEST
      if (newElem.seg2varFW[my_segId] >= 0)
        newElem.varsFW++;
    }
    else
    {
      newElemV[pid].seg2varFW[my_segId] = leader.seg2varBW[my_segId];

      if (newElemV[pid].seg2varFW[my_segId] >= 0)
        newElemV[pid].varsFW++;      
    }

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = HomCC * P3( &centers[3*my_segId] );p0 /= p0[2];
    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId2 = toPixId(p0);
      partner =  segImg2[pixId2];
    }
    segmentsareas1.push_back( std::pair<int,int> (my_segId, partner) );
  }
  assert( pid<0 || newElemV[pid].varsFW == leader.varsBW);
  assert( pid>=0  || newElem.varsFW == leader.varsBW);

  // why not project back all in the box2 and check, later consolidate anyway ?!
  for( int i =0; i< second.seg2segBW.size(); i++ )
  {
    int my_segId = second.seg2segBW[i].first; // each in the map is contained in the box2 by construction.
    if (pid <0 || newElemV.size() <=0 )
      newElem.seg2varBW[my_segId] = second.seg2varBW[my_segId];
    else
      newElemV[pid].seg2varBW[my_segId] = second.seg2varBW[my_segId];

    // single threaded
    if (pid <0 || newElemV.size() <=0 )
    {
      if (newElem.seg2varBW[my_segId] >= 0)
        newElem.varsBW++;
    }
    else
    {
      if (newElemV[pid].seg2varBW[my_segId] >= 0)
        newElemV[pid].varsBW++;      
    }

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = iHomCC * P3( &centers2[3*my_segId] );p0 /= p0[2];

    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId = toPixId(p0);
      partner =  segImg[pixId];
    }

    segmentsareas2.push_back( std::pair<int,int> (my_segId, partner) );
  }

  assert( pid<0   || newElemV[pid].varsBW == second.varsBW);
  assert( pid>=0  || newElem.varsBW == second.varsBW);

  if (pid <0 || newElemV.size() <=0 )
  {
    newElem.seg2segFW = segmentsareas1;
    newElem.seg2segBW = segmentsareas2;
    newElem.proposal = proposal;
    newElem.Hom  = Hom;
    newElem.Nom  = normal;
    newElem.iHom = iHom;
    newElem.iNom = vn_0;
  }
  else
  {
    newElemV[pid].seg2segFW = segmentsareas1;
    newElemV[pid].seg2segBW = segmentsareas2;
    newElemV[pid].proposal = proposal;
    newElemV[pid].Hom  = Hom;
    newElemV[pid].Nom  = normal;
    newElemV[pid].iHom = iHom;
    newElemV[pid].iNom = vn_0;    
  }
}

template<class Scalar>
int Datasegments<Scalar>::
initNewSeg( int propId, int proposal, int boxX, int boxY, genHomoG<Scalar>* gHom, int expansionArea, const std::vector<P3>* normals, int addToSeg2Vars, int addToAllSeg2Vars, int pid )
{
  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(nSegments);

  expansionArea =0;

  P3 p0;
  if ( nExpCenters >0 && expCenters != NULL)
  {
    assert(propId < nExpCenters);
    p0     = Kl * P3( &expCenters[3*propId] );// center is in camera coords
  }
  else
  {
    assert(propId < nSegments);
    p0     = Kl * P3( &centers[3*propId] );// center is in camera coords
  }
  p0 /= p0[2];

  P3 pc_proj(p0);

  int startX = int ( floor(p0[0] -0.5) ) - boxX;
  int endX   = int ( floor(p0[0] -0.5) ) + boxX;
  int startY = int ( floor(p0[1] -0.5) ) - boxY;
  int endY   = int ( floor(p0[1] -0.5) ) + boxY;

  ////// first project the bbox to get a box in the second image:
  M3 Hom  = gHom->getHomC0( proposal );
  M3 iHom = Hom;iHom.invert();

  P3 vn_0  = gHom->getViewNormalC0( proposal ); // no iiKlt2 * 

  p0    = Hom * P3(startX-expansionArea+1, startY-expansionArea+1, 1.);
  P3 p1 = Hom * P3(startX-expansionArea+1, endY+expansionArea+1, 1.);
  p0 /= p0[2];p1 /= p1[2];

  Scalar minX = min(p0[0], p1[0]);
  Scalar maxX = max(p0[0], p1[0]);
  Scalar minY = min(p0[1], p1[1]);
  Scalar maxY = max(p0[1], p1[1]);

  p0 = Hom * P3(endX+expansionArea+1, startY-expansionArea+1, 1.);
  p1 = Hom * P3(endX+expansionArea+1,   endY+expansionArea+1, 1.);
  p0 /= p0[2];p1 /= p1[2];

  minX = min(minX, p0[0]);minX = min(minX, p1[0]);minX = max(int(floor(minX-0.5)), 0);
  minY = min(minY, p0[1]);minY = min(minY, p1[1]);minY = max(int(floor(minY-0.5)), 0);
  maxX = max(maxX, p0[0]);maxX = min(maxX, p1[0]);maxX = min(int(floor(maxX-0.5)), w);
  maxY = max(maxY, p0[1]);maxY = min(maxY, p1[1]);maxY = min(int(floor(maxY-0.5)), h);

  startX =    int(startX/gridSize)*gridSize;
  endX   = (1+int(endX/gridSize))*gridSize;
  startY =    int(startY/gridSize)*gridSize;
  endY   = (1+int(endY/gridSize))*gridSize;

  P4i box1A = P4i (startX, startY, endX, endY);

#ifdef _withoutCuttingImprobable_
  P4i box2A = P4i (minX, minY, maxX, maxY);
#else
  pc_proj = Hom * pc_proj;pc_proj /= pc_proj[2];// -> projected center this time really projected one
  P4i box2A = P4i (pc_proj[0], pc_proj[1], pc_proj[0], pc_proj[1]);

  P4i box2B = P4i (pc_proj[0], pc_proj[1], pc_proj[0], pc_proj[1]);
  box2B.minimize( P4i (w,h,w,h) );// new test: TODO check
  box2B.maximize( P4i (0,0,0,0) );
#endif
  box1A.minimize( P4i (w,h,w,h) );
  box1A.maximize( P4i (0,0,0,0) );
  box2A.minimize( P4i (w,h,w,h) );
  box2A.maximize( P4i (0,0,0,0) );

  /////////////
#ifdef _DEBUG
  Scalar depthOwna    = 1. /(vn_0 | P3( &centers2[3*propId] ));// 2nd view
  Scalar depthOwnb    = 1./ ((*normals)[proposal] | P3( &centers[3*propId] ));// 2nd view
#endif
  // all those within the bounds are interesting for data term and variables:

  if (pid <0 || gs1V.size() <=0 )
  {
    newElem.clear();
    newElem.seg2varFW.resize(nSegments , -1);
    newElem.seg2varBW.resize(nSegments2, -1);
    newElem.varsFW =0;
    newElem.varsBW =0;
  }
  else
  {
    newElemV[pid].clear();
    newElemV[pid].seg2varFW.resize(nSegments , -1);
    newElemV[pid].seg2varBW.resize(nSegments2, -1);
    newElemV[pid].varsFW =0;
    newElemV[pid].varsBW =0;    
  }
  for (int i=0; i < projPatchCenters.size(); i++)
  {
    PatchCenter<Scalar> &pc = projPatchCenters[i];
    // those are of interest:
    if (!(pc.px >= box1A[0] && pc.py >= box1A[1] && pc.px <= box1A[2] && pc.py <= box1A[3]))
      continue;

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = Hom * pc.center;p0 /= p0[2];
    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId2 = toPixId(p0);

      partner =  segImg2[pixId2];
      // check if partner visible:

//      Scalar testD = 1./((inoms_pix [proposal]) | P3(projPatchCenters2[partner].px, projPatchCenters2[partner].py, 1.));

      if ( ((inoms_pix [proposal]) | P3(projPatchCenters2[partner].px, projPatchCenters2[partner].py, 1.)) < 0 )
      {
        box2A[0] = min( box2A[0], projPatchCenters2[partner].px-1);
        box2A[2] = max( box2A[2], projPatchCenters2[partner].px+1);
        box2A[1] = min( box2A[1], projPatchCenters2[partner].py-1);
        box2A[3] = max( box2A[3], projPatchCenters2[partner].py+1);
      }
    }
    segmentsareas1.push_back( std::pair<int,int> (pc.segId, partner) );

    if ( ((noms_pix [proposal]) | pc.center ) < 0 )
    {
      if (pid <0 || gs1V.size() <=0 )
        newElem.seg2varFW[pc.segId] = newElem.varsFW++ + addToAllSeg2Vars;
      else
        newElemV[pid].seg2varFW[pc.segId] = newElemV[pid].varsFW++ + addToAllSeg2Vars;
    }
  }

  // restrict maximum expansion size to .. 
  // TODO also consider depth here - e.g. parts in front of ?
  minX=box2A[0];    minY=box2A[1];maxX=box2A[2];    maxY=box2A[3];

  if ( (maxX-minX)*(maxY-minY) > 0 && (maxX-minX)*(maxY-minY) > (boxX*2+1)*(boxY*2+1)*resizer*resizer )
  {
    /// center maps to :
#ifdef _writeDebugOut_
    printf("Fixed Box: before %.1f>%.1f  (minX,maxX,..)  %.1f,%.1f,%.1f,%.1f ", (maxX-minX)*(maxY-minY), (boxX*2+1)*(boxY*2+1)*resizer*resizer, minX, maxX, minY, maxY );
#endif

    // 1.st plot depth's of 4 points:
    //Scalar dc  = 1./(noms_pix [proposal]  | P3(pc_proj[0], pc_proj[1], 1.));

#ifdef _withoutCuttingImprobable_
    // in quick hack already done
    pc_proj = Hom * pc_proj;pc_proj /= pc_proj[2];// -> projected center this time really projected one
#endif

    minX = max (minX, (pc_proj[0]) - resizer*boxX);
    maxX = min (maxX, (pc_proj[0]) + resizer*boxX);
    minY = max (minY, (pc_proj[1]) - resizer*boxY);
    maxY = min (maxY, (pc_proj[1]) + resizer*boxY);
    maxX = max(maxX , minX );minX = min(maxX , minX );
    maxY = max(maxY , minY );minY = min(maxY , minY );

#ifdef _writeDebugOut_
    printf("After: %.1f. center:%.1f  %.1f | (minX,maxX)  %.1f,%.1f,%.1f,%.1f\n", (maxX-minX)*(maxY-minY), pc_proj[0], pc_proj[1], minX, maxX, minY, maxY);
#endif
    //mexEvalString("drawnow");
  }
  minX =    floor(minX/gridSize) *gridSize;
  maxX = (1+floor(maxX/gridSize))*gridSize;
  minY =    floor(minY/gridSize) *gridSize;
  maxY = (1+floor(maxY/gridSize))*gridSize;
  box2A = P4i (minX, minY, maxX, maxY);

  // new:
  box2A.minimize( P4i (w,h,w,h) );
  box2A.maximize( P4i (0,0,0,0) ); 


  for (int i=0; i < projPatchCenters2.size(); i++)
  {
    PatchCenter<Scalar> pc = projPatchCenters2[i];
    // those in the box2 are looked at here 
    if (!(pc.px >= box2A[0] && pc.py >= box2A[1] && pc.px <= box2A[2] && pc.py <= box2A[3]))
      continue;

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = iHom * pc.center;p0 /= p0[2];

    int lol=0;
    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId = toPixId(p0);
      partner =  segImg[pixId];
    }

    if ( (pid <0 || gs1V.size() <=0) && partner != -1 && newElem.seg2varFW[partner] <0 ) continue;// could lead to strong shrinkage! e.g. oob stuff
    if (!(pid <0 || gs1V.size() <=0) && partner != -1 && newElemV[pid].seg2varFW[partner] <0 ) continue;// could lead to strong shrinkage! e.g. oob stuff

    // oob can be possible?! blocking not so bad though, keeps current solution if good: fine
#ifdef _FW_ONLY_
    continue;
#endif
    segmentsareas2.push_back( std::pair<int,int> (pc.segId, partner) );

    // TODO: ATTENTION: here the test compares canonical view to other - also further away views ! 
    // these must be treated differently as in theory displacements can accumulate; so 2* disp, etc: done 
    // in front of cam -> non-sense to add a variable: (could actually be even stricter: 1 as 1 means super close)
    if ( ((inoms_pix [proposal]) | pc.center ) < 0 && (p0-pc.center).sqrnorm() < maxMotSqr ) // only here since otherwise different var sets are picked
    {
      if (pid <0 || gs1V.size() <=0 )
        newElem.seg2varBW[pc.segId] = newElem.varsBW++  +  newElem.varsFW + addToSeg2Vars + addToAllSeg2Vars;
      else
        newElemV[pid].seg2varBW[pc.segId] = newElemV[pid].varsBW++  +  newElemV[pid].varsFW + addToSeg2Vars + addToAllSeg2Vars;
#ifdef _test2_
      // TODO: tighter box ?? Could also reduce box size here: all that pass make up the box2A - nothing more!
      box2B[0] = min( box2B[0], pc.px-1);
      box2B[2] = max( box2B[2], pc.px+1);
      box2B[1] = min( box2B[1], pc.py-1);
      box2B[3] = max( box2B[3], pc.py+1);
#endif
    }
  }

#ifdef _test2_
  // is the reason for smaller energy: omitting this before
  box2B[0] =    floor(Scalar(box2B[0]/gridSize)) *gridSize;
  box2B[2] = (1+floor(Scalar(box2B[2]/gridSize)))*gridSize;
  box2B[1] =    floor(Scalar(box2B[1]/gridSize)) *gridSize;
  box2B[3] = (1+floor(Scalar(box2B[3]/gridSize)))*gridSize;
  box2A = box2B;
#endif

  box2A.minimize( P4i (w,h,w,h) );
  box2A.maximize( P4i (0,0,0,0) );
  box1A.minimize( P4i (w,h,w,h) );
  box1A.maximize( P4i (0,0,0,0) );

  if (pid <0 || gs1V.size() <=0 )
  {
    newElem.seg2segFW = segmentsareas1;
    newElem.seg2segBW = segmentsareas2;
    newElem.proposal = proposal;
    newElem.Hom  = Hom;
    newElem.Nom  = iKlt * (*normals)[proposal];
    newElem.iHom = iHom;
    newElem.iNom = iKlt2 * vn_0;
    newElem.firstIdFW = 0;
    newElem.firstIdBW = newElem.varsFW + addToSeg2Vars;
    newElem.boxFW = box1A;
    newElem.boxBW = box2A;
    box1 = box1A;
    box2 = box2A;
  }
  else
  {
    newElemV[pid].seg2segFW = segmentsareas1;
    newElemV[pid].seg2segBW = segmentsareas2;
    newElemV[pid].proposal = proposal;
    newElemV[pid].Hom  = Hom;
    newElemV[pid].Nom  = iKlt * (*normals)[proposal];
    newElemV[pid].iHom = iHom;
    newElemV[pid].iNom = iKlt2 * vn_0;
    newElemV[pid].firstIdFW = 0;
    newElemV[pid].firstIdBW = newElemV[pid].varsFW + addToSeg2Vars;
    newElemV[pid].boxFW = box1A;
    newElemV[pid].boxBW = box2A;
    box1 = box1A;
    box2 = box2A;
  }
  // new last variable number
  if (pid <0 || gs1V.size() <=0 )
    return newElem.varsBW + addToSeg2Vars;
  else
    return newElemV[pid].varsBW + addToSeg2Vars;
}

template<class Scalar>
void
  Datasegments<Scalar>::
  requestNElements( int nProcs, int nProps )
{
  for (int i = 0;i<nProcs;i++)
  {
    gws1V.push_back( new genWarpS<Scalar>( w,h, gs1.getRefImgPtr() ) ); 
    gws2V.push_back( new genWarpS<Scalar>( w,h, gs2.getRefImgPtr() ) );
    gs1V.push_back( new genScore<Scalar> ( w,h, gs1.getDataThresh(), gs1.getOobThresh() ) );
    gs2V.push_back( new genScore<Scalar> ( w,h, gs2.getDataThresh(), gs2.getOobThresh() ) );
    gs1V[i]->setRefImage( gs1.getRefImgPtr() );
    gs2V[i]->setRefImage( gs2.getRefImgPtr() );

    gs1V[i]->setMaxDisp( gs1.getMaxDisp() );
    gs2V[i]->setMaxDisp( gs2.getMaxDisp() );
    gs1V[i]->setMaxMot( gs1.getMaxMot() );
    gs2V[i]->setMaxMot( gs2.getMaxMot() );

    gs1V[i]->setAutoScores( gs1.getAutoScores() ) ;
    gs2V[i]->setAutoScores( gs2.getAutoScores() );
  }
  newElemV.resize(nProcs);
  dataElements.resize( nProps );
  initPropStack( nProcs );
}

template<class Scalar>
void
  Datasegments<Scalar>::
  releaseNElements()
{
  for (int i = 0;i<gws1V.size();i++)
  { if (gws1V[i] != NULL) delete gws1V[i];gws1V[i] = NULL; }
  for (int i = 0;i<gws2V.size();i++)
  { if (gws2V[i] != NULL) delete gws2V[i];gws2V[i] = NULL; }
  for (int i = 0;i<gs1V.size();i++)
  { if (gs1V[i] != NULL) delete gs1V[i];gs1V[i] = NULL; }
  for (int i = 0;i<gs2V.size();i++)
  { if (gs2V[i] != NULL) delete gs2V[i];gs2V[i] = NULL; }

  gws1V.clear();
  gws2V.clear();
  gs1V.clear();
  gs2V.clear();
  newElemV.clear();
}

template<class Scalar>
void
  Datasegments<Scalar>::
  createWarps( int pid, int propId )
{
  assert (pid < newElemV.size() );
  if (pid <0)
    createWarps ( box1, newElem.Hom, box2, newElem.iHom );
  else
    createWarps ( newElemV[pid].boxFW, newElemV[pid].Hom, newElemV[pid].boxBW, newElemV[pid].iHom, pid, propId );
}

template<class Scalar>
void Datasegments<Scalar>::
  createWarps( const P4i& bigbox1, const M3& Hom, const P4i& bigbox2, const M3& iHom, int pid, int propId )
{
  if (pid <0 || pid > newElemV.size() )
  {
    gws1.warp_noOmp_patchBased( bigbox2, iHom );// warps image 1:
    gws2.warp_noOmp_patchBased( bigbox1,  Hom );// warps image 2:

    // gws1 receives left/ orig image. 
    gs1.computeFullScoresCensus3OMP_Box( gws1.getOrigShortImage(), gws2.getWarpedShortImage(), gws2.getIdx(), gws2.getIdy(), 
      nSegments, segImg, box1[0], box1[1], box1[2], box1[3], newElem.seg2varFW, occlusions1 );
    gs2.computeFullScoresCensus3OMP_Box( gws2.getOrigShortImage(), gws1.getWarpedShortImage(), gws1.getIdx(), gws1.getIdy(), 
      nSegments2, segImg2, box2[0], box2[1], box2[2], box2[3], newElem.seg2varBW, occlusions2 );

    // append and create full description of this exansion:
    std::vector<Scalar> resD = gs1.getScores();
    std::vector<Scalar> resF = gs1.getFreeScores();
    std::vector<int> resV    = gs1.getFreeVariables();

    for(int i=0; i<(newElem.seg2segFW).size(); i++ )
    {
#ifdef __depthControl__
      int sid = (newElem.seg2segFW)[i].first;
      Scalar depth0      = 1./(Centers[ sid ]|newElem.Nom); 
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas1[sid] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresFW ? adaptiveAutoScoresFW[sid] : 1.)); // special case segment behind cam -- thus invisible
        newElem.dataFW.push_back( impPen );
        newElem.freedataFW.push_back( 0 );
        newElem.freeFW.push_back( 0 );
      }
      else
#endif
      {
        newElem.dataFW.push_back( resD[ (newElem.seg2segFW)[i].first ] );
        newElem.freedataFW.push_back( resF[ (newElem.seg2segFW)[i].first ] );
        newElem.freeFW.push_back( resV[ (newElem.seg2segFW)[i].first ] );
      }
    }
    resD = gs2.getScores();
    resF = gs2.getFreeScores();
    resV = gs2.getFreeVariables();
    for(int i=0; i<(newElem.seg2segBW).size(); i++ )
    {
#ifdef __depthControl__
      int sid = (newElem.seg2segBW)[i].first;
      Scalar depth0  = 1./(iCenters[ sid ]|newElem.iNom);
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas2[ sid ] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresBW ? adaptiveAutoScoresBW[sid] : 1.)); // special case segment behind cam -- thus invisible
        newElem.dataBW.push_back( impPen );
        newElem.freedataBW.push_back( 0 );
        newElem.freeBW.push_back( 0 );
      }
      else
#endif
      {
        newElem.dataBW.push_back( resD[ (newElem.seg2segBW)[i].first ] );
        newElem.freedataBW.push_back( resF[ (newElem.seg2segBW)[i].first ] );
        newElem.freeBW.push_back( resV[ (newElem.seg2segBW)[i].first ] );
      }
    }

    dataElements.push_back( newElem );
  }
  else
  {
    // if -1 and no neigh is +1 do no evaluate newElemV[pid].seg2varFW[pc.segId] warp is per pixel u know ..
    // per pixel check might be more costly than doing nothing
    gws1V[pid]->warp_noOmp_patchBased( bigbox2, iHom );// warps image 1:
    gws2V[pid]->warp_noOmp_patchBased( bigbox1,  Hom );// warps image 2:

    // gws1 receives left/ orig image. 
    gs1V[pid]->computeFullScoresCensus3OMP_Box( gws1V[pid]->getOrigShortImage(), gws2V[pid]->getWarpedShortImage(), gws2V[pid]->getIdx(), gws2V[pid]->getIdy(), 
      nSegments, segImg, newElemV[pid].boxFW[0], newElemV[pid].boxFW[1], newElemV[pid].boxFW[2], newElemV[pid].boxFW[3], 
      newElemV[pid].seg2varFW, occlusions1 );
    gs2V[pid]->computeFullScoresCensus3OMP_Box( gws2V[pid]->getOrigShortImage(), gws1V[pid]->getWarpedShortImage(), gws1V[pid]->getIdx(), gws1V[pid]->getIdy(), 
      nSegments2, segImg2, newElemV[pid].boxBW[0], newElemV[pid].boxBW[1], newElemV[pid].boxBW[2], newElemV[pid].boxBW[3],
      newElemV[pid].seg2varBW, occlusions2 );

    // append and create full description of this exansion:
    std::vector<Scalar> resD = gs1V[pid]->getScores();
    std::vector<Scalar> resF = gs1V[pid]->getFreeScores();
    std::vector<int> resV    = gs1V[pid]->getFreeVariables();

    for(int i=0; i<(newElemV[pid].seg2segFW).size(); i++ )
    {
#ifdef __depthControl__
      int sid = (newElemV[pid].seg2segFW)[i].first;
//      Scalar depth0      = 1./(Centers[ sid ]|newElemV[pid].Nom); 
      Scalar depth0  = 1./((Centers[ sid ])|(noms_pix[newElemV[pid].proposal]));
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas1[ sid ] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresFW ? adaptiveAutoScoresFW[sid] : 1.)); // special case segment behind cam -- thus invisible
        newElemV[pid].dataFW.push_back( impPen );
        newElemV[pid].freedataFW.push_back( 0 );
        newElemV[pid].freeFW.push_back( 0 );
      }
      else
#endif
      {
        newElemV[pid].dataFW.push_back(     resD[ (newElemV[pid].seg2segFW)[i].first ] );
        newElemV[pid].freedataFW.push_back( resF[ (newElemV[pid].seg2segFW)[i].first ] );
        newElemV[pid].freeFW.push_back(     resV[ (newElemV[pid].seg2segFW)[i].first ] );
      }
    }
    resD = gs2V[pid]->getScores();
    resF = gs2V[pid]->getFreeScores();
    resV = gs2V[pid]->getFreeVariables();
    for(int i=0; i<(newElemV[pid].seg2segBW).size(); i++ )
    {
#ifdef __depthControl__
      int sid = (newElemV[pid].seg2segBW)[i].first;
//      Scalar depth0  = 1./(iCenters[ sid ]|newElemV[pid].iNom);
      Scalar depth0  = 1./((iCenters[ sid ])|(inoms_pix[newElemV[pid].proposal]));
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas2[ sid ] * Scalar(__depthControlAmplifier__ * impFactor);//*(xtraPen+ (autoScoresBW ? adaptiveAutoScoresBW[sid] : 1.)); // special case segment behind cam -- thus invisible
        newElemV[pid].dataBW.push_back( impPen );
        newElemV[pid].freedataBW.push_back( 0 );
        newElemV[pid].freeBW.push_back( 0 );
      }
      else
#endif
      {
        newElemV[pid].dataBW.push_back    ( resD[ (newElemV[pid].seg2segBW)[i].first ] );
        newElemV[pid].freedataBW.push_back( resF[ (newElemV[pid].seg2segBW)[i].first ] );
        newElemV[pid].freeBW.push_back    ( resV[ (newElemV[pid].seg2segBW)[i].first ] );
      }
    }
    dataElements[propId] = newElemV[pid];// not in order need to resize dataelements before + have propid as info
  }
}

template<class Scalar>
void 
Datasegments<Scalar>::
averageDataScore()
{
  if (dataElements.size() <=0) return;

  std::vector<int> allFreeFW   ( dataElements[0].seg2varFW.size(), 0 );
  std::vector<Scalar> allScoreFW  ( dataElements[0].seg2varFW.size(), 0. );
  std::vector<int> allFreeBW   ( dataElements[0].seg2varBW.size(), 0 );
  std::vector<Scalar> allScoreBW  ( dataElements[0].seg2varBW.size(), 0. );

  for(int i=0; i< dataElements.size();i++)
  {
      for(int j=0; j< dataElements[i].seg2segFW.size();j++)
      {
        allFreeFW[ dataElements[i].seg2segFW[j].first ] += dataElements[i].freeFW[j];
        allScoreFW[ dataElements[i].seg2segFW[j].first ] += dataElements[i].freedataFW[j];
      }
      for(int j=0; j< dataElements[i].seg2segBW.size();j++)
      {
        allFreeBW[ dataElements[i].seg2segBW[j].first ] += dataElements[i].freeBW[j];
        allScoreBW[ dataElements[i].seg2segBW[j].first ] += dataElements[i].freedataBW[j];
      }
  }

  // per free pixel (pixel in segment with correspondence, not oob)
  std::vector<Scalar> avScoreFW  (dataElements[0].seg2varFW.size(), 0. );
  std::vector<Scalar> avScoreBW  (dataElements[0].seg2varFW.size(), 0. );

  for(int i=0; i< avScoreFW.size();i++)
    avScoreFW[i] = allScoreFW[i] / max( 1, allFreeFW[i] );
  for(int i=0; i< avScoreBW.size();i++)
    avScoreBW[i] = allScoreBW[i] / max( 1, allFreeBW[i] );
  ////
  if ( adaptiveAutoScoresFW.size() != 0 )
  {
    for (int i=0;i< adaptiveAutoScoresFW.size(); i++)
        if ( allFreeFW[i] > 1 ) 
          adaptiveAutoScoresFW[i] = avScoreFW[i];// min( adaptiveAutoScoresFW[i], avScoreFW[i]);

      // dataElements are given per proposal ! not per segment
      for (int i=0;i< dataElements.size(); i++)
        dataElements[i].updateOobPenaltyFW ( adaptiveAutoScoresFW, areas1 );
  }

  if ( adaptiveAutoScoresBW.size() != 0 )
  {
    for (int i=0;i< adaptiveAutoScoresBW.size(); i++)
        if ( allFreeBW[i] > 1 ) 
          adaptiveAutoScoresBW[i] = avScoreBW[i];//min (adaptiveAutoScoresBW[i], avScoreBW[i]);

      for (int i=0;i< dataElements.size(); i++)
        dataElements[i].updateOobPenaltyBW ( adaptiveAutoScoresBW, areas2 );
  }
}


template<class Scalar>
Scalar Datasegments<Scalar>::
  getEnergy( ) 
{
  int dummy;
  Scalar energy(0);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < nSegments; i++ )
  {
    int segId = i;

    // currentsolution is not to lookup in the dataElements array !! 
    // currentSolution == normal/rotation; dataElement == normal/rot + center combiniation -> not the same amount
    int prop0   = currentSolution1[ segId ];
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// projects onto

#ifdef _simple_data_
    seg2_0 = -1;
#endif

    if (seg2_0 <0) // maps on no other segment:
      energy +=  currentSolutionE.dataFW[ segId ];
    else
    {
      int prop00 = currentSolution2[seg2_0];

      if( prop00 == prop0 ) // normal data term
        energy += currentSolutionE.dataFW[ segId ];
      else // must check depth in both versions:
        energy += getPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ segId ], areas1[ segId ],
        currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
        dummy, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. ) + areas1[ segId ] * addPotts;
    }
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < nSegments2; i++ )
  {
    int segId = i;

    int prop0   = currentSolution2[ segId ];
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment
#ifdef _simple_data_
    seg2_0 = -1;
#endif
    if (seg2_0 < 0 )
      energy += currentSolutionE.dataBW[ segId ];
    else
    {
      int prop00 = currentSolution1[seg2_0];

      if( prop00 == prop0 ) // normal data term
        energy += currentSolutionE.dataBW[ segId ];
      else // must check depth in both versions:
        energy += getPenalty ( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ segId ], areas2[ segId ],
        currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
        dummy, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. ) + areas2[ segId ] * addPotts;
    }///////////////////////////////////////////
  }/////////////////
#endif
  return energy;
}


template<class Scalar>
void  Datasegments<Scalar>::
getEnergyMap( std::vector<int> &energyMap, std::vector<int> &energyMap2 )
{
  energyMap.resize(nSegments,-1);
  energyMap2.resize(nSegments,-1);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < nSegments; i++ )
  {
    int segId = i;

    // currentsolution is not to lookup in the dataElements array !! 
    // currentSolution == normal/rotation; dataElement == normal/rot + center combiniation -> not the same amount
    int prop0   = currentSolution1[ segId ];
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// projects onto

    if (seg2_0 <0) // maps on no other segment:
      energyMap[i] = -1;
    else
    {
      int prop00 = currentSolution2[seg2_0];

      if( prop00 == prop0 ) // normal data term
        energyMap[i] = 0;
      else // must check depth in both versions:
        energyMap[i] = whichPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ segId ], ihoms_pix[prop00], iCenters[ seg2_0 ], prop00, prop0 );
    }
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < nSegments2; i++ )
  {
    int segId = i;

    int prop0   = currentSolution2[ segId ];
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment

    if (seg2_0 < 0 )
      energyMap2[i] = -1;
    else
    {
      int prop00 = currentSolution1[seg2_0];

      if( prop00 == prop0 ) // normal data term
        energyMap2[i] = 0;
      else // must check depth in both versions:
        energyMap2[i] = whichPenalty ( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ segId ], homs_pix[prop00], Centers[ seg2_0 ], prop00, prop0 );
    }///////////////////////////////////////////
  }/////////////////
#endif
}

template<class Scalar>
void
  Datasegments<Scalar>::
  getDataEnergyMap( std::vector<Scalar> &energyMap, std::vector<Scalar> &energyMap2 )
{
  energyMap.resize(nSegments,-1);
  energyMap2.resize(nSegments,-1);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < nSegments; i++ )
  {
    int segId = i;

    int prop0   = currentSolution1[ segId ];
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// projects onto

    if (seg2_0 < 0 )
      energyMap[i] = currentSolutionE.dataFW[ segId ];
    else
    {
      int prop00 = currentSolution2[seg2_0];
      int dummy=0;

      if( prop00 == prop0 ) // normal data term
        energyMap[i] = currentSolutionE.dataFW[ segId ];
      else
        energyMap[i] = getPenalty( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ segId ], areas1[ segId ],
        currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
        dummy, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. ) + areas1[ segId ] * addPotts;
    }

  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < nSegments2; i++ )
  {
    int segId = i;

    int prop0   = currentSolution2[ segId ];
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment

    if (seg2_0 < 0 )
      energyMap2[i] = currentSolutionE.dataBW[ segId ];
    else
    {
      int prop00 = currentSolution1[seg2_0];
      int dummy=0;

      if( prop00 == prop0 ) // normal data term
        energyMap2[i] = currentSolutionE.dataBW[ segId ];
      else
        energyMap2[i] = getPenalty( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ segId ], areas2[ segId ],
        currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
        dummy, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. ) + areas2[ segId ] * addPotts;

    }

  }/////////////////
#endif
}


template<class Scalar>
void
  Datasegments<Scalar>::
  buildFromCurrentSolution( const std::vector<int>& _currentSolution, const std::vector<int> &_currentSolution2 )
{
  // mapping each segment to a proposal
  currentSolution1 = _currentSolution;
  currentSolution2 = _currentSolution2;
  currentSolutionE.clear();

  std::vector<M3> _homs( nSegments);
  std::vector<M3> _ihoms(nSegments2);

  std::vector<int> segmentInvolved1( nSegments, 1);// all are involved
  std::vector<int> segmentInvolved2( nSegments2, 1);// all are involved

  std::vector< std::pair<int,int> > segmentsareas1; segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2; segmentsareas2.reserve(nSegments);
  /////////////

  // all those within the bounds are interesting for data term and variables:
  for (int i=0; i < currentSolution1.size(); i++)
  {
    PatchCenter<Scalar> pc(0,0, P3(0,0,0), i );
    int proposal = currentSolution1[i];
    M3 Hom  = homs_pix[proposal];// dataElements[proposal].Hom;
    M3 HomCC  = Hom * Kl;

    _homs[i] = Hom;

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = HomCC * P3( &centers[3*pc.segId] );p0 /= p0[2];
    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId2 = toPixId(p0);
      partner =  segImg2[pixId2];
    }
    segmentsareas1.push_back( std::pair<int,int> (pc.segId, partner) );
  }
  for (int i=0; i < currentSolution2.size(); i++)
  {
    PatchCenter<Scalar> pc(0,0, P3(0,0,0), i );
    int proposal = currentSolution2[i];

    // should be equal!
    M3 iHom   = ihoms_pix[proposal]; 
    M3 iHomCC = iHom * Kl2;
    _ihoms[i] = iHom;

    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = iHomCC * P3( &centers2[3*pc.segId] );p0 /= p0[2];
    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      int pixId = toPixId(p0);
      partner =  segImg[pixId];
    }

    segmentsareas2.push_back( std::pair<int,int> (pc.segId, partner) );
  }
  currentSolutionE.clear();
  currentSolutionE.seg2segFW = segmentsareas1;// !!! this might be the problem - check it at the end
  currentSolutionE.seg2segBW = segmentsareas2;
  currentSolutionE.proposal = -1;

  gws1.warp_noOmp_patchBased( _ihoms, segImg2, w,h );
  gws2.warp_noOmp_patchBased(  _homs, segImg,  w,h );

  Scalar sumError =0;Scalar sumVError =0;int sumFError =0;
  gs1.computeFullScoresCensus3OMP_Box( gws1.getOrigShortImage(), gws2.getWarpedShortImage(), gws2.getIdx(), gws2.getIdy(), 
    nSegments, segImg, 0, 0, w, h, segmentInvolved1, occlusions1 );
  gs2.computeFullScoresCensus3OMP_Box( gws2.getOrigShortImage(), gws1.getWarpedShortImage(), gws1.getIdx(), gws1.getIdy(), 
    nSegments2, segImg2, 0, 0, w, h, segmentInvolved2, occlusions2 );

  // append and create full description of this exansion:
  std::vector<Scalar> resD = gs1.getScores();
  std::vector<Scalar> resF = gs1.getFreeScores();
  std::vector<int> resV    = gs1.getFreeVariables();

  for(int i=0; i<(currentSolutionE.seg2segFW).size(); i++ )
  {
#ifdef __depthControl__
      int sid  = (currentSolutionE.seg2segFW)[i].first;
      int pid = currentSolution1[sid];
//      Scalar depth0  = 1./((Kl*Centers[ sid ])|noms_pix[pid]);
      Scalar depth0  = 1./((Centers[ sid ])|noms_pix[pid]);
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas1[ sid ] * Scalar(__depthControlAmplifier__ * impFactor);
        resD[ sid ] = impPen;
        resV[ sid ] = 0;
        resF[ sid ] = 0;
      }
#endif
    sumError  += resD[ (currentSolutionE.seg2segFW)[i].first ];
    sumFError += resV[ (currentSolutionE.seg2segFW)[i].first ];
    sumVError += resF[ (currentSolutionE.seg2segFW)[i].first ];
    currentSolutionE.dataFW.push_back( resD[ (currentSolutionE.seg2segFW)[i].first ] );
    currentSolutionE.freedataFW.push_back( resF[ (currentSolutionE.seg2segFW)[i].first ] );
    currentSolutionE.freeFW.push_back( resV[ (currentSolutionE.seg2segFW)[i].first ] );
  }

#ifdef _writeDebugOut_
  printf("--- Data-Error %f - %f, %d pix free\n", sumError, sumVError, sumFError);
#endif
  sumError =0;sumVError =0;sumFError =0;
  resD = gs2.getScores();
  resF = gs2.getFreeScores();
  resV = gs2.getFreeVariables();
  for(int i=0; i<(currentSolutionE.seg2segBW).size(); i++ )
  {
#ifdef __depthControl__
      int sid = (currentSolutionE.seg2segBW)[i].first;
      int pid = currentSolution2[sid];
//      Scalar depth0  = 1./((Kl2*iCenters[ sid ])|inoms_pix[pid]);
      Scalar depth0  = 1./((iCenters[ sid ])|inoms_pix[pid]);
      if ( depth0 > Scalar(-_minDepth_) )
      {
        Scalar impPen = areas2[ sid ] * Scalar(__depthControlAmplifier__ * impFactor);
        resD[ sid ] = impPen;
        resV[ sid ] = 0;
        resF[ sid ] = 0;
      }
#endif

    sumError += resD[ (currentSolutionE.seg2segBW)[i].first ];sumFError += resV[ (currentSolutionE.seg2segBW)[i].first ];sumVError += resF[ (currentSolutionE.seg2segBW)[i].first ];
    currentSolutionE.dataBW.push_back( resD[ (currentSolutionE.seg2segBW)[i].first ] );
    currentSolutionE.freedataBW.push_back( resF[ (currentSolutionE.seg2segBW)[i].first ] );
    currentSolutionE.freeBW.push_back( resV[ (currentSolutionE.seg2segBW)[i].first ] );
  }
#ifdef _writeDebugOut_
  printf("--- Data-Error v2 %f - %f, %d pix free\n", sumError, sumVError, sumFError);
#endif
  //      correctInFrontOfCam();
}


/// delivers the data term as edges and unaries, with the variable names given in the expansion set
template<class Scalar>
void
  Datasegments<Scalar>::
  expandRegion( int propId ) //, genHomoG<Scalar>* gHom )
{
  const dataElem<Scalar> &temp = dataElements[propId];
  int trial = temp.proposal;

#ifdef _DEBUG
  if ( inoms_pix[trial] != temp.iNom || noms_pix[trial] != temp.Nom )
  {
    printf("Mismatch in noms\n");        mexEvalString("drawnow");
  }
  if ( ihoms_pix[trial]*P3(1.,1.,1.) != temp.iHom*P3(1.,1.,1.) || homs_pix[trial]*P3(1.,1.,1.) != temp.Hom*P3(1.,1.,1.) )
  {
    printf("Mismatch in homs\n");        mexEvalString("drawnow");
  }
#endif

  unaries.clear();
  binaries.clear();
  //    unaries.assign ( temp.varsFW + temp.varsBW, P2(0,0) );  //    unaries.assign ( temp.dataFW.size() + temp.dataBW.size(), P2(0,0) );  
  binaries.reserve( 4*(temp.seg2segFW.size() + temp.seg2segBW.size()) );//, Binary<Scalar> (-1,-1));
  unaries.reserve(  2*(temp.seg2segFW.size() + temp.seg2segBW.size()) );//, could be FW,BW and 0,1

  // used to ensure that each binary term produced is unique - no doubles so we can skip merging edges !
  std::vector< std::pair<int,int> > edgeStoreMap (temp.varsFW, std::pair<int,int> (-1,-1));
  std::vector< std::pair<int,int> > edgeStoreMapB(temp.varsFW, std::pair<int,int> (-1,-1));

  motionPairsFW.clear();
  motionPairsBW.clear();
  if (computeMotionPairs)
  {
    motionPairsFW.reserve(4*(temp.seg2segFW.size() + temp.seg2segBW.size()));
    motionPairsBW.reserve(4*(temp.seg2segFW.size() + temp.seg2segBW.size()));
  }
  //motionPairs.push_back( MotionPairs( segI,segJ, varId1, varid2, v1=0/1 ) )
  ////////////////

  // part one: left to right:
  const std::vector< std::pair<int,int> > &fwMap = temp.seg2segFW;

  for(int i= 0; i < fwMap.size(); i++ )
  {
    int edgeCase1 = 0;
    int edgeCase2 = 0;

    // 1. assign a label 0 to segi in left view:
    int segId = fwMap[i].first;
    if (segId<0) continue;
    // so. check onto which segment the segment projects to and check if it is considered
    // how? here:
    int prop0   = currentSolution1[ segId ];

    int prop1   = trial;
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// lbl 0 in view1 maps to this segment
    int seg2_1  = fwMap[i].second;                           // lbl 1 in view1 maps to this segment
    int varId0x = -2;//does not project to anything. 
    int varId1x = -2;
    int varIdx  = temp.seg2varFW[ segId ];// id of the first variable

    // NEW _ TEST
    if ( varIdx < 0) continue;

    if (seg2_0 > -1 )
      varId0x = temp.seg2varBW[ seg2_0 ];// variable id (can be -1 == not considered)
    if (seg2_1 > -1 )
      varId1x = temp.seg2varBW[ seg2_1 ];// should be >=0
    // these are the local variables. 
    //      varId0x =-2;
    //      varId1x =-2;
    // first addres case 00 and 01:

    if ( prop0 == trial && (varId0x<=-1 || currentSolution2[seg2_0] == trial) ) continue; // no impact 

    // then binary, else 00 is a unary, 01 does not exist, no matter what 
    if (varId0x == -2) // maps on no other segment:
      unaries.push_back( Unary<Scalar>(varIdx, P2( currentSolutionE.dataFW[ segId ], 0) )) ;
    else // segment exists but no variable for it:, so 01 does not exist
    {
      // only need to consider 00:
      // same segi? or dame depth? or larger or smaller -> 3 possibilities

      int prop00 = currentSolution2[seg2_0];
      Scalar penalty =0;
      if( prop00 == prop0 ) // normal data term
        penalty = currentSolutionE.dataFW[ segId ];
      else // must check depth in both versions:
      {
        penalty = getPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ segId ], areas1[ segId ],
          currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
          edgeCase1, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
        penalty += addPotts * areas1[ segId ];
      }
      if ( varId0x == -1 )//other segment has no partner
        unaries.push_back( Unary<Scalar>(varIdx, P2( penalty, 0) ));
      else
      {
        Scalar pen00 = penalty;

        int prop01 = trial;
        if( prop01 == prop0 ) // normal data term
          penalty = currentSolutionE.dataFW[ segId ];
        else // must check depth in both versions:
        {
          penalty = getPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop01], Centers[ segId ], areas1[ segId ],
            currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
            edgeCase1, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
          penalty += addPotts * areas1[ segId ];
        }
        ////
        binaries.push_back( Binary<Scalar> (varIdx, varId0x, P4 (pen00, penalty ,0,0) ) );
        edgeStoreMap[(varIdx - temp.firstIdFW)] = std::pair<int,int> (varId0x, binaries.size()-1);

        // same for 01 case, an edge exists:

      }
    }///////////////////////////////////////////
    // case 10, 11
    if (varId1x == -2) // maps on no other segment - oob
      unaries.push_back( Unary<Scalar>( varIdx, P2( 0, temp.dataFW[ i ]) ) );
    else // segment exists but no variable for it:, so 01 does not exist
    {
      Scalar penalty =0;
      // only need to consider 10:
      // same segi? or dame depth? or larger or smaller -> 3 possibilities
      int prop10 = currentSolution2[seg2_1];
      if( prop10 == trial ) // normal data term
        penalty = temp.dataFW[ i ];
      else // must check depth in both versions:
      {
        penalty = getPenalty ( temp.Hom, temp.iNom, inoms_pix[prop10], Centers[ segId ], areas1[ segId ],
          temp.dataFW[ i ], temp.freedataFW[ i ], temp.freeFW[ i ], 
          edgeCase2, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
        penalty += addPotts * areas1[ segId ];
      }
      if ( varId1x == -1 ) // not encoded variable
        unaries.push_back( Unary<Scalar>(varIdx, P2( 0, penalty ) ));
      //          unaries[ varIdx ] += P2( 0, penalty );
      else  // is encoded
      {
        Scalar pen10 = penalty;

        int prop11 = trial;
        if( prop11 == prop1 ) // normal data term
          penalty = temp.dataFW[ i ];//currentSolutionE.dataFW[ segId ];
        else // can not enter this case: -- NONSENSE FROM HERE
        {
          assert( 0 );
          // get depth relationship:
          penalty = getPenalty ( temp.Hom, temp.iNom, inoms_pix[prop11], Centers[ segId ], areas1[ segId ],
            temp.dataFW[ i ], temp.freedataFW[ i ], temp.freeFW[ i ], 
            edgeCase2, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
          penalty += addPotts * areas1[ segId ];
        }
        ////
        if (varId0x != varId1x)
        {
          binaries.push_back( Binary<Scalar> (varIdx, varId1x, P4 (0, 0, pen10, penalty ) ) );
          edgeStoreMapB[varIdx - temp.firstIdFW] = std::pair<int,int> (varId1x, binaries.size()-1);
        }
        else // double entry:
        {
          binaries[binaries.size()-1].edge += P4 (0, 0, pen10, penalty );
          //edgeStoreMap[(varIdx - temp.firstIdFW)] = std::pair<int,int> (varId1x, binaries.size()-1);          
        }
        // same for 01 case, an edge exists:
      }
    }///////////////////////////////////////////

    if (computeMotionPairs)
    {
      // check needed iff the proposal id is different AND no occlusion OR implausible Situation occurs AND not oob
      if (!edgeCase1) seg2_0 = -1;
      if (!edgeCase2) seg2_1 = -1;
      // only if not occluded or implausible: and not oob: varId0x !=-2
      if ( seg2_0 >=0 || seg2_1 >=0 ) // edgeCase1==1 if checking is necessary
        motionPairsFW.push_back( MotionPairs( segId, varIdx, seg2_0, varId0x, seg2_1, varId1x ) );// also occlusion procedure!
    }
  }/////////////////
  ///////////////// BW: ////////////////////////
#ifndef _FW_ONLY_
  const std::vector< std::pair<int,int> > &bwMap = temp.seg2segBW;
  for(int i= 0; i < bwMap.size(); i++ )
  {
    int edgeCase1=0;
    int edgeCase2=0;

    // 1. assign a label 0 to segi in left view:
    int segId = bwMap[i].first;
    if (segId<0) continue;
    // so. check onto which segment the segment projects to and check if it is considered
    // how? here:
    int prop0   = currentSolution2[ segId ];

    int prop1   = trial;
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment
    int seg2_1  = bwMap[i].second;                           // 0 in view2 maps to this segment
    int varId0x = -2;//does not project to anything. 
    int varId1x = -2;
    int varIdx  = temp.seg2varBW[ segId ];// id of the first variable

    // NEW _ TEST
    if ( varIdx < 0) continue;

    if (seg2_0 > -1 )
      varId0x = temp.seg2varFW[ seg2_0 ];// variable id (can be -1 == not considered)
    if (seg2_1 > -1 )
      varId1x = temp.seg2varFW[ seg2_1 ];// should be >=0
    // these are the local variables. 

    if ( prop0 == trial && (varId0x<=-1 || currentSolution1[seg2_0] == trial) ) continue; // no impact 

    //      varId0x =-2;
    //      varId1x =-2;
    // first addres case 00 and 01:
    // then binary, else 00 is a unary, 01 does not exist, no matter what 
    if (varId0x == -2) // maps on no other segment:
      //        unaries[ varIdx ] += P2( currentSolutionE.dataBW[ segId ], 0);
      unaries.push_back( Unary<Scalar>(varIdx, P2( currentSolutionE.dataBW[ segId ], 0) ));
    else // segment exists but no variable for it:, so 01 does not exist
    {
      int prop00 = currentSolution1[seg2_0];
      Scalar penalty =0;
      if( prop00 == prop0 ) // normal data term
        penalty = currentSolutionE.dataBW[ segId ];
      else // must check depth in both versions:
      {
        penalty = getPenalty ( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ segId ], areas2[ segId ],
          currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
          edgeCase1, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
        penalty += addPotts * areas2[ segId ];
      }
      if ( varId0x == -1 ) // partner without variable index
        unaries.push_back( Unary<Scalar>(varIdx, P2( penalty, 0) ));
      //          unaries[ varIdx ] += P2( penalty , 0);
      else
      {
        Scalar pen00 = penalty;

        int prop01 = trial;
        Scalar penalty =0;
        if( prop01 == prop0 ) // normal data term
          penalty = currentSolutionE.dataBW[ segId ];
        else // must check depth in both versions
        {
          penalty = getPenalty ( ihoms_pix[prop0], noms_pix[prop0], temp.Nom, iCenters[ segId ], areas2[ segId ],
            currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
            edgeCase1, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
          penalty += addPotts * areas2[ segId ];
        }
        ////
        if ( edgeStoreMap[(varId0x - temp.firstIdFW)].first != varIdx ) 
        {
          if ( edgeStoreMapB[varId0x - temp.firstIdFW].first != varIdx )
            binaries.push_back( Binary<Scalar> (varId0x, varIdx, P4 (pen00, 0, penalty,0 ) ) );//new
          else
            binaries[ edgeStoreMapB[varId0x - temp.firstIdFW].second ].edge += P4 (pen00, 0, penalty,0 );//old
          //            binaries.push_back( Binary<Scalar> (varId0x, varIdx, P4 (pen00, 0, penalty,0 ) ) );
        }
        else
          binaries[ edgeStoreMap[varId0x - temp.firstIdFW].second ].edge += P4 (pen00, 0, penalty,0 );//old
        // same for 01 case, an edge exists:
      }
    }///////////////////////////////////////////
    // case 10, 11

    if (varId1x == -2) // maps on no other segment:
      unaries.push_back( Unary<Scalar>(varIdx, P2( 0, temp.dataBW[ i ] )) );
    else // segment exists but no variable for it:, so 01 does not exist
    {
      // only need to consider 10:
      // same segi? or dame depth? or larger or smaller -> 3 possibilities
      int prop10 = currentSolution1[seg2_1];
      Scalar penalty =0;
      if( prop10 == trial ) // normal data term
        penalty = temp.dataBW[ i ];
      else // must check depth in both versions:
      {
        penalty = getPenalty ( temp.iHom , temp.Nom, noms_pix[prop10], iCenters[ segId ], areas2[ segId ],
          temp.dataBW[ i ], temp.freedataBW[ i ], temp.freeBW[ i ], edgeCase2,
          autoScoresBW ? adaptiveAutoScoresBW[segId] : 1.);
        penalty += addPotts * areas2[ segId ];
      }
      if ( varId1x == -1 ) // not encoded variable
        unaries.push_back( Unary<Scalar>(varIdx, P2( 0, penalty )) );
      //          unaries[ varIdx ] += P2( 0, penalty );
      else  // is encoded
      {
        Scalar pen10 = penalty;

        int prop11 = trial;
        Scalar penalty =0;
        if( prop11 == prop1 ) // normal data term
          penalty = temp.dataBW[ i ];//currentSolutionE.dataBW[ segId ];
        else // must check depth in both versions: -- NONSENSE FROM HERE
        {
          assert( 0 );
          // get depth relationship:
          penalty = getPenalty ( temp.iHom, temp.Nom, noms_pix[prop11], iCenters[ segId ], areas2[ segId ],
            temp.dataBW[ i ], temp.freedataBW[ i ], temp.freeBW[ i ], edgeCase2, 
            autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
          penalty += addPotts * areas2[ segId ];
        }

        if ( edgeStoreMap[(varId1x - temp.firstIdFW)].first != varIdx ) 
        {
          if ( edgeStoreMapB[(varId1x - temp.firstIdFW)].first != varIdx ) 
          {
            if (varId0x != varId1x) // new:
              binaries.push_back( Binary<Scalar> (varId1x, varIdx, P4 (0, pen10, 0, penalty ) ) );
            else // double entry:
              binaries[binaries.size()-1].edge += P4 (0, pen10, 0, penalty ); // same position
          }
          else
            binaries[ edgeStoreMapB[(varId1x - temp.firstIdFW)].second ].edge += P4 (0, pen10, 0, penalty );
        }
        else
        {
          binaries[ edgeStoreMap[(varId1x - temp.firstIdFW)].second ].edge += P4 (0, pen10, 0, penalty );
        }
      }
    }///////////////////////////////////////////

    if (computeMotionPairs)
    {
      // check needed iff the proposal id is different AND no occlusion OR implausible Situation occurs AND not oob
      if (!edgeCase1) seg2_0 = -1;
      if (!edgeCase2) seg2_1 = -1;
      // only if not occluded or implausible: and not oob: varId0x !=-2
      if ( seg2_0 >=0 || seg2_1 >=0 ) // edgeCase1==1 if checking is necessary
        motionPairsBW.push_back( MotionPairs( segId, varIdx, seg2_0, varId0x, seg2_1, varId1x ) );// also occlusion procedure!
    }
  }/////////////////
#endif

  expandRegionLongRangeInteractions( propId );
}

template<class Scalar>
void
  Datasegments<Scalar>::
  updateSolutions( int propId, std::vector<int>& _currentSolution1, std::vector<int>& _currentSolution2, const std::vector<int>& solution01 )
{
  // dataElements: n normals/rotations. 
  // propId      : defines where (center) expansion takes place
  const dataElem<Scalar> &temp = dataElements[propId];
  int trial = temp.proposal;

  const std::vector< std::pair<int,int> > &fwMap = temp.seg2segFW;
  const std::vector< std::pair<int,int> > &bwMap = temp.seg2segBW;
  //
  // run along current solution, if something changed, update it: no 
  // run along last trial 
  int nn = 0;
  for (int i=0;i< fwMap.size(); i++,nn++)
  {
    // NEW _ TEST
    if (temp.seg2varFW[ fwMap[i].first ] < 0) continue;

    // a change:
    if ( solution01[ temp.seg2varFW[ fwMap[i].first ] ] )
      //      if ( solution01[nn] )
    {
      currentSolutionE.freeFW[ fwMap[i].first ] = temp.freeFW[i];
      currentSolutionE.dataFW[ fwMap[i].first ] = temp.dataFW[i];
      currentSolutionE.freedataFW[ fwMap[i].first ] = temp.freedataFW[i];
      currentSolutionE.seg2segFW[ fwMap[i].first ]  = temp.seg2segFW[i];

      currentSolution1[  fwMap[i].first ] = trial;
      _currentSolution1[ fwMap[i].first ] = trial;
    }
  }
#ifndef _FW_ONLY_
  for (int i=0;i < bwMap.size(); i++, nn++)
  {
    // NEW _ TEST
    if (temp.seg2varBW[ bwMap[i].first ] < 0) continue;

    // a change:
    if ( solution01[ temp.seg2varBW[ bwMap[i].first ] ] )
      //      if ( solution01[nn] )
    {
      currentSolutionE.freeBW[ bwMap[i].first ] = temp.freeBW[i];
      currentSolutionE.dataBW[ bwMap[i].first ] = temp.dataBW[i];
      currentSolutionE.freedataBW[ bwMap[i].first ] = temp.freedataBW[i];
      currentSolutionE.seg2segBW[ bwMap[i].first ]  = temp.seg2segBW[i];
      currentSolution2[  bwMap[i].first ] = trial;
      _currentSolution2[ bwMap[i].first ] = trial;
    }
  }
#endif
}


template<class Scalar>
void
  Datasegments<Scalar>::
  expandRegionLongRangeInteractions( int propId )
{
  const dataElem<Scalar> &temp = dataElements[propId];
  int trial = temp.proposal;
  // check all forwards, label 0, where they project to:
  for ( int i=0;i<currentSolution1.size();i++ )
  {
    int segId = i;
    // now if the segment is a variable : next
    if ( temp.seg2varFW[segId] > -1 ) continue;

    int prop0   = currentSolution1[ segId ];
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// the segi on which it projects to
    if (seg2_0 <0) continue;
    // now if the segment is not a variable : next
    int varId0x = temp.seg2varBW[ seg2_0 ];// variable id (can be -1 == not considered)
    if ( varId0x <= -1 ) continue;

    ///////////////////
    // so from currentsolution a segment not a variable 
    // falls onto a segment in 2nd view 
    // which is a variable:
    ///////////////////

    // otherwise: we need to generate a data term:
    int prop00  = currentSolution2[ seg2_0 ];
    int prop01  = trial;

    int edgeCase1(0), edgeCase2(0);

    // find both penalties
    Scalar pen00 = getPenalty ( homs_pix[prop0], inoms_pix[prop0], inoms_pix[prop00], Centers[ segId ], areas1[ segId ],
      currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
      edgeCase1, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. ) + (prop0 != prop00) * areas1[ segId ] * addPotts;
    Scalar pen01 = getPenalty ( homs_pix[prop0], inoms_pix[prop0], temp.iNom, Centers[ segId ], areas1[ segId ],
      currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], currentSolutionE.freeFW[ segId ], 
      edgeCase2, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. ) + (prop0 != trial ) * areas1[ segId ] * addPotts;

    if(pen00 !=pen01)
      unaries.push_back( Unary<Scalar>(varId0x, P2( pen00, pen01) )) ;

    if (computeMotionPairs)
    {
      // check needed iff the proposal id is different AND no occlusion OR implausible Situation occurs AND not oob
      int seg2_1 = seg2_0;//int varId1x = varId0x;
      if (!edgeCase1 || (prop0 == prop00)) {seg2_0 = -1;}
      if (!edgeCase2 || (prop0 == trial))  {seg2_1 = -1;}
      // only if not occluded or implausible: and not oob: varId0x !=-2
      if ( seg2_0 >=0 || seg2_1 >=0 ) // edgeCase1==1 if checking is necessary
        motionPairsFW.push_back( MotionPairs( segId, -1, seg2_0, varId0x, seg2_1, varId0x ) );// also occlusion procedure!
    }
  }
  // and backward
#ifndef _FW_ONLY_
  for ( int i=0;i<currentSolution2.size();i++ )
  {
    int segId = i;
    // now if the segment is a variable : next
    if ( temp.seg2varBW[segId] > -1 ) continue;

    int prop0   = currentSolution2[ segId ];
    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// the segi on which it projects to
    if (seg2_0 <0) continue;
    // now if the segment is not a variable : next
    int varId0x = temp.seg2varFW[ seg2_0 ];// variable id (can be -1 == not considered)
    if ( varId0x <= -1 ) continue;

    // otherwise: we need to generate a data term:
    int prop00  = currentSolution1[ seg2_0 ];
    int prop01  = trial;

    int edgeCase1(0), edgeCase2(0);
    // find both penalties
    Scalar pen00 = getPenalty ( ihoms_pix[prop0], noms_pix[prop0], noms_pix[prop00], iCenters[ segId ], areas2[ segId ],
      currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
      edgeCase1, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. ) + (prop0 != prop00) * areas2[ segId ] * addPotts;
    Scalar pen01 = getPenalty ( ihoms_pix[prop0], noms_pix[prop0], temp.Nom, iCenters[ segId ], areas2[ segId ],
      currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], currentSolutionE.freeBW[ segId ], 
      edgeCase2, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. ) + (prop0 != trial ) * areas2[ segId ] * addPotts;

    if(pen00 !=pen01)
      unaries.push_back( Unary<Scalar>(varId0x, P2( pen00, pen01) )) ;

    if (computeMotionPairs)
    {
      // check needed iff the proposal id is different AND no occlusion OR implausible Situation occurs AND not oob
      int seg2_1 = seg2_0;
      if (!edgeCase1 || (prop0 == prop00)) seg2_0 = -1;
      if (!edgeCase2 || (prop0 == trial))  seg2_1 = -1;
      // only if not occluded or implausible: and not oob: varId0x !=-2
      if ( seg2_0 >=0 || seg2_1 >=0 ) // edgeCase1==1 if checking is necessary
        motionPairsBW.push_back( MotionPairs( segId, -1, seg2_0, varId0x, seg2_1, varId0x ) );// also occlusion procedure!
    }
  }
#endif
}


/// does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled: p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
template<class Scalar>
Scalar
  Datasegments<Scalar>::
  // area: full area of segment, freePixel how many of these are actually 'free'
  getPenalty ( const M3& Hom, P3 Nom, P3 iNom, P3 centre, int area, Scalar data, Scalar freeData, Scalar freePix, int& edgeCase, Scalar autoScore)
{
  edgeCase = 0;
  if (simpleApproach)
  {
    edgeCase = 1;
    return data;
  }

  Scalar localOccPen = xtraPen+autoScore;// ?? 0.2 -> TODO ALSO FOR OOB

#ifdef _simple_data_
  return data; // what happens? : not working as well, but no energy increase
#endif

  P3 p = Hom * centre;p /= p[2];//projected point
  Scalar depth0      = 1./(p|Nom); 
  Scalar depthThresh = max ( maxDepthThresh, fabs(depth0  * __OCC_THRESH__) );
  Scalar depthDiff   = depth0 - 1./(p|iNom);// let depths by positive -> depthDiff>0 implies an occlusion case, sign change: < 0 occlusion case

  // necessary for the still a prop variant, i.e. the init is like this - should never appear else  
  // NEVER since the 2nd variable is == -1 (set in initNewSeg/initSegFromView)
  if ( depth0 > Scalar(-_minDepth_) )
#ifdef _inFrontPenalty_
    return max(data, data-freeData + area * Scalar(impFactor)*localOccPen) - area*addPotts; // special penalty for areas behind cam ?
#else
    return data-freeData + freePix*localOccPen - area*addPotts;// not here
#endif

  if ( depthDiff > depthThresh )// projection occludes other point but the other is visible: impossible
    if ( impossibleApproach )
      return max(data, data-freeData + area * Scalar(impFactor)*localOccPen) - area*addPotts;
    else
      return max(data, data-freeData + area * localOccPen) - area*addPotts;
  else // occlusion case:
    if ( depthDiff < -depthThresh ) // projection is occluded by other point: occlusion model
      return data-freeData + freePix*localOccPen - area*addPotts; // kills gross errors: NICE! EPE
    else
    {
      edgeCase = 1;// ok case
      return data;
    }
}
#endif