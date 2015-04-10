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

#ifndef __DataPerPixel__cpp
#define __DataPerPixel__cpp

#include "DataDefinitionsVC.h"
#include "DataSegments.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;

/////////////////////////////////////////////////////////

template<class Scalar>
void
  Datasegments<Scalar>::
  set_AutoScoresPix(Scalar* _autoScores1, Scalar* _autoScores2)
{
  adaptiveAutoScoresFW.clear();    adaptiveAutoScoresBW.clear();
  autoScoresFW = _autoScores1;autoScoresBW = _autoScores2;
  gs1.setAutoScores( autoScoresFW );
  gs2.setAutoScores( autoScoresBW );
  updateAutoScoresPix();
}

#ifndef _checkValiadityAuto_
template<class Scalar>
void
  Datasegments<Scalar>::
  updateAutoScoresPix( )
{
  if (autoScoresFW != NULL)
  {
    if ( adaptiveAutoScoresFW.size() == 0 )
      adaptiveAutoScoresFW.insert(adaptiveAutoScoresFW.end(), &autoScoresFW[0], &autoScoresFW[w*h]);
    else
      for (int i=0;i< adaptiveAutoScoresFW.size(); i++)
      {
        if (currentSolutionE.freeFW[i] > 0) 
          adaptiveAutoScoresFW[i] = std::min( adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i] );
      }
  }


  if (autoScoresBW != NULL)
  {
    if ( adaptiveAutoScoresBW.size() == 0 )
      adaptiveAutoScoresBW.insert(adaptiveAutoScoresBW.end(), &autoScoresBW[0], &autoScoresBW[w*h]);
    else
      for (int i=0;i< adaptiveAutoScoresBW.size(); i++)
      {
        if (currentSolutionE.freeBW[i] > 0) 
          adaptiveAutoScoresBW[i] = std::min( adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i] );
      }
  }
}
#else
template<class Scalar>
void
  Datasegments<Scalar>::
  updateAutoScoresPix( )
{
  if (autoScoresFW != NULL)
  {
    if ( adaptiveAutoScoresFW.size() == 0 )
      adaptiveAutoScoresFW.insert(adaptiveAutoScoresFW.end(), &autoScoresFW[0], &autoScoresFW[w*h]);
    else
      for (int i=0;i< adaptiveAutoScoresFW.size(); i++)
      {
        if (currentSolutionE.freeFW[i] > 0) 
        {
          int prop0   = segImg[ i ];
          int pix2_0  = currentSolutionE.seg2segFW[ i ].second;// projects onto
          
          if (pix2_0 >=0) // maps on other pixel:
          {
            int prop00 = segImg2[ pix2_0 ];
            if( prop00 != prop0 && whichPenalty ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop00].iNom, P3FromGlobal(i), dataElements[prop00].iHom, P3FromGlobal(pix2_0), prop00, prop0 ) != 3 )
             continue;
          }
          adaptiveAutoScoresFW[i] = std::min( adaptiveAutoScoresFW[i], currentSolutionE.freedataFW[i] );
        }
      }
  }
  if (autoScoresBW != NULL)
  {
    if ( adaptiveAutoScoresBW.size() == 0 )
      adaptiveAutoScoresBW.insert(adaptiveAutoScoresBW.end(), &autoScoresBW[0], &autoScoresBW[w*h]);
    else
      for (int i=0;i< adaptiveAutoScoresBW.size(); i++)
      {
        if (currentSolutionE.freeBW[i] > 0) 
        {
          int prop0   = segImg2[ i ];
          int pix2_0  = currentSolutionE.seg2segBW[ i ].second;// projects onto

          if (pix2_0 >= 0 )
          {
            int prop00 = segImg[ pix2_0 ];
            if( prop00 != prop0 && whichPenalty ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop00].Nom, P3FromGlobal(i), dataElements[prop00].Hom, P3FromGlobal(pix2_0), prop00, prop0 ) != 3 )
              continue;
          }
          adaptiveAutoScoresBW[i] = std::min( adaptiveAutoScoresBW[i], currentSolutionE.freedataBW[i] );
        }
      }
  }
}
#endif

// to be done for the  pairs with reference outgoing images
/// why not precomputing all data terms again in advance - because storage == 50*50==2500 * size image * 8 images ,etc. - no way
/// so simply precompute necessary stuff only no data until later
/// seg2seg is now pixel to pixel and seg2var is pixeltovar
/// needs also map from local to global pixe coords! (for smoothing)
template<class Scalar>
void
  Datasegments<Scalar>::
  preInitPixel( int proposal, genHomoG<Scalar>* gHom, const std::vector<P3>* normals ) // int boxX, int boxY, int expansionArea, 
{
  ////// first project the bbox to get a box in the second image:
  M3 Hom  = gHom->getHomC0( proposal );
  M3 iHom = Hom;iHom.invert();

  P3 vn_0  = gHom->getViewNormalC0( proposal ); // no iiKlt2 * 

  newElem.clear();
  newElem.varsFW =0;
  newElem.varsBW =0;

  // in cam coordinates, so not vnom * pixel, but vNom * Kl^-1*pixel
  if (vNoms_cam.size() <= proposal)
  {
    vNoms_cam.push_back(-(vn_0));// camera coords
    homs_pix.push_back ( Hom);
    ihoms_pix.push_back (iHom);
  }
  else
  {
    vNoms_cam[proposal] = -(vn_0); // operate on ?
    homs_pix[proposal]  = Hom;
    ihoms_pix[proposal]  = iHom;
  }

  newElem.proposal = proposal;
  newElem.Hom  = Hom;
  newElem.Nom  = iKlt * (*normals)[proposal];
  newElem.iHom = iHom;
  newElem.iNom = iKlt2 * vn_0; // operate an something different then above?
  newElem.w    = w;
  newElem.h    = h;

  dataElements.push_back( newElem ); // stored by proposal
}


template<class Scalar>
void
  Datasegments<Scalar>::
  preinitPixelFromView( int proposal, const dataElem<Scalar>& leader, const dataElem<Scalar>& second)
{
  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(nSegments);

  // map all elements form the leader to other view copying viewnormals, 

  ////// first project the bbox to get a box into the second image!:
  // this is pushing back to ref image then to RT image !
  M3 Hom  = second.Hom * leader.iHom;
  M3 iHom = Hom;iHom.invert();

  P3 normal  = leader.iNom;
  P3 vn_0    = second.iNom;

  /////////////
  // all those within the bounds are interesting for data term and variables:
  newElem.clear();
  newElem.varsFW =0;
  newElem.varsBW =0;

  newElem.proposal = proposal;
  newElem.Hom  = Hom;
  newElem.Nom  = normal;
  newElem.iHom = iHom;
  newElem.iNom = vn_0;
  newElem.w    = w;
  newElem.h    = h;

  // in cam coordinates, so not vnom * pixel, but vNom * Kl^-1*pixel
  if (vNoms_cam.size() <= proposal)
  {
    vNoms_cam.push_back(-(iiKlt2*vn_0));// which is this time the same as in newElem.Nom/iNom, while different above
    homs_pix.push_back(Hom);
    ihoms_pix.push_back(iHom);
  }
  else
  {
    vNoms_cam[proposal]  = -(iiKlt2*vn_0); // to turn from pixel to camera coords
    homs_pix[proposal]   = Hom;
    ihoms_pix[proposal]  = iHom;
  }
  dataElements.push_back( newElem ); // stored by proposal so thats not the caes TODO
}

/// which is the 'best' order for the expansion moves ? does not matter ?
template<class Scalar>
void
  Datasegments<Scalar>::
  sortInitNewPixel( std::vector<int>& order, std::vector<int>& seg2Prop )
{
  std:: vector <std::pair< Scalar, int >  > sortMe;
  for ( int i =0; i< seg2Prop.size(); i++ )
  {
    const dataElem<Scalar>& elem  = dataElements[ seg2Prop[i] ];
    // 
    Scalar dd = -1./(Centers[i] | elem.Nom);
    sortMe.push_back ( std::pair<Scalar, int> (dd, i) );
  }
  std::sort( sortMe.begin(), sortMe.end() );

  order.assign(seg2Prop.size(), 0);
  for ( int i=0;i<order.size(); i++ )
    //      order[i] = sortMe[i].second; 
    order[order.size()-1-i] = sortMe[i].second;
}

template<class Scalar>
bool
  Datasegments<Scalar>::
  testNewPixel( int centerId, int proposal, int boxX, int boxY, int expansionArea )
{
  int cSize = 3; 

  assert(proposal < dataElements.size());

  dataElem<Scalar> pot_NewElem = dataElements[proposal];
  P3 p0  = Kl * P3( &centers[3*centerId] );// center is in camera coords
  p0    /= p0[2];

  // this time really per pixel:  // do i really want to restrict this or handle it later ?
  int startX = int(floor ( p0[0] -0.5 )) - boxX;// changed restrict to boundary here
  int endX   = int(floor ( p0[0] -0.5 )) + boxX;
  int startY = int(floor ( p0[1] -0.5 )) - boxY;
  int endY   = int(floor ( p0[1] -0.5 )) + boxY;
  endX = max(endX , startX );startX = min(endX , startX );
  endY = max(endY , startY );startY = min(endY , startY );
  box1 = P4i (startX, startY, endX, endY);

  ////// first project the bbox to get a box in the second image:
  M3  Hom  = pot_NewElem.Hom;
  // like this since consider right border and right camera.
  // there are (a lot of) pixels where the partner is not in the image borders
  // of the left cam. So the uncut! box must be projected
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

  minX = min(minX, p0[0]);minX = min(minX, p1[0]);
  minY = min(minY, p0[1]);minY = min(minY, p1[1]);
  maxX = max(maxX, p0[0]);maxX = min(maxX, p1[0]);
  maxY = max(maxY, p0[1]);maxY = min(maxY, p1[1]);
  maxX = max(maxX , minX );minX = min(maxX , minX );
  maxY = max(maxY , minY );minY = min(maxY , minY );

  startX = max(startX-cSize, 0);
  startY = max(startY-cSize, 0);
  endX   = min(endX+cSize, w-1);
  endY   = min(endY+cSize, h-1);
  endX   = max(endX , startX );startX = min(endX , startX );
  endY   = max(endY , startY );startY = min(endY , startY );
  box1   = P4i (startX, startY, endX, endY);

  // make boxes 1 pixel larger for uniform handling
  box1 += P4i(0,0,1,1);

  /////////////
  // TODO: pixel totally oob are still treated - that is stupid
  // new:
  box1.minimize( P4i (w,h,w,h) );
  box1.maximize( P4i (0,0,0,0) );

  // 1. compare with elements under revision:
  for (int i=0;i<storedElem.size();i++)
    if (storedElem[i].proposal != -1) 
      if ( (box1[0] < storedElem[i].boxFW[2] && box1[2] > storedElem[i].boxFW[0]) && 
        (box1[1] < storedElem[i].boxFW[3] && box1[3] > storedElem[i].boxFW[1]) )
        return false;

  // seg2seg is now pixel to pixel and seg2var is pixeltovar
  int innerId=0;
  for (int i=box1[0]; i < box1[2]; i++)
    for (int j=box1[1]; j < box1[3]; j++, innerId++)
    {
      P3 p0 = Hom * P3( i+1, j+1,1. );p0 /= p0[2];

      if ( isinBounds(p0) )
      {
        minX = min(p0[0], minX);
        maxX = max(p0[0], maxX);
        minY = min(p0[1], minY);
        maxY = max(p0[1], maxY);
      }
    }

    minX = max(int(floor(minX-0.5))-cSize, 0);
    minY = max(int(floor(minY-0.5))-cSize, 0);
    maxX = min(int(floor(maxX-0.5))+cSize, w-1);
    maxY = min(int(floor(maxY-0.5))+cSize, h-1);
    maxX = max(maxX , minX );minX = min(maxX , minX );
    maxY = max(maxY , minY );minY = min(maxY , minY );
    box2 = P4i (minX, minY, maxX, maxY);

    // 2. compare with elements under revision (this time 2nd box):
    //    newElem.boxBW;
    // return if no overlap
    for (int i=0;i<storedElem.size();i++)
      if (storedElem[i].proposal != -1) 
        if ( (box2[0] < storedElem[i].boxBW[2] && box2[2] > storedElem[i].boxBW[0]) && 
          (box2[1] < storedElem[i].boxBW[3] && box2[3] > storedElem[i].boxBW[1]) )
          return false;

    return true;
}


/// seg2seg is now pixel to pixel and seg2var is pixeltovar segId init the position of the expansion, nothing else
template<class Scalar>
int
  Datasegments<Scalar>::
  initNewPixel( int centerId, int proposal, int boxX, int boxY, int expansionArea, int addToSeg2Vars, int addTo2AllVars )//, genHomoG<Scalar>* gHom, const std::vector<P3>* normals, 
{
  int cSize = 3; 
  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(nSegments);

  newElem = dataElements[proposal];
  newElem.proposal = proposal;

  P3 p0  = Kl * P3( &centers[3*centerId] );// center is in camera coords
  p0    /= p0[2];

#ifdef _DEBUG
  int segmentId0 = -1;
  int inBounds = isinBounds(p0);
  if ( inBounds )
  {
    int pixId2 = toPixId( p0 );
    segmentId0 = segImg[ pixId2 ];
  }
#endif

  // this time really per pixel:  // do i really want to restrict this or handle it later ?
  int startX = int(floor ( p0[0] -0.5 )) - boxX;// changed restrict to boundary here
  int endX   = int(floor ( p0[0] -0.5 )) + boxX;
  int startY = int(floor ( p0[1] -0.5 )) - boxY;
  int endY   = int(floor ( p0[1] -0.5 )) + boxY;
  endX = max(endX , startX );startX = min(endX , startX );
  endY = max(endY , startY );startY = min(endY , startY );
  box1 = P4i (startX, startY, endX, endY);

  ////// first project the bbox to get a box in the second image:
  M3  Hom  = newElem.Hom;
  M3 iHom  = newElem.iHom;
  /// act on camera coords, others on image coords 
  M3 HomCC  = Hom * Kl;
  M3 iHomCC = iHom * Kl2;

  // like this since consider right border and right camera.
  // there are (a lot of) pixels where the partner is not in the image borders
  // of the left cam. So the uncut! box must be projected
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

  // there are cases when all is moved outside, so no back data term; but the motion is gigantic

  minX = min(minX, p0[0]);minX = min(minX, p1[0]);//minX = max(int(floor(minX-0.5)), 0);
  minY = min(minY, p0[1]);minY = min(minY, p1[1]);//minY = max(int(floor(minY-0.5)), 0);
  maxX = max(maxX, p0[0]);maxX = min(maxX, p1[0]);//maxX = min(int(floor(maxX-0.5)), w-1);
  maxY = max(maxY, p0[1]);maxY = min(maxY, p1[1]);//maxY = min(int(floor(maxY-0.5)), h-1);
  maxX = max(maxX , minX );minX = min(maxX , minX );
  maxY = max(maxY , minY );minY = min(maxY , minY );
  //    box2 = P4i (minX, minY, maxX, maxY);

  startX = max(startX-cSize, 0);
  startY = max(startY-cSize, 0);
  endX   = min(endX+cSize, w-1);
  endY   = min(endY+cSize, h-1);
  endX   = max(endX , startX );startX = min(endX , startX );
  endY   = max(endY , startY );startY = min(endY , startY );
  box1   = P4i (startX, startY, endX, endY);

  // make boxes 1 pixel larger for uniforma handling
  box1 += P4i(0,0,1,1);

  /////////////
  // TODO: pixel totally oob are still treated - that is stupid
  // new:
  box1.minimize( P4i (w,h,w,h) );
  box1.maximize( P4i (0,0,0,0) );

  newElem.seg2varFW.resize( (box1[2]-box1[0])*(box1[3]-box1[1]), -1);// pix 2 variable here
  newElem.varsFW =0;
  newElem.varsBW =0;

  pix_idxFW.clear();pix_idxFW.resize( newElem.seg2varFW.size() );
  pix_idyFW.clear();pix_idyFW.resize( newElem.seg2varFW.size() );

  // seg2seg is now pixel to pixel and seg2var is pixeltovar
  int innerId=0;
  for (int i=box1[0]; i < box1[2]; i++)
    for (int j=box1[1]; j < box1[3]; j++, innerId++)
    {
      // int innerId = i*(box1[3]-box1[1]) +j;
      // if segImg in large coordiantes == proposal: skip

      int pixId = int( i )*h +int ( j );
      int pixId2 = -1;
      // lookup partner:
      int partner = -1; // -1: no partner
      P3 p0 = Hom * P3( i+1, j+1,1. );p0 /= p0[2];
      pix_idxFW[innerId] = p0[0];pix_idyFW[innerId] = p0[1];

      int inBounds = isinBounds(p0);
      if ( inBounds )
      {
        pixId2 = toPixId(p0);
        partner =  segImg2[pixId2]; // thats the solution - i need the pixel
        // a pixel is a pixel large, so if expansion area is > 0.5 no problemo
        minX = min(p0[0], minX);
        maxX = max(p0[0], maxX);
        minY = min(p0[1], minY);
        maxY = max(p0[1], maxY);
      }

      // not interesting, though box should include the pixel, since in second image - or other image ..
      segmentsareas1.push_back( std::pair<int,int> (pixId, pixId2) );

      // NEW: TODO check if sane - should deliver speed up as well : collapses though
#ifdef _Shortcut_Speedup_
      if ( (i-box1[0] >= cSize || i < cSize ) && (box1[2]-i > cSize  || i >= w-cSize ) && (j-box1[1] >= cSize || j < cSize ) && (box1[3]-j > cSize || j >= h-cSize ) )
      {
        if (segImg[ pixId ] == proposal) 
          continue;
        //          if ( 1./(newElem.Nom|P3( i+1, j+1,1. )) > -1. ) // too close
        //            continue;
      }
#endif
      // ---------------------------------------

      if ( (i-box1[0] >= cSize || i < cSize ) && (box1[2]-i > cSize  || i >= w-cSize ) && (j-box1[1] >= cSize || j < cSize ) && (box1[3]-j > cSize || j >= h-cSize ) )
        newElem.seg2varFW[innerId] = addTo2AllVars+newElem.varsFW++; // all in box are variables
    }

    minX = max(int(floor(minX-0.5))-cSize, 0);
    minY = max(int(floor(minY-0.5))-cSize, 0);
    maxX = min(int(floor(maxX-0.5))+cSize, w-1);
    maxY = min(int(floor(maxY-0.5))+cSize, h-1);
    maxX = max(maxX , minX );minX = min(maxX , minX );
    maxY = max(maxY , minY );minY = min(maxY , minY );
    box2 = P4i (minX, minY, maxX, maxY);

    // restrict maximum expansion size to ..
    if ((maxX-minX)*(maxY-minY) > 0 && (maxX-minX)*(maxY-minY) > (boxX*2+1)*(boxY*2+1)*resizer*resizer )
    {
      /// center maps to :
      p0  = HomCC * P3( &centers[3*centerId] );p0/=p0[2];
      p0[0] = int(floor(p0[0]-0.5));p0[1] = int(floor(p0[1]-0.5));

#ifdef _writeDebugOut_
      printf("Fixed Box: before %.1f  (minX,maxX,..)  %.1f,%.1f,%.1f,%.1f ", (maxX-minX)*(maxY-minY), minX, maxX, minY, maxY );
#endif

      minX = max (minX, (p0[0]) - resizer*boxX);
      maxX = min (maxX, (p0[0]) + resizer*boxX);
      minY = max (minY, (p0[1]) - resizer*boxY);
      maxY = min (maxY, (p0[1]) + resizer*boxY);
      maxX = max(maxX , minX );minX = min(maxX , minX );
      maxY = max(maxY , minY );minY = min(maxY , minY );

#ifdef _writeDebugOut_
      printf("After: %.1f. center: [%.1f %.1f] (minX,maxX) %.1f,%.1f,%.1f,%.1f\n", (maxX-minX)*(maxY-minY), p0[0], p0[1], minX, maxX, minY, maxY);        //mexEvalString("drawnow");
#endif
      //       printf("After: %.1f. \n", (maxX-minX)*(maxY-minY));        mexEvalString("drawnow");
      box2 = P4i (minX, minY, maxX, maxY);
    }

    // technical reason ?!
    // make boxes 1 pixel larger for uniform handling
    box2 += P4i(0,0,1,1);

    // new:
    box2.minimize( P4i (w,h,w,h) );
    box2.maximize( P4i (0,0,0,0) );

    newElem.seg2varBW.resize( (box2[2]-box2[0])*(box2[3]-box2[1]), -1);
    pix_idxBW.clear();pix_idxBW.resize( newElem.seg2varBW.size() );
    pix_idyBW.clear();pix_idyBW.resize( newElem.seg2varBW.size() );
    newElem.boxFW = box1;

    newElem.firstIdFW=0;
    newElem.firstIdBW = newElem.varsFW + addToSeg2Vars;

    innerId=0;
    for (int i=box2[0]; i < box2[2]; i++ )
      for (int j=box2[1]; j < box2[3]; j++, innerId++)
      {
        //        newElem.seg2varBW[innerId] = -1; // see above
        int pixId  = i*h + j;
        int oob = 0;

        if (i<0 || i>=w || j < 0 || j>=h) oob=1;//pixId = -1;// could also continue right here !

        int pixId2 = -1;
        // lookup partner:
        int partner = -1; // -1: no partner
        P3 p0 = iHom * P3( i+1, j+1,1. );p0 /= p0[2];
        pix_idxBW[innerId] = p0[0];pix_idyBW[innerId] = p0[1];

        int inBounds = isinBounds(p0);
        if ( inBounds )
        {
          pixId2 = toPixId( p0 );
          partner =  segImg[pixId2];
        }
        // TODO ON/OFF ??? now all in box are variables
        // not in box: could also do self check if not projecting onto itself - skip - but a variables should be produced though not neceessarily a binary
        //        if (partner != -1 && newElem.seg2varFW[newElem.getBoxIdFW(partner)] <0 ) continue;// could lead to strong shrinkage! e.g. oob stuff

        // oob can be possible?! blocking not so bad though, keeps current solution if good: fine
#ifdef _FW_ONLY_
        continue;
#endif
#ifdef _breakCircles
        // now only for partner search - therefore no binary, just unary - micro less ndef
        int bid = newElem.getBoxIdFW(pixId2);
        if (bid != -1 && segmentsareas1[bid].second != pixId )
          pixId2 = -1;
#endif

        segmentsareas2.push_back( std::pair<int,int> (pixId, pixId2) );// global2global

        // NEW: TODO check if sane - should deliver speed up as well
#ifdef _Shortcut_Speedup_
        if ( !oob && (i-box2[0] >= cSize || i < cSize ) && (box2[2]-i > cSize  || i >= w-cSize ) && (j-box2[1] >= cSize || j < cSize ) && (box2[3]-j > cSize || j >= h-cSize ) )
        {
          if (segImg2[ pixId ] == proposal) continue;

          //        if ( 1./(newElem.iNom|P3( i+1, j+1,1. )) > -1. ) // too close
          //          continue;
        }
#endif
        // -------------------------------

        // new: variables are only those inside box w.r.t. census size (need pixel ids for this and lookup positions as well)
        if ( !oob && (i-box2[0] >= cSize || i < cSize ) && (box2[2]-i > cSize  || i >= w-cSize ) && (j-box2[1] >= cSize || j < cSize ) && (box2[3]-j > cSize || j >= h-cSize ) )
          newElem.seg2varBW[innerId] = newElem.varsBW++  +  newElem.varsFW + addToSeg2Vars + addTo2AllVars;
      }

      newElem.seg2segFW = segmentsareas1;
      newElem.seg2segBW = segmentsareas2;
      newElem.boxFW = box1;
      newElem.boxBW = box2;
      // new last variable number
      return newElem.varsBW + addToSeg2Vars;
}

/*! normals here are stored vNormals in other view!!! these must exist - could also use a dummy version? providing vnormals nad homos already.
copies the seg2var assignments from BOTH inputs !!, centerId delivers the position! where to expand the box
*/
template<class Scalar>
void
  Datasegments<Scalar>::
  initNewPixelFromView( int centerId, int proposal, const dataElem<Scalar>& leader, const dataElem<Scalar>& second)//, int expansionArea, const P4i& _box1, const P4i& _box2, 
{
  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(nSegments);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(nSegments);
  assert(proposal < dataElements.size());
  //    assert(centerId < nSegments);
  assert(centerId < Centers.size());
  newElem = dataElements[proposal];

  newElem.firstIdFW = leader.firstIdBW;
  newElem.firstIdBW = second.firstIdBW;

  // map all elements form the leader to other view copying viewnormals, 

  ////// first project the bbox to get a box into the second image!:
  // this is pushing back to ref image then to RT image !
  M3 Hom  = newElem.Hom;//second.Hom * leader.iHom;
  M3 iHom = newElem.iHom;//Hom;iHom.invert();

  M3 HomCC  = Hom * Kl;
  M3 iHomCC = iHom * Kl2;

  P3 normal  = newElem.Nom;//leader.iNom;
  P3 vn_0    = newElem.iNom;//second.iNom;

  box1 = leader.boxBW;// ! still needed
  box2 = second.boxBW;// !

  // pointless to check data if oob
  //box1[0] = min(w, max(box1[0], 0));
  //box1[2] = min(w, max(box1[2], 0));
  //box1[1] = min(h, max(box1[1], 0));
  //box1[3] = min(h, max(box1[3], 0));

  /////////////
  // all those within the bounds are interesting for data term and variables:
  newElem.clear();
  newElem.seg2varFW.resize(leader.seg2varBW.size() , -1);
  newElem.seg2varBW.resize(second.seg2varBW.size(),  -1);
  newElem.varsFW =0;
  newElem.varsBW =0;
  newElem.boxFW = box1;
  pix_idxFW.clear();pix_idxFW.resize( newElem.seg2varFW.size() );
  pix_idyFW.clear();pix_idyFW.resize( newElem.seg2varFW.size() );
  pix_idxBW.clear();pix_idxBW.resize( newElem.seg2varBW.size() );
  pix_idyBW.clear();pix_idyBW.resize( newElem.seg2varBW.size() );

  int innerId = 0;
  for( int i =0; i< leader.seg2segBW.size(); i++, innerId++ )
  {
    int gPixId  = leader.seg2segBW[i].first;
    int gPixId2 = -1; 

    newElem.seg2varFW[i] = leader.seg2varBW[i];
    if (leader.seg2varBW[i]>=0) 
      newElem.varsFW++;// just copy

    P3 q = leader.getImageCoordBW( i );
    int inBounds1 = isinBounds(q);
#ifdef _DEBUG
    if ( !inBounds1 && leader.seg2varBW[i]>=0)
    {
      int weird  = leader.seg2segBW[i].first;
      int weird2 = leader.seg2segBW[i].second;
    }
#endif
    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = Hom * q;p0 /= p0[2];
    pix_idxFW[innerId] = p0[0];pix_idyFW[innerId] = p0[1];

    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      gPixId2 = toPixId( p0 );
      partner =  segImg2[gPixId2];
    }
    segmentsareas1.push_back( std::pair<int,int> (gPixId, gPixId2) );
  }
  assert(newElem.varsFW == leader.varsBW);

  innerId = 0;
  // why not project back all in the box2 and check, later consolidate anyway ?!
  for( int i =0; i< second.seg2segBW.size(); i++, innerId++ )
  {
    int gPixId  = second.seg2segBW[i].first;
    int gPixId2 = -1; 

    newElem.seg2varBW[i] = second.seg2varBW[i];
    if (second.seg2varBW[i]>=0) 
      newElem.varsBW++;

    P3 q = second.getImageCoordBW( i );
    int inBounds1 = isinBounds(q);
#ifdef _DEBUG
    if ( !inBounds1 && second.seg2varBW[i]>=0)
    {
      int weird  = leader.seg2segBW[i].first;
      int weird2 = leader.seg2segBW[i].second;
    }
#endif
    // lookup partner:
    int partner = -1; // -1: no partner
    P3 p0 = iHom * second.getImageCoordBW( i );p0 /= p0[2];
    pix_idxBW[innerId] = p0[0];pix_idyBW[innerId] = p0[1];

    int inBounds = isinBounds(p0);
    if ( inBounds )
    {
      gPixId2 = toPixId( p0 );
      partner =  segImg[gPixId2];
    }
    segmentsareas2.push_back( std::pair<int,int> (gPixId, gPixId2) );
  }
  assert( newElem.varsBW == second.varsBW );

  newElem.boxBW = box2;
  newElem.boxFW = box1;
  newElem.seg2segFW = segmentsareas1;
  newElem.seg2segBW = segmentsareas2;
  newElem.proposal = proposal;
  newElem.w    = w;
  newElem.h    = h;
}


template<class Scalar>
void
  Datasegments<Scalar>::
  createWarpsPixel( std::vector< Scalar >& idxFW, std::vector< Scalar >& idyFW, std::vector< Scalar >& idxBW, std::vector< Scalar >& idyBW )
{
  gws1.warp_noOmp_patchBased( idxBW, idyBW, newElem.seg2segFW );// actually should be increased by # of census border-pixel: done
  gws2.warp_noOmp_patchBased( idxFW, idyFW, newElem.seg2segBW );

  // per pixel stuff: ids maps to global ids. lIds: local ids, just -1: ignore, consider if >=0 -- saves into locScores: value per pixel in box!, global ids are now: seg2var.first
  if (idxFW.size() > 0)
    gs1.computeFullScoresMappedCensusNew_short( gws1.getOrigShortImage(), gws2.getWarpedShortImage(), &(idxFW[0]), &(idyFW[0]), 
    newElem.seg2segFW, newElem.seg2varFW, newElem.boxFW[3]-newElem.boxFW[1], newElem.boxFW[2]-newElem.boxFW[0], occlusions1 );
  else
    gs1.clearScores();

  if (idxBW.size() > 0)
    gs2.computeFullScoresMappedCensusNew_short( gws2.getOrigShortImage(), gws1.getWarpedShortImage(), &(idxBW[0]), &(idyBW[0]), 
    newElem.seg2segBW, newElem.seg2varBW, newElem.boxBW[3]-newElem.boxBW[1], newElem.boxBW[2]-newElem.boxBW[0], occlusions2 );
  else
    gs2.clearScores();

  newElem.dataFW = gs1.getScores() + gs1.getPenScores();
  newElem.dataBW = gs2.getScores() + gs2.getPenScores();
  newElem.freedataFW = gs1.getScores();
  newElem.freedataBW = gs2.getScores();
  newElem.freeFW.assign( newElem.dataFW.size() , 1);
  newElem.freeBW.assign( newElem.dataBW.size() , 1);
  // check box corners -- if violation check all newElem.seg2VarFW is a vector of pixel ids P3FromGlobal(ids) newElem.proposal  dataElements[xx].iNom
#ifdef __depthControl__

  P3 nom  = dataElements[ newElem.proposal ].Nom;
  P3 inom = dataElements[ newElem.proposal ].iNom;
  Scalar d0 ( 1./(nom|P3( newElem.boxFW[0]+1, newElem.boxFW[1]+1,1. )));
  Scalar d1 ( 1./(nom|P3( newElem.boxFW[2]+1, newElem.boxFW[1]+1,1. )));
  Scalar d2 ( 1./(nom|P3( newElem.boxFW[0]+1, newElem.boxFW[3]+1,1. )));
  Scalar d3 ( 1./(nom|P3( newElem.boxFW[2]+1, newElem.boxFW[3]+1,1. )));

 if ( d0 > Scalar(-_minDepth_) || d1 > Scalar(-_minDepth_) || d2 > Scalar(-_minDepth_) || d3 > Scalar(-_minDepth_) )
   for (int j=0;j<newElem.seg2segFW.size();j++)
   {
     int pid = newElem.seg2segFW[j].first;
     Scalar dd = Scalar(1.)/(nom|P3FromGlobal(pid));
     if ( dd > Scalar(-_minDepth_) )
     {
       newElem.dataFW [j] = Scalar(__depthControlAmplifier__ * impFactor);
       newElem.freedataFW [j] = 0;
       newElem.freeFW[j] = 0;
     }
   }

  d0 =( 1./(inom|P3( newElem.boxBW[0]+1, newElem.boxBW[1]+1,1. )));
  d1 =( 1./(inom|P3( newElem.boxBW[2]+1, newElem.boxBW[1]+1,1. )));
  d2 =( 1./(inom|P3( newElem.boxBW[0]+1, newElem.boxBW[3]+1,1. )));
  d3 =( 1./(inom|P3( newElem.boxBW[2]+1, newElem.boxBW[3]+1,1. )));

 if ( d0 > Scalar(-_minDepth_) || d1 > Scalar(-_minDepth_) || d2 > Scalar(-_minDepth_) || d3 > Scalar(-_minDepth_) )
   for (int j=0;j<newElem.seg2segBW.size();j++)
   {
     int pid   = newElem.seg2segBW[j].first;
     Scalar dd = Scalar(1.)/(inom|P3FromGlobal(pid));
     if ( dd > Scalar(-_minDepth_) )
     {
       newElem.dataBW [j] = Scalar(__depthControlAmplifier__ * impFactor);
       newElem.freedataBW [j] = 0;
       newElem.freeBW[j] = 0;
     }
   }
#endif
}

template<class Scalar>
void
  Datasegments<Scalar>::
  getEnergyMapPixel( std::vector<int> &energyMap, std::vector<int> &energyMap2 )
{
  energyMap.resize(w*h,-1);
  energyMap2.resize(w*h,-1);
  Scalar energy(0);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;

    int prop0   = segImg[ pixId ];
    int pix2_0  = currentSolutionE.seg2segFW[ pixId ].second;// projects onto

    if (pix2_0 <0) // maps on no other segment:
      energyMap[i] = -1;
    else
    {
      int prop00 = segImg2[pix2_0];

      if( prop00 == prop0 ) // normal data term
        energyMap[i] = 0;
      else // must check depth in both versions:
        energyMap[i] = whichPenalty ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop00].iNom, P3FromGlobal(pixId), dataElements[prop00].iHom, P3FromGlobal(pix2_0), prop00, prop0 );
    }
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;

    int prop0   = segImg2[ i ];
    int pix2_0  = currentSolutionE.seg2segBW[ i ].second;// 0 in view2 maps to this segment

    if (pix2_0 < 0 )
      energyMap2[i] = -1;
    else
    {
      int prop00 = segImg[ pix2_0 ];

      if( prop00 == prop0 ) // normal data term
        energyMap2[i] = 0;
      else // must check depth in both versions:
        energyMap2[i] = whichPenalty ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop00].Nom, P3FromGlobal(i), dataElements[prop00].Hom, P3FromGlobal( pix2_0), prop00, prop0 );
    }///////////////////////////////////////////
  }/////////////////
#endif
}

template<class Scalar>
Scalar Datasegments<Scalar>::
getEnergyPixel( ) 
{
  int dummy(0);
  Scalar energy(0);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;

    int prop0   = segImg[ pixId ];
    int pix2_0  = currentSolutionE.seg2segFW[ pixId ].second;// projects onto

#ifdef _simple_data_
    pix2_0 = -1;
#endif
    if (pix2_0 <0) // maps on no other segment:
      energy +=  currentSolutionE.dataFW[ pixId ];
    else
    {
      int prop00 = segImg2[pix2_0];

      if( prop00 == prop0 ) // normal data term
        energy += currentSolutionE.dataFW[ i ];
      else // must check depth in both versions:
        energy += getPenalty ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop00].iNom, P3FromGlobal(pixId), 1, //P3( (pixId/h) +1., (pixId%h) +1., 1.), 1,
        currentSolutionE.dataFW[ pixId ], currentSolutionE.freedataFW[ pixId ], 1, dummy, 1. );
    }
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;

    int prop0   = segImg2[ i ];
    int pix2_0  = currentSolutionE.seg2segBW[ i ].second;// 0 in view2 maps to this segment

#ifdef _simple_data_
    pix2_0 = -1;
#endif
    if (pix2_0 < 0 )
      energy += currentSolutionE.dataBW[ i ];
    else
    {
      int prop00 = segImg[ pix2_0 ];

      if( prop00 == prop0 ) // normal data term
        energy += currentSolutionE.dataBW[ i ];
      else // must check depth in both versions:
        energy += getPenalty ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop00].Nom, P3FromGlobal(i), 1,// P3( (i/h) +1., (i%h) +1., 1.), 1,
        currentSolutionE.dataBW[ i ], currentSolutionE.freedataBW[ i ], 1, dummy, 1.);
    }///////////////////////////////////////////
  }/////////////////
#endif
  return energy;
}

/// compare the solutions from all 0's and the one with 1's set as suggested - where are the differences?
template<class Scalar>
Scalar
  Datasegments<Scalar>::
  getEnergyPixelLocal( std::vector<int>& solution01, int trial ) 
{
  // at first construct a mapping from pixel to variables! and back. s.t. I know when a flipped pixel is involved
  std::vector<int> changedPixel (w*h, 0);
  std::vector<int> changedPixel2(w*h, 0);
  std::vector<int> changedPixelList;
  for ( int i=0;i< newElem.seg2segFW.size(); i++ )
  {
    if (newElem.seg2varFW[i] >=0 && solution01[ newElem.seg2varFW[i] ] )
    {
      changedPixel [newElem.seg2segFW[i].first] = 1;
      changedPixelList.push_back( newElem.seg2segFW[i].first );
    }
  }
  for ( int i=0;i< newElem.seg2segBW.size(); i++ )
  {
    if (newElem.seg2varBW[i] >=0 && solution01[ newElem.seg2varBW[i] ] )
    {
      changedPixel2 [newElem.seg2segBW[i].first] = 1;
      changedPixelList.push_back( newElem.seg2segBW[i].first );
    }
  }
  ///////////////////////////////////////////////////////////////////////////////
  Scalar energy(0), energyA(0), energy0(0), energy1(0);
  ////////////////
  // part one: left to right:
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;
    int whatcase=0;
    int prop0   = segImg[ pixId ];
    int pix2_0  = currentSolutionE.seg2segFW[ pixId ].second;// projects onto

    if (! ( changedPixel[i] || ( (pix2_0>0) && changedPixel2[pix2_0] ) ) )
      continue;

#ifdef _simple_data_
    pix2_0 = -1;
#endif
    ///// energy 0:
    if (pix2_0 <0) // maps on no other segment:
      energy =  currentSolutionE.dataFW[ pixId ];
    else
    {
      int prop00 = segImg2[pix2_0];

      if( prop00 == prop0 ) // normal data term
        energy = currentSolutionE.dataFW[ i ];
      else // must check depth in both versions:
        energy = getPenalty ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop00].iNom, P3FromGlobal(pixId), 1, //P3( (pixId/h) +1., (pixId%h) +1., 1.), 1,
        currentSolutionE.dataFW[ pixId ], currentSolutionE.freedataFW[ pixId ], 1, 1. );
    }
    energy0 += energy;
    ///////////////////////// part2: oh oh, cases 01,10,11 !!!!
    // 3 cases!
    if (changedPixel[i]) // cases 10, 11
    {
      // first find pixel we project onto:
      int boxId = newElem.getBoxIdFW(i);
      assert( newElem.seg2segFW[boxId].first == i);
      int pixel2 = newElem.seg2segFW[boxId].second;

      if ( changedPixel2[ pixel2 ] )
        energyA = newElem.dataFW[boxId];//case 11
      else // case 10:
      {
        int prop10 = segImg2[pixel2];whatcase=1;
        energyA = getPenalty ( newElem.Hom, newElem.iNom, dataElements[prop10].iNom, P3FromGlobal(pixId), 1,
          newElem.dataFW[ boxId ], newElem.freedataFW[ boxId ], 1, 1. );
        if (energy < energyA)
          int breakHere =1;
      }
    }
    else // case 01 since changedPixel2[pix2_0]==1;
    {
      assert (pix2_0 >=0 && changedPixel2[pix2_0]==1);
      int boxId = newElem.getBoxIdFW(i);whatcase=2;
      energyA = getPenalty ( dataElements[prop0].Hom, dataElements[prop0].iNom, newElem.iNom, P3FromGlobal(pixId), 1,
        currentSolutionE.dataFW[ pixId ], currentSolutionE.freedataFW[ pixId ], 1, 1. );
      if (energy < energyA)
        int breakHere =1;
    }
    energy1 += energyA;

    if (energy < energyA)
      int breakHere =1;
  }
#ifndef _FW_ONLY_
  ///////////////// BW: ////////////////////////
  for(int i= 0; i < w*h; i++ )
  {
    int pixId = i;
    int whatcase=0;
    int prop0   = segImg2[ i ];
    int pix2_0  = currentSolutionE.seg2segBW[ i ].second;// 0 in view2 maps to this segment

    if (! ( changedPixel2[i] || ( (pix2_0>0) && changedPixel[pix2_0] ) ) )
      continue;


#ifdef _simple_data_
    pix2_0 = -1;
#endif
    if (pix2_0 < 0 )
      energy = currentSolutionE.dataBW[ i ];
    else
    {
      int prop00 = segImg[ pix2_0 ];

      if( prop00 == prop0 ) // normal data term
        energy = currentSolutionE.dataBW[ i ];
      else // must check depth in both versions:
        energy = getPenalty ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop00].Nom, P3FromGlobal(i), 1,// P3( (i/h) +1., (i%h) +1., 1.), 1,
        currentSolutionE.dataBW[ i ], currentSolutionE.freedataBW[ i ], 1, 1.);
    }///////////////////////////////////////////
    energy0 += energy;
    // 3 cases now

    ///////////////////////// part2: oh oh, cases 01,10,11 !!!!
    // 3 cases!
    if (changedPixel2[i]) // cases 10, 11
    {
      // first find pixel we project onto:
      int boxId = newElem.getBoxIdBW(i);
      assert( newElem.seg2segBW[boxId].first == i);
      int pixel2 = newElem.seg2segBW[boxId].second;

      if ( changedPixel[ pixel2 ] )
        energyA = newElem.dataBW[boxId];//case 11
      else // case 10:
      {
        int prop10 = segImg[pixel2];whatcase=1;
        energyA = getPenalty ( newElem.iHom, newElem.Nom, dataElements[prop10].Nom, P3FromGlobal(pixId), 1,
          newElem.dataBW[ boxId ], newElem.freedataBW[ boxId ], 1, 1. );
        if (energy < energyA)
          int breakHere =1;
      }
    }
    else // case 01 since changedPixel2[pix2_0]==1;
    {
      assert (pix2_0 >=0 && changedPixel[pix2_0] ==1);
      int boxId = newElem.getBoxIdBW(i);whatcase=2;
      energyA = getPenalty ( dataElements[prop0].iHom, dataElements[prop0].Nom, newElem.Nom, P3FromGlobal(pixId), 1,
        currentSolutionE.dataBW[ pixId ], currentSolutionE.freedataBW[ pixId ], 1, 1. );
      if (energy < energyA)
        int breakHere =1;
    }
    energy1 += energyA;

    if (energy < energyA)
      int breakHere =1;

  }/////////////////
#endif

  Scalar diff = energy1-energy0;
  return diff;
}


/// needs as input only the segment mapping, which I have here, should also prepare zbuffer
template<class Scalar>
void
  Datasegments<Scalar>::
  buildFromCurrentSolutionPixel( )
{
  /// the current solution and its pixel positions
  std::vector< Scalar > sol_pix_idxFW; 
  std::vector< Scalar > sol_pix_idyFW; 
  std::vector< Scalar > sol_pix_idxBW; 
  std::vector< Scalar > sol_pix_idyBW;

  /////////////////
  sol_pix_idxFW.clear();sol_pix_idxFW.resize(w*h);
  sol_pix_idyFW.clear();sol_pix_idyFW.resize(w*h);
  sol_pix_idxBW.clear();sol_pix_idxBW.resize(w*h);
  sol_pix_idyBW.clear();sol_pix_idyBW.resize(w*h);
  // fill these: later do more like create a z-buffer 

  std::vector<int> segmentInvolved1( w*h, 1);// all are involved

  std::vector< std::pair<int,int> > segmentsareas1;    segmentsareas1.reserve(w*h);
  std::vector< std::pair<int,int> > segmentsareas2;    segmentsareas2.reserve(w*h);
  /////////////

  int innerId=0;
  for (int i=0; i < w; i++ )
    for (int j=0; j < h; j++, innerId++)
    {
      int pixId = i*h +j;
      int pixId2 = -1;
      int prop0  =  segImg[pixId];
      // lookup partner:
      int partner = -1; // -1: no partner
      M3 Hom = dataElements[prop0].Hom;
      P3 p0 = Hom * P3( i+1, j+1, 1. );p0 /= p0[2];
      sol_pix_idxFW[innerId] = p0[0];sol_pix_idyFW[innerId] = p0[1];

      int inBounds = isinBounds(p0);
      if ( inBounds )
      {
        pixId2 = toPixId( p0 );
        partner =  segImg2[pixId2]; // thats the solution - i need the pixel
      }
      segmentsareas1.push_back( std::pair<int,int> (pixId, pixId2) );
    }

    innerId=0;
    for (int i=0; i < w; i++ )
      for (int j=0; j < h; j++, innerId++)
      {
        int pixId = i*h + j;
        int pixId2 = -1;
        int prop0  =  segImg2[pixId];
        // lookup partner:
        int partner = -1; // -1: no partner
        M3 Hom = dataElements[prop0].iHom;

        P3 p0 = Hom * P3( i+1, j+1, 1. );p0 /= p0[2];
        sol_pix_idxBW[innerId] = p0[0];sol_pix_idyBW[innerId] = p0[1];
        int inBounds = isinBounds(p0);
        if ( inBounds )
        {
          pixId2 = toPixId( p0 );
          partner =  segImg[pixId2]; // thats the solution - i need the pixel
          // since these parts are occluders its fine, but assume an occlusion here - so no corr in other image
          // now 
        }
        segmentsareas2.push_back( std::pair<int,int> (pixId, pixId2) );
      }
      currentSolutionE.seg2segFW = segmentsareas1;
      currentSolutionE.seg2segBW = segmentsareas2;
      ///////////////////

      gws1.warp_noOmp_patchBased( sol_pix_idxBW, sol_pix_idyBW, currentSolutionE.seg2segFW );// actually should be increased by # of census border-pixel
      gws2.warp_noOmp_patchBased( sol_pix_idxFW, sol_pix_idyFW, currentSolutionE.seg2segBW );

      /// !!!!!!!! must be ,h,w not other way round !!!!
      //      per pixel stuff: ids maps to global ids. lIds: local ids, just -1: ignore, consider if >=0 -- saves into locScores: value per pixel in box!, global ids are now: seg2var.first
      gs1.computeFullScoresMappedCensusNew_short( gws1.getOrigShortImage(), gws2.getWarpedShortImage(), gws2.getIdx(), gws2.getIdy(), //&(sol_pix_idxFW[0]), &(sol_pix_idyFW[0]), 
        currentSolutionE.seg2segFW, segmentInvolved1, h,w, occlusions1 );
      gs2.computeFullScoresMappedCensusNew_short( gws2.getOrigShortImage(), gws1.getWarpedShortImage(), gws1.getIdx(), gws1.getIdy(), //&(sol_pix_idxBW[0]), &(sol_pix_idyBW[0]), 
        currentSolutionE.seg2segBW, segmentInvolved1, h,w, occlusions2 );

      currentSolutionE.dataFW = gs1.getScores() + gs1.getPenScores();// penalty is data - free + occlusion if in image!?
      currentSolutionE.dataBW = gs2.getScores() + gs2.getPenScores();
      currentSolutionE.freedataFW = gs1.getScores();
      currentSolutionE.freedataBW = gs2.getScores();

      currentSolutionE.freeFW.assign( currentSolutionE.dataFW.size() , 1);
      currentSolutionE.freeBW.assign( currentSolutionE.dataBW.size() , 1);

#ifdef __depthControl__
    for (int i=0; i < w; i++ )
      for (int j=0; j < h; j++)
      {
        int pixId = i*h + j;
        int prop0FW  =  segImg [pixId];
        int prop0BW  =  segImg2[pixId];
        // lookup partner:
        int partner = -1; // -1: no partner

        P3 nom  = dataElements[ prop0FW ].Nom;
        P3 inom = dataElements[ prop0BW ].iNom;

        Scalar dfw = Scalar(1.) / ( nom  | P3( i+1, j+1, 1. ));
        Scalar dbw = Scalar(1.) / ( inom | P3( i+1, j+1, 1. ));

        if ( dfw > Scalar(-_minDepth_))
        {
          currentSolutionE.dataFW[pixId]     = Scalar(__depthControlAmplifier__ * impFactor);
          currentSolutionE.freedataFW[pixId] = 0;
          currentSolutionE.freeFW[pixId]     = 0;
        }
        if ( dbw > Scalar(-_minDepth_))
        {
          currentSolutionE.dataBW[pixId]     = Scalar(__depthControlAmplifier__ * impFactor);
          currentSolutionE.freedataBW[pixId] = 0;
          currentSolutionE.freeBW[pixId]     = 0;
        }
     }
#endif
   
   Scalar sumError(0), sumVError(0);
      for (int i =0;i<currentSolutionE.dataFW.size(); i++)
      {
        sumError  += currentSolutionE.dataFW[i];sumVError += currentSolutionE.freedataFW[i];
        if (currentSolutionE.dataFW[i] != currentSolutionE.freedataFW[i])
          currentSolutionE.freeFW[i] =0;
      }
#ifdef _writeDebugOut_
      printf("--- Data-Error    %f - %f \n", sumError, sumVError);
#endif
      sumError= 0; sumVError = 0;
      for (int i =0;i<currentSolutionE.dataBW.size(); i++)
      {
        sumError  += currentSolutionE.dataBW[i];sumVError += currentSolutionE.freedataBW[i];
        if (currentSolutionE.dataBW[i] != currentSolutionE.freedataBW[i])
          currentSolutionE.freeBW[i] =0;
      }
#ifdef _writeDebugOut_
      printf("--- Data-Error v2 %f - %f \n", sumError, sumVError);
#endif
      //
      std::vector<M3> homsFW( dataElements.size() );
      std::vector<M3> homsBW( dataElements.size() );
      for (int i=0;i<dataElements.size(); i++)
      {
        homsFW[i] = dataElements[i].Hom;
        homsBW[i] = dataElements[i].iHom;
      }
#ifdef _expansion_ON_
      if ( doOccBuffering )
      {
        oMapBufFW.initOccBuffer( homsFW );
        oMapBufBW.initOccBuffer( homsBW );
      }
#endif
}


/// delivers the data term as edges and unaries, with the variable names given in the expansion set
template<class Scalar>
void
  Datasegments<Scalar>::
  expandRegionPixel( int trial )
{
  int edgeCase =0;

  const dataElem<Scalar> &temp = newElem;
  //    int trial = newElem.proposal;
  unaries.clear();
  binaries.clear();
  //    unaries.assign ( temp.varsFW + temp.varsBW, P2(0,0) );  //    unaries.assign ( temp.dataFW.size() + temp.dataBW.size(), P2(0,0) );  
  binaries.reserve( 4*(temp.seg2segFW.size() + temp.seg2segBW.size()) );//, Binary<Scalar> (-1,-1));
  unaries.reserve(  2*(temp.seg2segFW.size() + temp.seg2segBW.size()) );//, could be FW,BW and 0,1
  ////////////////
  //std::vector<int> edgeMap ((temp.varsFW) * (temp.varsBW), -1);

  std::vector< std::pair<int,int> > edgeStoreMap ((temp.varsFW), std::pair<int,int> (-1,-1));
  std::vector< std::pair<int,int> > edgeStoreMapB((temp.varsFW), std::pair<int,int> (-1,-1));

  // part one: left to right:
  const std::vector< std::pair<int,int> > &fwMap = temp.seg2segFW;
  for(int i= 0; i < fwMap.size(); i++ )
  {
    int varIdx  = temp.seg2varFW[ i ];// id of the first variable
    if (varIdx<0) continue;
    // 1. assign a label 0 to segi in left view:
    int segId = fwMap[i].first; // global pixel id
    // so. check onto which segment the segment projects to and check if it is considered
    // how? here:
    int prop0   = segImg[ segId ];
    //      if ( prop0 == trial ) continue; // no impact - right
    int prop1   = trial;
    int seg2_0  = currentSolutionE.seg2segFW[ segId ].second;// global pixel id, variable id by global2local mapping
    int seg2_1  = fwMap[i].second;                           // global pixel id, variable id by same index
    int varId0x = -2;//does not project to anything. 
    int varId1x = -2;

    if (seg2_0 > -1 )
      varId0x = temp.getLocalIdBW( seg2_0 );// variable id (can be -1 == not considered)
    if (seg2_1 > -1 )
      varId1x = temp.getLocalIdBW( seg2_1 );// should be >=0 - NO
    // these are the local variables. 
#ifdef _simple_data_
    varId0x =-2;
    varId1x =-2;
#endif

    if ( prop0 == trial && (varId0x<=-1 || segImg2[ seg2_0 ] == trial) ) continue; // no impact 

    P3 ownImageCoord;
    if (varId0x != -2 || varId1x!=-2 )
      ownImageCoord = temp.getImageCoordFW(i);

    // first addres case 00 and 01:
    // then binary, else 00 is a unary, 01 does not exist, no matter what 
    if (varId0x == -2) // maps on no other pixel:
      unaries.push_back( Unary<Scalar>(varIdx, P2( currentSolutionE.dataFW[ segId ], 0) )) ;
    else // segment exists but no variable for it:, so 01 does not exist
    {
      // only need to consider 00:
      // same segi? or dame depth? or larger or smaller -> 3 possibilities

      int prop00 = segImg2[seg2_0];
      Scalar penalty =0;
      if( prop00 == prop0 ) // normal data term
        penalty = currentSolutionE.dataFW[ segId ];
      else // must check depth in both versions:
      {
        penalty = getPenalty_x2 ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop00].iNom, seg2_0, ownImageCoord, 1,
          currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1.);
        penalty += addPotts;
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
          penalty = getPenalty_x2 ( dataElements[prop0].Hom, dataElements[prop0].iNom, dataElements[prop01].iNom, seg2_0, ownImageCoord, 1,
            currentSolutionE.dataFW[ segId ], currentSolutionE.freedataFW[ segId ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
          penalty += addPotts;
        }

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
      int prop10 = segImg2[seg2_1];
      if( prop10 == trial ) // normal data term
        penalty = temp.dataFW[ i ];
      else // must check depth in both versions:
      {
        penalty = getPenalty_x2 ( temp.Hom, temp.iNom, dataElements[prop10].iNom, seg2_1, ownImageCoord, 1,
          temp.dataFW[ i ], temp.freedataFW[ i ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
        penalty += addPotts;
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
        else // can not enter this case: -- IMPOSSIBLE FROM HERE
        {
          assert( 0 );
          penalty = getPenalty_x2 ( temp.Hom, temp.iNom, dataElements[prop11].iNom, seg2_1, ownImageCoord, 1,
            temp.dataFW[ i ], temp.freedataFW[ i ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[segId] : 1. );
          penalty += addPotts;
        }
        if (varId0x != varId1x)
        {
          binaries.push_back( Binary<Scalar> (varIdx, varId1x, P4 (0, 0, pen10, penalty ) ) );
          edgeStoreMapB[(varIdx - temp.firstIdFW)] = std::pair<int,int> (varId1x, binaries.size()-1);
        }
        else // double entry:
        {
          binaries[binaries.size()-1].edge += P4 (0, 0, pen10, penalty );
          //            edgeStoreMap[(varIdx - temp.firstIdFW)] = std::pair<int,int> (varId1x, binaries.size()-1);          
        }
      }
    }///////////////////////////////////////////
#ifdef _expansion_ON_
    expandRegionLongRangeInteractionsPixel( segId, varIdx, trial, 1 );
#endif
  }/////////////////
  ///////////////// BW: ////////////////////////
#ifndef _FW_ONLY_
  const std::vector< std::pair<int,int> > &bwMap = temp.seg2segBW;
  for(int i= 0; i < bwMap.size(); i++ )
  {
    int varIdx  = temp.seg2varBW[ i ];// id of the first variable
    if (varIdx<0) continue;
    // 1. assign a label 0 to segi in left view:
    int segId = bwMap[i].first;
    // so. check onto which segment the segment projects to and check if it is considered
    // how? here:
    int prop0   = segImg2[ segId ];
    int prop1   = trial;

    int seg2_0  = currentSolutionE.seg2segBW[ segId ].second;// 0 in view2 maps to this segment
    int seg2_1  = bwMap[i].second;                           // 0 in view2 maps to this segment
    int varId0x = -2;//does not project to anything. 
    int varId1x = -2;

#ifdef _breakCircles
    // check for pixel with label 0 its partner
    if (seg2_0 > -1  && currentSolutionE.seg2segFW[ seg2_0 ].second != segId )
      seg2_0 = -1;
#endif
    if (seg2_0 > -1 )
      varId0x = temp.getLocalIdFW( seg2_0 );// variable id (can be -1 == not considered)
    if (seg2_1 > -1 )
      varId1x = temp.getLocalIdFW( seg2_1 );// should be >=0
    // these are the local variables. 
#ifdef _simple_data_
    varId0x =-2;
    varId1x =-2;
#endif

    if ( prop0 == trial && (varId0x<=-1 || segImg[ seg2_0 ] == trial) ) continue; // no impact 

    // first addres case 00 and 01:
    // then binary, else 00 is a unary, 01 does not exist, no matter what 
    if (varId0x == -2) // maps on no other segment:
      unaries.push_back( Unary<Scalar>(varIdx, P2( currentSolutionE.dataBW[ segId ], 0) ));
    else // segment exists but no variable for it:, so 01 does not exist
    {
      int prop00 = segImg[seg2_0];
      Scalar penalty =0;
      if( prop00 == prop0 ) // normal data term
        penalty = currentSolutionE.dataBW[ segId ];
      else // must check depth in both versions:
      {
        penalty = getPenalty_x2 ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop00].Nom, seg2_0, temp.getImageCoordBW(i), 1,
          currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
        penalty += addPotts;
      }
      if ( varId0x == -1 ) // partner without variable index
        unaries.push_back( Unary<Scalar>(varIdx, P2( penalty, 0) ));
      else
      {
        Scalar pen00 = penalty;

        int prop01 = trial;
        Scalar penalty =0;
        if( prop01 == prop0 ) // normal data term
          penalty = currentSolutionE.dataBW[ segId ];
        else // must check depth in both versions
        {
          penalty = getPenalty_x2 ( dataElements[prop0].iHom, dataElements[prop0].Nom, dataElements[prop01].Nom, seg2_0, temp.getImageCoordBW(i), 1,
            currentSolutionE.dataBW[ segId ], currentSolutionE.freedataBW[ segId ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
          penalty += addPotts;
        }
        if ( edgeStoreMap[(varId0x - temp.firstIdFW)].first != varIdx ) 
        {
          if ( edgeStoreMapB[(varId0x - temp.firstIdFW)].first != varIdx )
            binaries.push_back( Binary<Scalar> (varId0x, varIdx, P4 (pen00, 0, penalty,0 ) ) );//new
          else
            binaries[ edgeStoreMapB[(varId0x - temp.firstIdFW)].second ].edge += P4 (pen00, 0, penalty,0 );//old
        }
        else
          binaries[ edgeStoreMap[(varId0x - temp.firstIdFW)].second ].edge += P4 (pen00, 0, penalty,0 );//old
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
      int prop10 = segImg[seg2_1];
      Scalar penalty =0;
      if( prop10 == trial ) // normal data term
        penalty = temp.dataBW[ i ];
      else // must check depth in both versions:
      {
        penalty = getPenalty_x2 ( temp.iHom , temp.Nom, dataElements[prop10].Nom, seg2_1, temp.getImageCoordBW(i), 1,
          temp.dataBW[ i ], temp.freedataBW[ i ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
        penalty += addPotts;
      }
      if ( varId1x == -1 ) // not encoded variable
        unaries.push_back( Unary<Scalar>(varIdx, P2( 0, penalty )) );
      else  // is encoded
      {
        Scalar pen10 = penalty;

        int prop11 = trial;
        Scalar penalty =0;
        if( prop11 == prop1 ) // normal data term
          penalty = temp.dataBW[ i ];//currentSolutionE.dataBW[ segId ];
        else // must check depth in both versions: does not happen -- NONSENSE FROM HERE
        {
          assert( 0 );
          penalty = getPenalty_x2 ( temp.iHom, temp.Nom, dataElements[prop11].Nom, seg2_1, temp.getImageCoordBW(i), 1,
            temp.dataBW[ i ], temp.freedataBW[ i ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[segId] : 1. );
          penalty += addPotts;
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
          binaries[ edgeStoreMap[(varId1x - temp.firstIdFW)].second ].edge += P4 (0, pen10, 0, penalty );
      }
    }
#ifdef _expansion_ON_
    expandRegionLongRangeInteractionsPixel( segId, varIdx, trial, 0 );
#endif
  }/////////////////
#endif
}

//  void updateSolutionsPixel( int trial, std::vector<int>& _currentSolution1, std::vector<int>& _currentSolution2, const std::vector<int>& solution01 )
/// must also update the zbuffer information
template<class Scalar>
void
  Datasegments<Scalar>::
  updateSolutionsPixel( int trial, const std::vector<int>& solution01, int fromPotElem )
{
  //    printf("Update solution: %d prop Now: %d\n", fromPotElem, storedElem[fromPotElem].proposal );

  assert(fromPotElem  < storedElem.size());
  const dataElem<Scalar> &temp = (fromPotElem < 0) ? newElem : storedElem[fromPotElem] ;
  const std::vector< std::pair<int,int> > &fwMap = temp.seg2segFW;// its pixel to pixel
  const std::vector< std::pair<int,int> > &bwMap = temp.seg2segBW;
  //
  // run along current solution, if something changed, update it: no 
  // run along last trial 
  for (int i=0;i< fwMap.size(); i++)
  {
    int localId = temp.seg2varFW[ i ] ;
    if (localId<0) continue;

    if ( solution01[ localId ] )
    {
      Scalar CSdataFW_curr = currentSolutionE.dataFW[ fwMap[i].first ];
      Scalar CSdataFW_temp = temp.dataFW[i];

      if (CSdataFW_curr < CSdataFW_temp)
        CSdataFW_temp = CSdataFW_temp;  
      currentSolutionE.dataFW[ fwMap[i].first ]     = temp.dataFW[i];
      currentSolutionE.freedataFW[ fwMap[i].first ] = temp.freedataFW[i];
      currentSolutionE.seg2segFW[  fwMap[i].first ] = temp.seg2segFW[i];//pixel to pixel mapping

      segImg[ fwMap[i].first ] = trial;
#ifdef _expansion_ON_
      oMapBufFW.update( fwMap[i].first, temp.Hom );
#endif
    }
  }
#ifndef _FW_ONLY_
  for (int i=0;i < bwMap.size(); i++)
  {
    int localId = temp.seg2varBW[ i ] ;
    if (localId<0) continue;

    if ( solution01[ localId ] )
    {
      Scalar CSdataBW_curr = currentSolutionE.dataBW[ bwMap[i].first ];
      Scalar CSdataBW_temp = temp.dataBW[i];

      if (CSdataBW_curr < CSdataBW_temp)
        CSdataBW_temp = CSdataBW_temp;

      currentSolutionE.dataBW[ bwMap[i].first ] = temp.dataBW[i];
      currentSolutionE.freedataBW[ bwMap[i].first ] = temp.freedataBW[i];
      currentSolutionE.seg2segBW[ bwMap[i].first ]  = temp.seg2segBW[i];//pixel to pixel mapping
      segImg2[ bwMap[i].first ] = trial;
#ifdef _expansion_ON_
      oMapBufBW.update( bwMap[i].first, temp.iHom );
#endif
    }
  }
#endif

  assert(fromPotElem  < storedElem.size());
  if (fromPotElem >= 0) storedElem[fromPotElem].proposal =-1;
}


/*! idea already knowledge of the pixel, call right there in expansion pixel, 
* find pixel in other view which might be occluded by this pixel.
* Example for all FW pixel with a variable we need to check if a pixel 
* in the other view not in the variable list falls onto it.
* If so, add a nuary term, if assigning label0/1 changes the score of that pixel.
* so: e.g. FW variable varIdx is pixel pixId. Check all pixel in 2nd view (BW).
* which fall onto this one (pixId). Use OmapBufBW! to get all pixel in BW view falling
* onto it.
*/
template<class Scalar>
void
  Datasegments<Scalar>::
  expandRegionLongRangeInteractionsPixel( int pixId, int varIdx, int trial, int FW )
{
  int edgeCase(0);
  if (FW)
  {
    int prop00 = segImg[pixId];//pixId in FW map, find all bw pixel falling onto it
    const typename std::vector< int >& potOList = oMapBufBW.get_Partners( pixId );
    for (int i=0; i<potOList.size(); i++ )
    {
      int id = potOList[i];
      int prop = segImg2[ id ];

      int varId0 = newElem.getLocalIdBW( id );//variable id of the occluder:
      if (varId0 >= 0) continue;

#ifdef _DEBUG
      // debug test:
      P3 dest = dataElements[prop].iHom * P3FromGlobal(id);dest /= dest[2];// P3((id/h) + 1., (id%h) + 1., 1.);dest /= dest[2];
      //        int projectsTo = dest[0]*h + dest[1];
      int projectsTo = toPixId(dest);//int( dest[0]-0.5 )*h +int ( dest[1]-0.5 );
      assert(projectsTo == pixId);
#endif
      // find both penalties, map the BW pixel to FW view, check if occluded/or other stuff
      Scalar pen00 = getPenalty_x2 ( dataElements[prop].iHom, dataElements[prop].Nom, dataElements[ prop00 ].Nom, pixId, P3FromGlobal(id), 1, //P3((id/h) + 1., (id%h) + 1., 1.), 1,
        currentSolutionE.dataBW[ id ], currentSolutionE.freedataBW[ id ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[id] : 1.) + (prop != prop00) * addPotts;
      Scalar pen01 = getPenalty_x2 ( dataElements[prop].iHom, dataElements[prop].Nom, newElem.Nom, pixId, P3FromGlobal(id),1,//P3((id/h) + 1., (id%h) + 1., 1.), 1,
        currentSolutionE.dataBW[ id ], currentSolutionE.freedataBW[ id ], 1, edgeCase, autoScoresBW ? adaptiveAutoScoresBW[id] : 1.) + (prop != trial) * addPotts;
      if(pen00 !=pen01)
        unaries.push_back( Unary<Scalar>(varIdx, P2( pen00, pen01) ));
    }
  }
  else
  {
#ifndef _FW_ONLY_
    int prop00 = segImg2[pixId];//pixId in BW map, find all pixel falling onto it
    const typename std::vector< int >& potOList = oMapBufFW.get_Partners( pixId );
    for (int i=0; i < potOList.size(); i++ )
    {
      int id = potOList[i];
      int prop = segImg[ id ];

      int varId0 = newElem.getLocalIdFW( id );//variable id of the occluder:
      if (varId0 >= 0) continue;

#ifdef _DEBUG
      // debug test:
      P3 dest = dataElements[prop].Hom * P3FromGlobal(id);dest /= dest[2];// P3((id/h) + 1., (id%h) + 1., 1.);dest /= dest[2];
      //        int projectsTo = dest[0]*h + dest[1];
      int projectsTo = toPixId( dest );//int( dest[0]-0.5 )*h +int ( dest[1]-0.5 );
      assert(projectsTo == pixId);
#endif
      // find both penalties
      Scalar pen00 = getPenalty_x2 ( dataElements[prop].Hom, dataElements[prop].iNom, dataElements[ prop00 ].iNom, pixId, P3FromGlobal(id),1,//P3((id/h) + 1., (id%h) + 1., 1.), 1,
        currentSolutionE.dataFW[ id ], currentSolutionE.freedataFW[ id ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[id] : 1. ) + (prop != prop00) * addPotts;
      Scalar pen01 = getPenalty_x2 ( dataElements[prop].Hom, dataElements[prop].iNom, newElem.iNom, pixId, P3FromGlobal(id),1,//P3((id/h) + 1., (id%h) + 1., 1.), 1,
        currentSolutionE.dataFW[ id ], currentSolutionE.freedataFW[ id ], 1, edgeCase, autoScoresFW ? adaptiveAutoScoresFW[id] : 1. ) + (prop != trial) * addPotts;
      if(pen00 !=pen01)
        unaries.push_back( Unary<Scalar>(varIdx, P2( pen00, pen01) )) ;
    }
#endif
  }
}


/// does it use [(R|t) * (P3d,1)] | (N_new,-1) = 0 - must be fulfilled: p^t * (R|t)^t * (R|t)^-t * (N_old,-1) = p^t * N_old = 0 
template<class Scalar>
Scalar Datasegments<Scalar>::
getPenalty_x2 ( const M3& Hom, const P3& Nom, const P3& iNom, int pixId, P3 centre, 
  int area, Scalar data, Scalar freeData, Scalar freePix, int& edgeCase, Scalar autoScore)
{
  P3  c2 = P3FromGlobal( pixId );
  edgeCase = 0;

  if (simpleApproach)
  {
    edgeCase = 1;
    return data;
  }
#ifdef _simple_data_
  return data; // what happens? : not working as well, but no energy increase
#endif
  Scalar localOccPen = xtraPen+autoScore;// ?? 0.2 -> TODO ALSO FOR OOB
  P3 p = Hom * centre;p /= p[2];//projected point

  // on/off ??
  //    p[0] = trunc(p[0]+0.5);// change No 1
  //    p[1] = trunc(p[1]+0.5);

  Scalar depth0      = 1./(p|Nom); 
  Scalar depthThresh = max ( maxDepthThresh, fabs(depth0  * __OCC_THRESH__) );
  Scalar depthDiff   = depth0 - 1./(p|iNom);// let depths by positive -> depthDiff>0 implies an occlusion case, sign change: < 0 occlusion case

  // OOB case - or in front of image plane
  if ( depth0 > Scalar(-_minDepth_) )
#ifdef _inFrontPenalty_
    return max(data, data-freeData + area * Scalar(impFactor)*localOccPen) - area*addPotts; // in front of cam : impossible
#else
    return data-freeData + freePix*localOccPen - area*addPotts; // in front of cam: out of bounds
#endif

  if ( depthDiff > depthThresh ) // projection occludes other point but the other is visible: impossible
  {
    if ( impossibleApproach )
      return max(data, data-freeData + area * Scalar(impFactor)*localOccPen) - area*addPotts; // impFactor not too large: deadlocks
    else
      return max(data, data-freeData + area * localOccPen) - area*addPotts; // ok that is a bad idea need the extra treatment above
  }
  else
  {
    if ( depthDiff < -depthThresh ) // projection is occluded by other point: occlusion model
      return data-freeData + freePix*localOccPen - area*addPotts; // occlusion penalty

    edgeCase = 1;// normal case
    return data;
  }
}

#endif