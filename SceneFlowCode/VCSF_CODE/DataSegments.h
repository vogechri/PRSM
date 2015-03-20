/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

This software is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this software.
*/

///////////////////////////////////////////////////////////
//////      Computing the (binary) data score       ///////
///////////////////////////////////////////////////////////


#ifndef __DATASegments__h
#define __DATASegments__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "Math/Mat2x2T.h"

#include <algorithm>
//#include <functional>

#include "DataDefinitionsVC.h"

#include "genHom.h"
#include "genWarp.h"
#include "genScore.h"
#include "OcclusionExpansionBuffer.h"
#include "DataContainers.h"
using namespace std;
using namespace Math;


template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.resize(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   result.begin(), std::plus<T>());
    return result;
}

///////////////////////////////////////
/*! Provides the binary data term and associations of variables
in expansion areas and segments/pixels.
Keeps track of the current solution.
Provides multi-processor operations (for parallel evaluation)
*/
template<typename Scalar> class Datasegments
{
public:

  typedef Math::Mat2x2T<Scalar>     M2;
  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 4>  P4;
  typedef Math::VectorT<int, 4>     P4i;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
//  typedef Math::VectorT<Scalar, 5>  P5;

  typedef std::pair <int,int>       Seg_ij;

  Datasegments() : segImg(NULL), segImg2(NULL), nSegments(0), nSegments2(0), 
      w(0), h(0), occPenalty(0.8), impPenalty(8),
      centers(NULL), centers2(NULL), expCenters(NULL), nExpCenters(0), gws1( 0,0, NULL ), gws2( 0,0, NULL ), 
      gs1(0,0), gs2(0,0), occlusions1(NULL), occlusions2(NULL), oMapBufFW(0,0,NULL), oMapBufBW(0,0,NULL), 
      simpleApproach(false), impossibleApproach(true), doOccBuffering(false),
      addPotts(0.25), maxDepthThresh(0.1), impFactor(_impFactor_),
      computeMotionPairs(false), scale3D(1), __OCC_THRESH__(0.1), resizer(2.5),
      autoScoresFW(NULL), autoScoresBW(NULL), ctct(0,0,0,0), xtraPen(0.1), gridSize(16), maxMotSqr(1000*1000)
  {};

  Datasegments( Scalar* img1, Scalar *img2, std::vector<bool>& _oobOutside1, std::vector<bool>& _oobOutside2 )
    : segImg(NULL), segImg2(NULL), nSegments(0), nSegments2(0), 
      w(0), h(0), occPenalty(0.8), impPenalty(8),
      centers(NULL), centers2(NULL), expCenters(NULL), nExpCenters(0), gws1( 0,0, img1 ), gws2( 0,0, img2 ), 
      gs1(0,0), gs2(0,0), occlusions1(_oobOutside1), occlusions2(_oobOutside2), 
      oMapBufFW(0,0,NULL), oMapBufBW(0,0,NULL), 
      simpleApproach(false), impossibleApproach(true), doOccBuffering(false),
      addPotts(0.25), maxDepthThresh(0.1), impFactor(_impFactor_), 
      computeMotionPairs(false), scale3D(1), __OCC_THRESH__(0.1), resizer(2.5),
      autoScoresFW(NULL), autoScoresBW(NULL), ctct(0,0,0,0), xtraPen(0.1), gridSize(16), maxMotSqr(1000*1000)
  {};


  Datasegments(int _w, int _h, Scalar dtaThresh, Scalar oobThresh, 
    Scalar* img1, Scalar *img2, int* _segImg, int _nSegments, int* _segImg2, int _nSegments2, 
    std::vector<bool>& _oobOutside1, std::vector<bool>& _oobOutside2, bool _doOccBuffering=true )
    : segImg(_segImg), segImg2(_segImg2), nSegments(_nSegments), nSegments2(_nSegments2), 
      w(_w), h(_h), occPenalty(0.8), impPenalty(8),
      centers(NULL), centers2(NULL), expCenters(NULL), nExpCenters(0), gws1( _w,_h, img1 ), gws2( _w,_h, img2 ), 
      gs1(_w,_h, dtaThresh, oobThresh ), gs2(_w,_h, dtaThresh, oobThresh ),
      occlusions1(_oobOutside1), occlusions2(_oobOutside2),  
      oMapBufFW(_w,_h,_segImg), oMapBufBW(_w,_h,_segImg2), 
      simpleApproach(false), impossibleApproach(true), doOccBuffering(_doOccBuffering),
      addPotts(0.25), maxDepthThresh(0.1), impFactor(_impFactor_), 
      computeMotionPairs(false), scale3D(1), __OCC_THRESH__(0.1), resizer(2.5),
      autoScoresFW(NULL), autoScoresBW(NULL), ctct(0,0,0,0), xtraPen(0.1), gridSize(16), maxMotSqr(_w*_h)
  {
      gs1.setRefImage( img1 );
      gs2.setRefImage( img2 );
      setAreas();
  };

  ~Datasegments(){};

  void setOccThresh(Scalar occThresh){__OCC_THRESH__ = occThresh;};
  void setPotts (Scalar _addPotts)  {addPotts = _addPotts;};
//  void setMVPDistance (Scalar _dst) {mvp_dst = _dst*_dst;};
//  void setMNDistance (Scalar _dst)  {maxNormDiff = _dst;};
  void setScale3D(Scalar _val)      {scale3D = _val;};

  void setMotionPairs(bool _onOff) {computeMotionPairs = _onOff;};

  void setSimple(bool onOff = true) {simpleApproach = onOff;};
  void setImpossible(bool onOff = true) {impossibleApproach = onOff;};
  void setOccBuffering(bool onOff = true) {doOccBuffering = onOff;};

  /// info about the camera pair considered
  void setTimeCam(P4i ct_ct_) { ctct = ct_ct_;};
  /// return info about the camera pair considered
  P4i  getTimeCam()           { return ctct;};

  /// the edge length of a grid cell
  void setGridSize(int sz) {gridSize = sz;}

  void set_AutoScores(Scalar* _autoScores1, Scalar* _autoScores2);

  void updateAutoScores( );
  void averageDataScore( );

  void set_AutoScoresPix(Scalar* _autoScores1, Scalar* _autoScores2);

  void updateAutoScoresPix( );

  /// create debug output for visualization
  void buildDepthMap2( std::vector<Scalar>& depthMap, bool remap = false );
  /// create debug output for visualization
  void buildDepthMap( std::vector<Scalar>& depthMap, bool remap = false );

  void setExpCenters (Scalar* expCenters_, int nExpCenters_) { expCenters = expCenters_;nExpCenters = nExpCenters_;}

  void setCentersKs( Scalar* Kl_, Scalar* Kl2_, Scalar* _centers1, Scalar* _centers2, int nCenters=0 );

  P4i getBox1() {return box1;};
  P4i getBox2() {return box2;};

  P4i getBox1(int pid) {return newElemV[pid].boxFW;};
  P4i getBox2(int pid) {return newElemV[pid].boxBW;};

  const std::vector<P3>* getvNormals( ) {return &vNoms_cam;}; // applies on normal for second view
  const std::vector<M3>& getvHoms( )    {return  homs_pix;};   // trafo for first view into second
  const std::vector<M3>& getviHoms( )   {return ihoms_pix;};   // trafo for second view into first
  const std::vector<P3>& getpNoms( )    {return  noms_pix;};   // trafo for first view into second
  const std::vector<P3>& getpiNoms( )   {return inoms_pix;};   // trafo for second view into first


  void setOcclusionMaps( const std::vector<bool>& occlusionsR1, const std::vector<bool>& occlusionsR2 )
  {
    occlusions1 = occlusionsR1;occlusions2 = occlusionsR2;
  }

  /// tough too large bad results, too small takes ages!
  void setOccPen (Scalar _pen) {occPenalty = _pen;impPenalty = occPenalty;};
  void setMaxDisp(Scalar maxDisp_) {gs1.setMaxDisp(maxDisp_);gs2.setMaxDisp(-maxDisp_);maxMotSqr = maxDisp_*maxDisp_;};
  void setMaxMot (Scalar maxMot_)  {gs1.setMaxMot(maxMot_);gs2.setMaxMot(maxMot_);maxMotSqr = maxMot_*maxMot_;};

  /// could increase that over the iterations, to allow for soft optimization first (no deadlocks any more for high penalties)
  void setImpPen (Scalar _pen) {impPenalty = _pen;}//+addPotts;};//

  void setImpFactor (Scalar _impFactor) {impFactor = _impFactor;}
  void incImpFactor () {impFactor += 1;}

  /// generate list of 2d segment centers -- for test whether segment in expansion area 
  void init()
  {
    buildList ( Kl, centers, nSegments, projPatchCenters );
    buildList ( Kl2, centers2, nSegments2, projPatchCenters2 );
  }
  ////////////////////////////////////////////////////////////////////////////
  // functions which are called several times
  ////////////////////////////////////////////////////////////////////////////

  const dataElem<Scalar>& getElement    ()                            {return newElem;};
  const dataElem<Scalar>& getElement    (int pid)                     {return newElemV[pid];};

  const dataElem<Scalar>& getCurrentSol ()                            {return currentSolutionE;};
  const typename std::vector < dataElem<Scalar> >& getDataElements () {return dataElements;};


  int getElemVarsFirst(){return newElem.seg2segFW.size();};
  int getElemVarsPixel(){return newElem.varsFW;};

  int getElemVarsFirst(int pid){return newElemV[pid].seg2segFW.size();};
  int getElemVarsPixel(int pid){return newElemV[pid].varsFW;};

  /// init homographies and normals per view
  void initHoms( int normalId, genHomoG<Scalar>* gHom, const std::vector<P3>* normals );

  /// init homographies and normals  of pairs not involving the canonical view - using the output form these pairs
  void initHoms( int normalId, const std::vector<M3>& ihoms1, const std::vector<M3>& homs2, 
                 const std::vector<P3>& inoms1, const std::vector<P3>& inoms2 );

  // idea use a mapping from view ref to view c and ref to view b (b is neighbor of ref, c is not directly)
  // this already defines bboxes and involved segments.
  // some of these might not have a partner !! i need to model this ? or not 
  // what if i project these as now forward backward.
  // later consolidate the variable ids by going through map and adjusting seg2var mapping to be consistent,
  // this means that there can be jumps! not that unused segments receive a variable, they remain at -1
  // the variables can just be 1,2,5,6,7 in one viewpair.
  // i need a uniform mapping to variable ids from 2 views. i get 2 lists of involved segments from view a and view b. 


  /*! group proposal areas together, reduces time if proposals occur several times
  idea is t remove double occurences of normal/rotation pairs and join them make size dependent on variable size if joined
  */
  void reorganizeProposals (std::vector<int>& prop2Plane, int& nProposalExpansions, 
                            Scalar* expCenters, std::vector<Scalar>& propCenters,  //std::vector<dataElem<Scalar>::P3>& propCenters, 
                            const int boxX, const int boxY, std::vector<typename dataElem<Scalar>::P2i>& patchXY,
                            bool doNormalCheck = true);

  /*! 
  This function is to be called from pairs of the canonical to all other views.
  Thus we assign variable ids to all segments involved in the expansion areas - these must be unique.
  computes the expansion area in the other view(s) and other parts needed to compute photometric error.
   Also ids of variables are computed, needed for optimization and assigment variable,segment
   in : segi and proposals to expand, max box width, height
   result: list of segments involvde, 2 boxes in both views the warps have to be evaluated from. 
  */
  int initNewSeg( int propId, int proposal, int boxX, int boxY, genHomoG<Scalar>* gHom, int expansionArea, const std::vector<P3>* normals, int addToSeg2Vars = 0, int addToAllSeg2Vars = 0, int pid =-1 );

  /*! 
  Defines the mapping from segments to variables. the areas to 
  expand, homographies, normals per view, etc.
  Uses the output from initNewSeg - from other view pairs. 
  */
  void initNewSegFromView( int propId, int proposal, const P4i& _box1, const P4i& _box2, int expansionArea, const dataElem<Scalar>& leader, const dataElem<Scalar>& second, int pid =-1);

  /// release storage before per segment photometric data term is evaluated
  void requestNElements( int nProcs, int nProps );

  /// release storage after per segment photometric data term was evaluated completely
  void releaseNElements();

  /// compute and store photometric data term per segment
  void createWarps( int pid = -1, int propId = -1);

  /// evaluate photometric data cost for the area in the box, needs newElemV (or Hom, iHom) to be initialized
  void createWarps( const P4i& bigbox1, const M3& Hom, const P4i& bigbox2, const M3& iHom, int pid =-1, int propId = -1 );

  /// OBSOLETE idea: for all segments, get (sets) the propId (mvp-id), the assoc segment, and the 
  void countPropPerFrame();

  /// obsolete: can be used to assign the mvp stored in the current segment as expansion proposal
  void initSeg2Prop( int segId, std::vector<int>& seg2Prop );

  /// init function for paralel per pixel expansion/refinement
  void initElemStack(int amount)
  {
    storedElem.resize(amount);
    for (int i=0;i< storedElem.size();i++)
      storedElem[i].proposal = -1;
  }

  /// init function for parallel per segment expansion
  void initPropStack(int amount)
  {
    storedProp.resize(amount);
    for (int i=0;i< storedProp.size();i++)
      storedProp[i] = -1;
  }

  /*! multi processor: a test whether the proposal can be used
  expansion areas may not overlap if evaluated in parallel.
  In fact occlusions between areas expanded in parallel should also not exist.
  Segment version.
  */
  bool testNewProp( int proposal );

  // multi processor, pid: processor id
  const dataElem<Scalar>& getStoredElem(int pid) 
  {
    if (pid >= storedElem.size() )
      printf("not in range of storeElem: %d, %d", pid, storedElem.size());
    return storedElem[pid];
  }

  /// multi processor, pid: processor id
  void storeElem(int pid) 
  {
    if (pid >= storedElem.size() )
    {
      printf("not in range of storeElem: %d, %d", pid, storedElem.size());
      return;
    }
    storedElem[pid] = newElem;
  }

  /// multi processor, pid: processor id; if propId >=0 store for this id, else release the spot
  void storeReleaseProp(int pid, int propId = -1)
  {
    if (pid >= storedProp.size() )
    {
      printf("not in range of storeElem: %d, %d", pid, storedProp.size());
      return;
    }
    storedProp[pid] = propId;
  }

  ///////////////////////////////////////////////////////////////////////
  // PIXEL STUFF //
  ///////////////////////////////////////////////////////////////////////

  /// intializes newElem with info about new area to be expanded
  void preInitPixel( int proposal, genHomoG<Scalar>* gHom, const std::vector<P3>* normals );

    /*! normals here are stored vNormals in other view!!! these must exist - could also use a dummy version? providing vnormals nad homos already.
  copies the seg2var assignments from BOTH inputs !!
  */
  void preinitPixelFromView( int proposal, const dataElem<Scalar>& leader, const dataElem<Scalar>& second);

  /// determines the expansion order, areas far away are evaluated first
  void sortInitNewPixel( std::vector<int>& order, std::vector<int>& seg2Prop );

  /*! multi processor: a test whether the proposal can be used
  expansion areas may not overlap if evaluated in parallel.
  In fact occlusions between areas expanded in parallel should also not exist.
  Segment version.
  */
  bool testNewPixel( int centerId, int proposal, int boxX, int boxY, int expansionArea =1 );

  /*! to be called to generate expansion areas  (see initNewSeg)
  seg2seg is now pixel to pixel and seg2var is pixeltovar segId init the position of the expansion, nothing else
  */
  int initNewPixel( int centerId, int proposal, int boxX, int boxY, int expansionArea, int addToSeg2Vars = 0, int addTo2AllVars = 0  );

  /*! normals here are stored vNormals in other view!!! these must exist - could also use a dummy version? providing vnormals nad homos already.
  copies the seg2var assignments from BOTH inputs !!, centerId delivers the position! where to expand the box
  */
  void initNewPixelFromView( int centerId, int proposal, const dataElem<Scalar>& leader, const dataElem<Scalar>& second);

  void createWarpsPixel( )
   {
     if ( pix_idxFW.size() >0 || pix_idxBW.size() >0 )
       createWarpsPixel( pix_idxFW, pix_idyFW, pix_idxBW, pix_idyBW);
   }

   /// computes the photometric data cost per pixel in exaplnsion area (call initNewPixel/initNewPixelFromView first)
   void createWarpsPixel( std::vector< Scalar >& idxFW, std::vector< Scalar >& idyFW, 
                          std::vector< Scalar >& idxBW, std::vector< Scalar >& idyBW );
  ///////////////////////////////////////////////////////////////////////
  // END PIXEL STUFF //
  ///////////////////////////////////////////////////////////////////////

  /// compute the energy as a number of current solution (per segment step)
  Scalar getEnergy( ) ;

  /// visualization, delivers the identified class per pixel: occluded, oob, imp. etc.
  void getEnergyMap( std::vector<int> &energyMap, std::vector<int> &energyMap2 );

  /// visualization of the photometric data error per segment
  void getDataEnergyMap( std::vector<Scalar> &energyMap, std::vector<Scalar> &energyMap2 );

  /// visualization
  void getAutoScores( std::vector<Scalar> &energyMap, std::vector<Scalar> &energyMap2 );

  /// for visualization, data energy as an image
  void getEnergyMapPixel( std::vector<int> &energyMap, std::vector<int> &energyMap2 );

  /// return the energy from data term
  Scalar getEnergyPixel( );

  /// Debugging: compare the solutions from all 0's and the one with 1's set as suggested - where are the differences?
  Scalar getEnergyPixelLocal( std::vector<int>& solution01, int trial );

 /// debug function
  void checkSolution(const std::vector<int>& _currentSolution, const std::vector<int> &_currentSolution2 );

  void buildFromCurrentSolution( const std::vector<int>& _currentSolution, const std::vector<int> &_currentSolution2 );

  /*! Given the curent solution (local variable), the full cost is computed
      for all pixel in the image. 
      Used for expansion moves below, also consistently updated (see below)
  */
  void buildFromCurrentSolutionPixel( );
  // !!!!!!!!!!
  // currentsolution is not to lookup in the dataElements array !! 
  // currentSolution == normal/rotation; dataElement == normal/rot + location ( == expCenter) combiniation -> not even the same amount
  // !!!!!!!!!!

  std::vector<MotionPairs>& getMotionPairsFW() {return motionPairsFW;};
  std::vector<MotionPairs>& getMotionPairsBW() {return motionPairsBW;};

  /*! Delivers the data term as edges and unaries per segment, 
      with the variable names given in the expansion set.
      According to trial a region get expanded for graph cut minimization.
      Photometric data cost per pixel were computed way before, stored in dataElem.
      Here resulting binary costs are computed, together with variable ids.
  */
  void expandRegion( int propId );

  /*! Delivers the data term per pixel as edges and unaries, 
      with the variable names given in the expansion set.
      According to trial a region get expanded for graph cut minimization.
      Photometric data cost per pixel were computed before.
      Here resulting binary costs are computed, together with variable ids.
  */
  void expandRegionPixel( int trial );

  /// ceck the depth maps induced
  void expandRegionJointGeometry( const dataElem<Scalar> element, std::vector<int*> segmentations);

   /// check the depth maps induced by segment, element to get the seg2var combination, currSolE2, currSol1A for prop=0 lookup
  void expandSegRegionJointGeometry( const dataElem<Scalar> element, int propId, 
                                     std::vector<int>& currSol1A, std::vector<int>& currSol2A);

  Scalar getEnergySegRegionJointGeometry( std::vector<int>& currSol1A, std::vector<int>& currSol2A);

 ////////////////////////////////////
 // needs also dataElements[prop0] and currentSolutionE vectors
 void expandRegionJointSegmentation( std::vector <const dataElem<Scalar> > elements, 
                                     std::vector <const dataElem<Scalar> > currSolVec, 
                                     std::vector< const std::vector <dataElem<Scalar> > >  dtaElemVec, 
                                     std::vector<int*> segmentations, int& lastVar);

  /// according to the labeling from graphcut update solution internally
  void updateSolutions( int propId, std::vector<int>& _currentSolution1, std::vector<int>& _currentSolution2, const std::vector<int>& solution01 );

  void countNoSolutions(std::vector<int>& cns, std::vector<int>& counter );

  //  void updateSolutionsPixel( int trial, std::vector<int>& _currentSolution1, std::vector<int>& _currentSolution2, const std::vector<int>& solution01 )
  /// must also update the zbuffer information
  void updateSolutionsPixel( int trial, const std::vector<int>& solution01, int fromPotElem =-1 );

  /// for parallel stuff
  void invalidateElem( int fromPotElem ) 
  {     
    assert(fromPotElem  < storedElem.size());
    if (fromPotElem >= 0 && fromPotElem < storedElem.size()) storedElem[fromPotElem].proposal =-1; 
  };

  const dataElem<Scalar>& getDataElem ( int propId ) {return dataElements[propId];};
  typename std::vector< Binary<Scalar> > &getBinariesCons(){return binariesCons;};
  typename std::vector<Unary<Scalar> >   &getUnariesCons() {return unariesCons;};
  typename std::vector< Binary<Scalar> > &getBinaries(){return binaries;};
//  std::vector<P2>  &getUnaries() {return unaries;};
  typename std::vector<Unary<Scalar> >  &getUnaries() {return unaries;};
  std::vector<int> &getSegmentsInvolvedFW(int trial) {return dataElements[trial].seg2varFW;};
  std::vector<int> &getSegmentsInvolvedBW(int trial) {return dataElements[trial].seg2varBW;};
  int getNVars(int trial) {return dataElements[trial].dataFW.size() + dataElements[trial].dataBW.size(); }
  std::vector<int> &getSegmentsInvolvedFWPixel( )    {return newElem.seg2varFW;};
  std::vector<int> &getSegmentsInvolvedBWPixel(int trial) {return newElem.seg2varBW;};
  int getNVarsPixel(int trial) {return newElem.dataFW.size() + newElem.dataBW.size(); }

private:

  /*! idea already knowledge of the segment, call right there in expansion per segment, 
  * find segment in other view which might be occluded by this segment.
  * Example for all FW segments with a variable we need to check if a segment
  * in the other view not in the variable list falls onto it.
  * If so, add a unary term, if assigning label0/1 changes the score of that pixel.
  * so: e.g. FW variable varIdx is segment segId. Check all pixel in 2nd view (BW).
  * which fall onto this one (pixId). Use OmapBufBW! to get all pixel in BW view falling
  * onto it.
  */
  void expandRegionLongRangeInteractions( int propId );

  /*! idea already knowledge of the pixel, call right there in expansion pixel, 
  * find pixel in other view which might be occluded by this pixel.
  * Example for all FW pixel with a variable we need to check if a pixel 
  * in the other view not in the variable list falls onto it.
  * If so, add a nuary term, if assigning label0/1 changes the score of that pixel.
  * so: e.g. FW variable varIdx is pixel pixId. Check all pixel in 2nd view (BW).
  * which fall onto this one (pixId). Use OmapBufBW! to get all pixel in BW view falling
  * onto it.
  */
//  expandRegionLongRangeInteractionsPixel( pixId, varIdx, trial, FW );
  void expandRegionLongRangeInteractionsPixel( int pixId, int varIdx, int trial, int FW=1 );

  /// first project all centers and sort them - actually no projection needed
  void buildList ( M3& KMat, Scalar* centers, int nSegments, typename std::vector< PatchCenter<Scalar> >& projPatchCenters);

  void setKs ( Scalar* Kl_, Scalar* Kl2_  );

  void setAreas();

  inline bool isinBounds(P3 p0)
  {
    //        return (int( p0[0]-0.5 ) >= 0) && (int( p0[0]-0.5 ) < w) && (int( p0[1]-0.5 ) >= 0) && (int( p0[1]-0.5 ) < h);
    return (int( floor(p0[0]-0.5) ) >= 0) && (int( floor(p0[0]-0.5) ) < w) && (int( floor(p0[1]-0.5) ) >= 0) && (int( floor(p0[1]-0.5) ) < h);
  }

  inline int toPixId( P3 p0 )
  {
    //return int( p0[0]-0.5 )*h +int ( p0[1]-0.5 );
    assert(p0[2] == 1.);
    return int( floor (p0[0]-0.5) )*h +int ( floor(p0[1]-0.5) );
  }

  P3 P3FromGlobal(int i) {return P3(int(i/h)+1., (i%h)+1., 1.);};
//  int GlobalFromP3(P3 p) {assert(p[2] == 1.);return int( floor(p[0]-0.5) )*h +int ( floor(p[1]-0.5) );};

  /*! compute penalty for binary data term, area: full area of segment, freePixel how many of these are actually 'free'
  iNom: currently stored normal (2nd view), Nom the normal projected into the 2nd view 
  */
  Scalar getPenalty ( const M3& Hom, P3 Nom, P3 iNom, P3 centre, int area, 
                      Scalar data, Scalar freeData, Scalar freePix, 
                      int& edgeCase, Scalar autoScore);

  int whichPenalty ( M3 Hom, P3 Nom, P3 iNom, P3 centre, M3& iHom, P3 centre2, int prop1, int prop2 );

    // Actually used ?! on Brutus its the double check - which is similar to large tolerance - inmportant for per Segment not so much here
  // getPenalty_x2_asITShouldBe 
  Scalar getPenalty_x2 ( const M3& Hom, const P3& Nom, const P3& iNom, int pixId, P3 centre, 
                         int area, Scalar data, Scalar freeData, Scalar freePix, int& edgeCase, Scalar autoScore);

  // projection matrices needed for some stuff
  M3 Kl;
  M3 Kl2;
  M3 iKl;
  M3 iKlt;
  M3 iKl2;
  M3 iKlt2;
//  M3 Klt;
//  M3 Klt2;

  /// somehow absurd but this returns the original normal from Nom and iNom parts.
  M3 iiKlt2;

  /// for the trial only temporarily installed 
  std::vector< Scalar > pix_idxFW; 
  std::vector< Scalar > pix_idyFW; 
  std::vector< Scalar > pix_idxBW; 
  std::vector< Scalar > pix_idyBW;

  /// just for testing -- what if no view consistency is applied
  bool simpleApproach;

  bool impossibleApproach;

  /// turn off if not needed to save memory ?
  bool doOccBuffering;

  /// idea: these define the expansion areas, and can be totally unrelated to anything else, also being only visible in 2nd view, etc
  Scalar* expCenters;
  int nExpCenters;

  Scalar* centers;
  Scalar* centers2;
  /// the amount of segmetns in the image 1
  int nSegments;
  /// the amount of segmetns in the image 2
  int nSegments2;
  typename std::vector< PatchCenter<Scalar> > projPatchCenters;
  typename std::vector< PatchCenter<Scalar> > projPatchCenters2;

  P4i box1;
  P4i box2;

  genWarpS<Scalar> gws1, gws2;

  genScore<Scalar> gs1, gs2;

  typename std::vector< genWarpS<Scalar>* > gws1V, gws2V;

  typename std::vector< genScore<Scalar>* > gs1V, gs2V;

  std::vector<bool>& occlusions1;
  std::vector<bool>& occlusions2;

  // for single processor
  dataElem<Scalar> newElem;
  /// for parallel/multi processor
  typename std::vector< dataElem<Scalar> > newElemV;
  
  /// for parallel per pixel processing stuff
  typename std::vector < dataElem<Scalar> >  storedElem;

  /// for parallel paer segment processing stuff
  typename std::vector < int >  storedProp;  

  typename std::vector < dataElem<Scalar> >  dataElements;

  // !!!!!!!!!!
  // currentsolution is not to lookup in the dataElements array !! 
  // currentSolution == normal/rotation; dataElement == normal/rot + location ( == expCenter) combiniation -> not the same amount
  // !!!!!!!!!!

  /// stores the current solution(s)
  dataElem<Scalar> currentSolutionE;

  // mapping each segment to a proposal (to noraml/rotation); here a prop also contains location == expCenter
  std::vector<int> currentSolution1;
  std::vector<int> currentSolution2;

  /// stores # of pixel in the segments
  std::vector<int>  areas1;
  std::vector<int>  areas2;

  typename std::vector< Binary<Scalar> > binariesSegi;

  typename std::vector< Binary<Scalar> > binariesCons;
  typename std::vector< Unary<Scalar>  > unariesCons;

  /// merge in advance??? not possible
  typename std::vector< Binary<Scalar> > binaries;
  typename std::vector< Unary<Scalar>  > unaries;
  /// unaries : for each local variable one
  //  std::vector<P2> unaries;

  std::vector<P3>  Centers;// patch centers in view1
  std::vector<P3> iCenters;// patch centers in view2

  /// width
  int w;
  /// height
  int h;

  /// the segment ids for later test w.r.t. occlusions
  int* segImg;
  /// the segment2 ids for later test w.r.t. occlusions
  int* segImg2;

  /// makes sense to precompute: for each proposal one homo/and inverse one
  std::vector<M3> homs_pix;  // homographies in pixel coordinates - in contrast to pixel coords
  std::vector<M3> ihoms_pix; // homographies in pixel coordinates
  std::vector<P3> noms_pix;  // view normals in pixel coordinates - in contrast to pixel coords
  std::vector<P3> inoms_pix; // view normals in pixel coordinates

  std::vector<P3> vNoms_cam;// view normals, in camera coordiantes to compute the depth at that pixel quickly

  /// maximal size of a box in a projected view
  const Scalar resizer;

  Scalar occPenalty;
  Scalar impPenalty;

  Scalar impFactor;

  Scalar __OCC_THRESH__;
  /// penalty if mapping is between different planes
  Scalar addPotts;

  /// max distance allowed for mathcing moving planes
//  Scalar mvp_dst;

  /// maximum depth distance for occlusion checking. if below we say patches match anyway
  Scalar maxDepthThresh;

  /// for checking the normal differnece of matching elements
//  Scalar maxNormDiff;

  /// the difference in 3d which is allowed: e.g. depth diff: small, motion in coupling: 'larger'
  Scalar scale3D;
  /// buffer for per pixel expansion handling:
  OcclusionExpansionBuffer<Scalar> oMapBufFW;
  OcclusionExpansionBuffer<Scalar> oMapBufBW;

  /// stores the pairs of segments/assignments for temporal smoothness
  std::vector<MotionPairs> motionPairsFW;
  std::vector<MotionPairs> motionPairsBW;

  /// makes sense for motion image pairs only -- not in this version
  bool computeMotionPairs;

  Scalar* autoScoresFW;
  Scalar* autoScoresBW;

  std::vector<Scalar> adaptiveAutoScoresFW;
  std::vector<Scalar> adaptiveAutoScoresBW;

  /// info about view pair, camera id and time step are stored
  P4i ctct;
  Scalar xtraPen;
  Scalar gridSize;

  // maybe exit more early leads to less computation ?? 
  Scalar maxMotSqr;
};

#if defined(INCLUDE_TEMPLATES) && !defined(__DataPerSegment__cpp)
#include "DataPerSegment.cpp"
#endif

#if defined(INCLUDE_TEMPLATES) && !defined(__DataPerPixel__cpp)
#include "DataPerPixel.cpp"
#endif

#if defined(INCLUDE_TEMPLATES) && !defined(__DataGeneral__cpp)
#include "DataSegmentGeneral.cpp"
#endif

#endif // __OcclusionMapping__h
