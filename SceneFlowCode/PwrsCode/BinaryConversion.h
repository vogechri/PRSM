/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
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
License along with this software.*/

////////////////////////////////////////////////
// Convert higher order terms to binary terms //
////////////////////////////////////////////////

#ifndef __BinaryConversion__h
#define __BinaryConversion__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <set>
#include <limits>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include "OcclusionMappingBufferSegments.h" 
#include "BinaryConversion_HeapHelper.h" 
#include "Templates/HeapT.h"

using namespace std;
using namespace Math;

template<typename Scalar>
struct Unaries
{
  Unaries( )
    : label(-1), segId(-1), score(0)
  {};

  Unaries( int _label, int _segId, Scalar _score )
    : label(_label), segId(_segId), score(_score)
  {};
  //
  int label;
  int segId;
  Scalar score;
};

template<typename Scalar>
struct higherOrders
{
  typedef std::pair <int,int>  SegLabelPair;
  higherOrders( ) : score(0)
  {vars.reserve(5);};

  higherOrders( Scalar _score ) : score (_score)
  {vars.reserve(5);};

  ~higherOrders() {};
  ////////////////////////////////////////////

  bool operator==(const higherOrders<Scalar>& _rhs) const 
  {
    if ( vars.size() != _rhs.vars.size() )
      return false;
    for (int i =0; i<vars.size(); i++)
      if (vars[i].first != _rhs.vars[i].first || vars[i].second != _rhs.vars[i].second)
        return false;

    return true;
  }

  bool isSubMod()
  {
    assert(vars.size() > 1);
    if ( vars.size() == 2 && score>=0 && vars[0].second != vars[1].second )
      return true;
    if ( vars.size() == 2 && score<=0 && vars[0].second == vars[1].second )
      return true;
    if ( vars.size() > 2  && score < 0)
    {
      for(int k =1;k<vars.size();k++)
        if (vars[0].second != vars[k].second )
          return false;
      return true;
    }

    return false;
  }

  /// stored are pairs: variable, assignment
  std::vector < SegLabelPair > vars;

  Scalar score;
};

template<typename Scalar>
struct literalTerm
{
  literalTerm() : score(0) {};
  literalTerm(Scalar _score) : score(_score) {};
  ~literalTerm(){};

  std::set < int > literals;
  Scalar score;
};

// idea:
// build list of positive terms and negative terms
// contructed by ausmultiplizieren, and gathering
// the lists here are used later after combining all n (n view pairs) lists
// replacing the literal with the highest amount of occurences
//
// also gather parts which are sub-modular:
// all positive (neagtive) literals and negative weight
//
// replacement: goal only positive literals in term:
//
// replace (1-x_i)*ti by  ti - x_i*ti
//
// repeat until all negative and value negative or all positive
// value negative: on sub-mod list
// value positive: on replacement list
// 
// replacement list from all terms: gather, apply reduction
// sub-mod list: turn to sub mod binary terms with old ideas
//
// 0 step: 
// gather scores: occlusion-data and occlusion term as unary, note that some views are more tricky:
// merging as unifying both occlusion variables
// 1st step: 
// replace double occurences, mark already submodular parts
// 2.nd step:
// 
///////////////////////////////////////
//use :
// BinaryConverter<Scalar> bc_0( _occPenalty, std::vector< occSegList<Scalar> >& _occList, std::vector<int>& _freePix0, std::vector<int>& _freePix1  );
// BinaryConverter<Scalar> bc_1( _occPenalty, std::vector< occSegList<Scalar> >& _occList, std::vector< occSegList<Scalar> >& _occList2, std::vector<int>& _freePix0, std::vector<int>& _freePix1  )
// bc_0.setPenaltiesPerSegment( std::vector<Scalar>& _fixedPenalties, std::vector<Scalar>& _dataPen0, std::vector<Scalar>& _dataPen1 )
// bc_0.generateHigherOrderTerms();

//HigherOrderConverter<Scalar> hoc(nSegments);
// hoc.clear();
// hoc.addNonSubs ( bc_0.get_nSubmods() );
// hoc.addSubs ( bc_0.get_Submods(), bc_0.getUnaries() );
// hoc.convert()

template<typename Scalar> class BinaryConverter
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;

  typedef std::pair <int,int>       Seg_ij;
  typedef std::pair <int,int>       SegLbl;

  //  typedef std::pair <int,Scalar>    Unaries;

  // std::vector<SegLbl>& _fullOcclusions
  BinaryConverter( Scalar _occPenalty, std::vector< occSegList<Scalar> >& _occList, 
    std::vector<int>& _freePix0, std::vector<int>& _freePix1 )
    : occPenalty(_occPenalty), occlusionList(_occList) // , fullOcclusions(_fullOcclusions)
  {
    assert ( _freePix0.size() == _freePix1.size() );
    nSegments = _freePix0.size();
    setFreePixPerSegment( _freePix0, _freePix1 );
  }

  BinaryConverter( Scalar _occPenalty, std::vector< occSegList<Scalar> >& _occList, std::vector< occSegList<Scalar> >& _occList2, 
    std::vector<int>& _freePix0, std::vector<int>& _freePix1 )
    : occPenalty(_occPenalty)
  {
    assert ( _freePix0.size() == _freePix1.size() );
    nSegments = _freePix0.size();
    setFreePixPerSegment( _freePix0, _freePix1 );
    mergeLists( _occList, _occList2, occlusionList );
  }

  // special case per pixel:
  BinaryConverter( Scalar _occPenalty, 
    std::list< occPixelList< Scalar > >& _occList0, 
    std::list< occPixelList< Scalar > >& _occList1, 
    std::vector<Scalar>& _penPix) //, std::vector<Scalar>& _penPix1 )
    : occPenalty(_occPenalty)
  {
    // convert into this:
    // std::vector< occSegList<Scalar> > occList
//    assert ( _penPix.size() );
    nSegments = _penPix.size()/4;
    std::vector<int> _freePix0(nSegments,0);
    std::vector<int> _freePix1(nSegments,0); 

    for (int i = 0;i<nSegments;i++)
    {
       if ( _penPix[i +   nSegments]<=Scalar(0) )
         _freePix1[i] = 1;
       if ( _penPix[i + 3*nSegments]<=Scalar(0) )
         _freePix0[i] = 1;
    }
    setFreePixPerSegment( _freePix0, _freePix1 );

    // fill this:
    occlusionList.reserve( _occList0.size() + _occList1.size() );
    fillFromPixelList ( occlusionList, _occList0, 0 );
    fillFromPixelList ( occlusionList, _occList1, 1 );
    std::sort( occlusionList.begin(), occlusionList.end() );
  }

    // special case per pixel:
  BinaryConverter( Scalar _occPenalty, 
    std::list< occPixelList< Scalar > >& _occList0, 
    std::list< occPixelList< Scalar > >& _occList1, 
    std::list< occPixelList< Scalar > >& _occList0B, 
    std::list< occPixelList< Scalar > >& _occList1B, 
    std::vector<Scalar>& _penPix)//0, std::vector<Scalar>& _penPix1 )
    : occPenalty(_occPenalty)
  {
    // convert into this:
    // std::vector< occSegList<Scalar> > occList
    nSegments = _penPix.size()/4;
    std::vector<int> _freePix0(nSegments,0);
    std::vector<int> _freePix1(nSegments,0); 

    for (int i = 0;i<nSegments;i++)
    {
       if ( _penPix[i +   nSegments]<=Scalar(0) )
         _freePix1[i] = 1;
       if ( _penPix[i + 3*nSegments]<=Scalar(0) )
         _freePix0[i] = 1;
    }
    setFreePixPerSegment( _freePix0, _freePix1 );

    // fill this:
    std::vector< occSegList<Scalar> > occlusionList1;
    occlusionList1.reserve( _occList0.size() + _occList1.size() );
    std::vector< occSegList<Scalar> > occlusionList2;
    occlusionList2.reserve( _occList0B.size() + _occList1B.size() );

    fillFromPixelList ( occlusionList1, _occList0, 0 );
    fillFromPixelList ( occlusionList1, _occList1, 1 );
    std::sort( occlusionList1.begin(), occlusionList1.end() );
    fillFromPixelList ( occlusionList2, _occList0, 0 );
    fillFromPixelList ( occlusionList2, _occList1, 1 );
    std::sort( occlusionList2.begin(), occlusionList2.end() );

    occlusionList.reserve( occlusionList2.size() + occlusionList1.size() );
    mergeLists( occlusionList1, occlusionList2, occlusionList );
  }


  ~BinaryConverter(){};

  void addOcclusionListOuterPixels( std::list< occPixelList< Scalar > >& _occListOuter );

  void addOcclusionListOuterPixels( std::list< occPixelList< Scalar > >& _occListOuter, 
     std::list< occPixelList< Scalar > >& _occListOuter2);

  std::vector< Unaries<Scalar> >& getUnaries()        {return unaries;};
  std::vector< higherOrders<Scalar> >& get_Submods()  {return submods;};
  std::vector< higherOrders<Scalar> >& get_nSubmods() {return nonsubmods;};

  /// to be called before merging - the problem: free pixel for labels0 and labels 1 needed
  void setFreePixPerSegment( std::vector<int>& _freePix0, std::vector<int>& _freePix1 )
  {
    freePix.resize( _freePix0.size() + _freePix1.size() ); 
    std::copy(_freePix0.begin(), _freePix0.end(), freePix.begin() );
    std::copy(_freePix1.begin(), _freePix1.end(), freePix.begin() + _freePix0.size() );
  }

  std::vector<SegLbl>& getFullOcc() {return fullOcclusions;}
  /// list of all occlusions observed for the view
  std::vector< occSegList<Scalar> >& getOccList() {return  occlusionList;};

  /// those penalties are unaffected from occlusions, since these means unresonable behaviour
  void setPenaltiesPerSegment( std::vector<Scalar>& _dataPen0, std::vector<Scalar>& _dataPen1 ) //std::vector<Scalar>& _fixedPenalties, 
  { 
//    fixedPenalties = _fixedPenalties;
    dataPen0       = _dataPen0;
    dataPen1       = _dataPen1;
  }

  /// those penalties are unaffected from occlusions, since these means unresonable behaviour, special case per pixel:
  void setPenaltiesPerSegment( std::vector<Scalar>& _dataPen)//, std::vector<Scalar>& _dataPen1 )
  { 
    //    fixedPenalties = dataPen0; // anything
    assert(nSegments = _dataPen.size()/4);
    dataPen0.resize(nSegments);
    dataPen1.resize(nSegments);

    std::copy(_dataPen.begin(), _dataPen.begin()+nSegments, dataPen0.begin() );
    std::copy(_dataPen.begin()+2*nSegments, _dataPen.begin()+3*nSegments, dataPen1.begin() );
  }

  /// build unaries for drawing - only label 0 if occluded by label 0
  Scalar generateUnaries0L0()
  {
    Scalar fullPenalty (0);
    unaries.clear();
    for (int i = 0; i < occlusionList.size(); i++)
    {
      int sid = occlusionList[i].segId;

      int label = occlusionList[i].label;
      Scalar penalty;

      if ( !occlusionList[i].isLabel0OccludedbyLabel0() ) continue;

      penalty = occPenalty * freePix[sid + label*nSegments] - dataPen0[sid];
      fullPenalty += penalty;
      // now unary:
      unaries.push_back (Unaries<Scalar> (label, sid, penalty) );
    }
    return fullPenalty;
  }

  /// the same for the outer occlusion list !
  void generateHigherOrderTerms_outer( std::vector<Scalar>& globalDataScores );

  /// generate the unaries to be added, split up into submodular elements and non-submods
  void generateHigherOrderTerms();

  /// free/delete the store
  void freeLists() {};
  void addList() {};

private:

  void fillFromPixelList ( std::vector< occSegList<Scalar> >&occlusionList, std::list< occPixelList< Scalar > >& _occList0, int ownlbl, int startId=1 );

  /// merge 2 lists as two views combine, MUST check if no free Vars: ignore !
  void mergeLists ( std::vector< occSegList<Scalar> >& oListsa, std::vector< occSegList<Scalar> >& oListsb, 
                    std::vector< occSegList<Scalar> >& mergedList );

  /// merge 2 lists as two views combine, MUST check if no free Vars: ignore !
  void mergeListsOuter ( std::vector< occSegList<Scalar> >& oListsa, std::vector< occSegList<Scalar> >& oListsb, 
                         std::vector< occSegList<Scalar> >& mergedList );

  ///////////////////////////////////////////////

  /// the amount of segmetns in the images
  int nSegments;

  Scalar occPenalty;

  /// number of free pixel per segment, first part 0 label, 2nd part label 1
  std::vector<int> freePix;

  /// data penalties - to be replaced by # free pixel times the occlusion penalty
  std::vector<Scalar> dataPen0;
  std::vector<Scalar> dataPen1;

  /// vector of all segments which are anyway occluded
  std::vector<SegLbl> fullOcclusions;
  /// list of all occlusions observed for the view
  std::vector< occSegList<Scalar> > occlusionList;

  /// in the per pixel case the outer list stores occlusions for pixel not under evaluation
  std::vector< occSegList<Scalar> > occlusionListOuter;

  /// the altered unaries only:
  std::vector<Unaries<Scalar> > unaries;

  /// 
  std::vector< higherOrders<Scalar> > submods;

  /// the list will be created from this one
  std::vector< higherOrders<Scalar> > nonsubmods;

};


// use:
//HigherOrderConverter<Scalar> hoc(nSegments);
// clear();
// addSubs ( std::vector< higherOrders<Scalar> > _outer_submods, std::vector< Unaries<Scalar> > _unaries )
// convert()
// a set for the literals, recursive procedure for converting to positive only literals
template<typename ScalarT> class HigherOrderConverter
{
public:

  typedef ScalarT Scalar;
  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;

  typedef std::pair <int,int>       Seg_ij;
  typedef std::pair <int,int>       SegLbl;

  HigherOrderConverter(int _nSegments) : maxRecursionDepth(10) {setNSegments(_nSegments);};

  ~HigherOrderConverter(){};

  void setNSegments(int _nSegments) {nSegments = _nSegments;nVars = nSegments;};

  void clear() {
    outer_unaries.clear(); outer_submods.clear();nVars = nSegments;
    bl_00.clear();bl_01.clear();
    bl_10.clear();bl_11.clear();
    b0.clear();b1.clear();
    subMods.clear();
    nsubMods.clear();//positions.clear();terms.clear();
  };

  int get_nVars()  {return nVars;};
  int get_nEdges() {return bl_00.size() + bl_01.size() + bl_10.size() + bl_11.size();};

  std::vector< edgePenalty<Scalar> >& get_b00() {return bl_00;}; 
  std::vector< edgePenalty<Scalar> >& get_b01() {return bl_01;};
  std::vector< edgePenalty<Scalar> >& get_b10() {return bl_10;};
  std::vector< edgePenalty<Scalar> >& get_b11() {return bl_11;};

  std::vector< std::pair<int, Scalar> >& get_b0() {return b0;};
  std::vector< std::pair<int, Scalar> >& get_b1() {return b1;};

  void convert()
  {
//    printf("# nsm:%d ", nsubMods.size() );mexEvalString("drawnow"); 
    convertNonSubs();
//    printf("# osm:%d, sm:%d ", outer_submods.size(), subMods.size() );mexEvalString("drawnow"); 
    convertSubs();
  }

  void testNonSubs();

  /*! in: gathered non-submodular terms from any view pair, 
  * all build a large list of non-submodular subterms. 
  * And also submodular terms resultig from the conversion.
  */
  void addNonSubs ( std::vector< higherOrders<Scalar> > &nonsubmods );

  /// here the submodular parts are added as well:
  void addSubs ( std::vector< higherOrders<Scalar> >& _outer_submods, std::vector< Unaries<Scalar> >& _unaries );

  void addMergeNonSubs_Subs ( std::vector< higherOrders<Scalar> >& nonsubmods1, std::vector< higherOrders<Scalar> >& nonsubmods2,
                              std::vector< higherOrders<Scalar> >& submods1, std::vector< higherOrders<Scalar> >& submods2);

  void addUnaries ( std::vector< Unaries<Scalar> >& _unaries );

  private:

  // this is done term by term by the known substitution 2.3.2
  void convertSubs();

  void convertNonSubs();

  /// nsm already contains a vector of seg,label pairs and a score now convert those to positive literal terms
  void convertToTerms ( higherOrders<Scalar> nsm, std::vector<literalTerm<Scalar> > &subMods, std::vector<literalTerm<Scalar> > &nsubMods );

  void recursiveConversion ( literalTerm<Scalar> lTerm, higherOrders<Scalar>& nsm, int position, std::vector<literalTerm<Scalar> > &vec_lit );


    void getSmallerIterator( 
    typename std::vector< higherOrders<Scalar> >::iterator& nsm_it, 
    typename std::vector< higherOrders<Scalar> >::iterator& sm_it, 
    typename std::vector< higherOrders<Scalar> >::iterator& nsm_end, 
    typename std::vector< higherOrders<Scalar> >::iterator& sm_end, 
    typename std::vector< higherOrders<Scalar> >::iterator*& it, 
    typename std::vector< higherOrders<Scalar> >::iterator*& end);

  /// amount of variables without extra added
  int nSegments;
  int nVars;

  int maxRecursionDepth;

  /// gathered only positive literal terms with negative scalar multiplier: submodular (or just unaries)
  std::vector<literalTerm<Scalar> >  subMods;
  /// gathered only positive literal terms with positive scalar multiplier: non-submodular
  std::vector<literalTerm<Scalar> > nsubMods;

  /// all positive literals only
//  std::vector<Scalar>               unaries;

  std::vector< Unaries<Scalar> >      outer_unaries;
  std::vector< higherOrders<Scalar> > outer_submods;

  std::vector< edgePenalty<Scalar> > bl_00; // stores energy f00 = scalar
  std::vector< edgePenalty<Scalar> > bl_01; // stores energy f01 = scalar
  std::vector< edgePenalty<Scalar> > bl_10; // stores energy f10 = scalar
  std::vector< edgePenalty<Scalar> > bl_11; // stores energy f11 = scalar

  std::vector< std::pair<int, Scalar> > b1; // stores energy f1
  std::vector< std::pair<int, Scalar> > b0; // stores energy f0
};

#if defined(INCLUDE_TEMPLATES) && !defined(__BinaryConversion_CPP_)
#include "BinaryConversion.cpp"
#endif

#if defined(INCLUDE_TEMPLATES) && !defined(__HigherOrderConverter__cpp)
#include "HigherOrderConverter.cpp"
#endif

#endif // __OcclusionMapping__h
