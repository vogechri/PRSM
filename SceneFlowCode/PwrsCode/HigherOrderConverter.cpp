/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

////////////////////////////////////////////////
// Convert higher order terms to binary terms //
////////////////////////////////////////////////

#ifndef __HigherOrderConverter__cpp
#define __HigherOrderConverter__cpp

#include "BinaryConversion.h"

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

template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  testNonSubs()
{
  std::vector< higherOrders<Scalar> > my_nonsubmods;
  higherOrders<Scalar> nonSub1(88);
  higherOrders<Scalar> nonSub2(88);

  nonSub1.vars.push_back( higherOrders<Scalar>::SegLabelPair (1811,0) );
  nonSub1.vars.push_back( higherOrders<Scalar>::SegLabelPair (1737,0) );
  nonSub1.vars.push_back( higherOrders<Scalar>::SegLabelPair (1757,1) );
  nonSub1.vars.push_back( higherOrders<Scalar>::SegLabelPair (1782,0) );
  nonSub1.vars.push_back( higherOrders<Scalar>::SegLabelPair (1797,1) );

  nonSub2.vars.push_back( higherOrders<Scalar>::SegLabelPair (1811,1) );
  nonSub2.vars.push_back( higherOrders<Scalar>::SegLabelPair (1737,0) );
  nonSub2.vars.push_back( higherOrders<Scalar>::SegLabelPair (1757,1) );
  nonSub2.vars.push_back( higherOrders<Scalar>::SegLabelPair (1782,0) );
  nonSub2.vars.push_back( higherOrders<Scalar>::SegLabelPair (1797,1) );

  my_nonsubmods.push_back( nonSub1  ); 
  my_nonsubmods.push_back( nonSub2  ); 
  addNonSubs ( my_nonsubmods );
  convert();
}

/*! in: gathered non-submodular terms from any view pair, 
* all build a large list of non-submodular subterms. 
* And also submodular terms resultig from the conversion.
*/
template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  addNonSubs ( std::vector< higherOrders<Scalar> > &nonsubmods )
{
  // append to terms
  // update list so far
  for (int i =0; i< nonsubmods.size(); i++)
  {
    higherOrders<Scalar> nsm = nonsubmods[i];
    convertToTerms ( nsm, subMods, nsubMods );
    // and ?
  }
}

/// here the submodular parts are added as well:
template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  addSubs ( std::vector< higherOrders<Scalar> >& _outer_submods, std::vector< Unaries<Scalar> >& _unaries )
{
  // append
  if ( _outer_submods.size() >0 )
  {
    int end = outer_submods.size();
    if ( end >0 )
    {
      outer_submods.resize(end + _outer_submods.size());
      std::copy(_outer_submods.begin(), _outer_submods.end(), &(outer_submods[end]) );
    }
    else
      outer_submods = _outer_submods;
  }
  // append
  if (_unaries.size() > 0)
  {
    int end = outer_unaries.size();
    if (end >0)
    {
      outer_unaries.resize(end + _unaries.size());
      std::copy(_unaries.begin(), _unaries.end(), &(outer_unaries[end]) );
    }
    else
      outer_unaries = _unaries;
  }
}


template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  addMergeNonSubs_Subs ( std::vector< higherOrders<Scalar> >& nonsubmods1, std::vector< higherOrders<Scalar> >& nonsubmods2,
  std::vector< higherOrders<Scalar> >& submods1, std::vector< higherOrders<Scalar> >& submods2)
{
  // append to terms
  // update list so far
  typename std::vector< higherOrders<Scalar> >::iterator nsm_it1(nonsubmods1.begin()), nsm_it2(nonsubmods2.begin()), sm_it1(submods1.begin()), sm_it2(submods2.begin());
  typename std::vector< higherOrders<Scalar> >::iterator nsm_end1(nonsubmods1.end()), nsm_end2(nonsubmods2.end()), sm_end1(submods1.end()), sm_end2(submods2.end());

  //    printf("#b nsm:%d ", nsubMods.size() );mexEvalString("drawnow"); 
  //    printf("#b osm:%d, sm:%d ", outer_submods.size(), subMods.size() );mexEvalString("drawnow"); 

  while( nsm_it1 != nsm_end1 || nsm_it2 != nsm_end2 || sm_it1 != sm_end1 || sm_it2 != sm_end2 )
  {
    if ( nsm_it1 != nsm_end1 && nsm_it1->vars.size()<2) //!=2 ) 
    {
      convertToTerms ( *nsm_it1, subMods, nsubMods );
      nsm_it1++;
      continue;
    }
    if ( nsm_it2 != nsm_end2 && nsm_it2->vars.size()<2)//!=2 ) 
    {
      convertToTerms ( *nsm_it2, subMods, nsubMods );
      nsm_it2++;
      continue;
    }
    if ( sm_it1 != sm_end1 && sm_it1->vars.size()<2)//!=2 ) 
    {
      outer_submods.push_back ( *sm_it1 );
      sm_it1++;
      continue;
    }
    if ( sm_it2 != sm_end2 && sm_it2->vars.size()<2)//!=2 ) 
    {
      outer_submods.push_back ( *sm_it2 );
      sm_it2++;
      continue;
    }
    ///////////////////
    // now proceed with the smaller ones:
    typename std::vector< higherOrders<Scalar> >::iterator *it1(NULL), *itend1(NULL), *it2(NULL), *itend2(NULL); 
    getSmallerIterator( nsm_it1, sm_it1, nsm_end1, sm_end1, it1, itend1 );
    getSmallerIterator( nsm_it2, sm_it2, nsm_end2, sm_end2, it2, itend2 );

    if ( *it1 == *itend1 ) 
    {
      if ( (*it2)->isSubMod() )
        outer_submods.push_back ( *(*it2) ) ;
      else
        convertToTerms ( *(*it2), subMods, nsubMods );
      (*it2)++;
      continue;
    }
    if ( (*it2) == (*itend2) ) 
    {
      if ( (*it1)->isSubMod() )
        outer_submods.push_back ( *(*it1) ) ;
      else
        convertToTerms ( *(*it1), subMods, nsubMods );
      (*it1)++;
      continue;
    }
    // both are not sm:
    if (*(*it1) == *(*it2))
    {
      higherOrders<Scalar> nsm = *(*it1);
      nsm.score += (*it2)->score;
      if ( nsm.isSubMod() )
        outer_submods.push_back ( nsm ) ;
      else
        convertToTerms ( nsm, subMods, nsubMods );

      (*it1)++;(*it2)++;
      continue;
    }

    // do not care / both are not at the end:
    typename std::vector< higherOrders<Scalar> >::iterator *it(NULL), *itend(NULL);
    getSmallerIterator( *it1, *it2, *itend1, *itend2, it, itend );

    if ( (*it)->isSubMod() )
      outer_submods.push_back ( *(*it) ) ;
    else
      convertToTerms ( *(*it), subMods, nsubMods );
    (*it)++;

    //      if isSubMod()
    //      outer_submods.push_back ( .. ) ;
    //      else
    //         convertToTerms ( nsm, subMods, nsubMods );

    // lone nsm:
    //higherOrders<Scalar> nsm = nonsubmods[i];
    //convertToTerms ( nsm, subMods, nsubMods );
    //outer_submods.push_back ( .. ) ;
  }
  if (nsubMods.size() + outer_submods.size() + subMods.size()  > 200000 )
  {
    printf("#a nsm:%d ", nsubMods.size() );//mexEvalString("drawnow");
    printf("#a osm:%d, sm:%d ", outer_submods.size(), subMods.size() );//mexEvalString("drawnow");
  }
}

template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  addUnaries ( std::vector< Unaries<Scalar> >& _unaries )
{
  // append
  if (_unaries.size() > 0)
  {
    int end = outer_unaries.size();
    if (end >0)
    {
      outer_unaries.resize(end + _unaries.size());
      std::copy(_unaries.begin(), _unaries.end(), &(outer_unaries[end]) );
    }
    else
      outer_unaries = _unaries;
  }
#ifdef _DEBUG
  for( int i=0; i < outer_unaries.size(); i++ )
  {
    assert(outer_unaries[i].segId < nVars );
  }
#endif
}

// this is done term by term by the known substitution 2.3.2
template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  convertSubs()
{
  //-------------------------------------
  for (int i =0; i< subMods.size(); i++ )
  {
    std::set<int> lits = subMods[i].literals;
    Scalar score = subMods[i].score;

    if (fabs(score) < 0.0001) continue;

    if (lits.size() == 1)
      b1.push_back(std::pair<int, Scalar>( *(lits.begin()), score) );

    if (lits.size() == 2)
    {
      std::set<int>::const_iterator s_it = lits.begin();
      int v1 = *s_it;s_it++;
      int v2 = *s_it;
      bl_11.push_back(edgePenalty<Scalar>( v1, v2, score ) );
      assert (v1<nVars);
      assert (v2<nVars);
    }

    if (lits.size() > 2)
    {
      std::set<int>::const_iterator s_it( lits.begin() ), s_end ( lits.end() );
      int newVar = nVars++;
      b1.push_back( std::pair<int, Scalar>( newVar, -score * (lits.size()-1) ) );

      for ( ; s_it != s_end; s_it++ )
      {
        bl_11.push_back(edgePenalty<Scalar>( *s_it, newVar, score ) );
        assert (*s_it<nVars);
        assert (newVar<nVars);
      }
    }
  }
  // --------------------
  for( int i=0; i < outer_unaries.size(); i++ )
  {
    assert(outer_unaries[i].segId < nVars );
    if ( outer_unaries[i].label == 1 && outer_unaries[i].score != Scalar(0) )
      b1.push_back( std::pair<int, Scalar>( outer_unaries[i].segId, outer_unaries[i].score ) );
    if ( outer_unaries[i].label == 0 && outer_unaries[i].score != Scalar(0) )
      b0.push_back( std::pair<int, Scalar>( outer_unaries[i].segId, outer_unaries[i].score ) );
  }
  //----------------------     
  for( int i=0; i < outer_submods.size(); i++ )
  {
    std::vector < typename higherOrders<Scalar>::SegLabelPair >& vars = outer_submods[i].vars;
    Scalar score = outer_submods[i].score;

    if (fabs(score) < 0.0001) continue;

    if (vars.size() == 1)
    {
      assert( vars[0].first < nVars );
      if (vars[0].second == 1 )
        b1.push_back(std::pair<int, Scalar>( vars[0].first, score) );
      else
        b0.push_back(std::pair<int, Scalar>( vars[0].first, score) );
    }

    if (vars.size() == 2)
    {
      if (vars[0].second == 1 && vars[1].second == 1 )
        bl_11.push_back(edgePenalty<Scalar>( vars[0].first, vars[1].first, score ) );
      if (vars[0].second == 0 && vars[1].second == 1 )
        bl_01.push_back(edgePenalty<Scalar>( vars[0].first, vars[1].first, score ) );
      if (vars[0].second == 1 && vars[1].second == 0 )
        bl_10.push_back(edgePenalty<Scalar>( vars[0].first, vars[1].first, score ) );
      if (vars[0].second == 0 && vars[1].second == 0 )
        bl_00.push_back(edgePenalty<Scalar>( vars[0].first, vars[1].first, score ) );

      assert (vars[0].first<nVars);
      assert (vars[1].first<nVars);
    }

    // 2 cases only: all 1's all 0's
    if (vars.size() > 2)
    {         
      if ( vars[0].second == 1 ) // familiar case
      {
        int newVar = nVars++;
        b1.push_back( std::pair<int, Scalar>( newVar, -score * (vars.size()-1) ) );

        for ( int k =0; k < vars.size(); k++ )
        {
          bl_11.push_back(edgePenalty<Scalar>( vars[k].first, newVar, score ) );
          assert (newVar<nVars);
          assert (vars[k].first<nVars);
        }
      }

      if ( vars[0].second == 0 ) // other case
      {
        int newVar = nVars++;
        b0.push_back( std::pair<int, Scalar>( newVar, -score * (vars.size()-1) ) );

        for ( int k =0; k < vars.size(); k++ )
          bl_00.push_back(edgePenalty<Scalar>( vars[k].first, newVar, score ) );         
      }
    }
  }
}

template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  convertNonSubs()
{

  // 1st gather statistics: which variables to substitute first 

  std::vector< std::set<int> > lit2term(nSegments);

  for (int i =0; i< nsubMods.size(); i++ )
  {
    if (nsubMods[i].literals.size() == 1)
    {
      if ( nsubMods[i].score != Scalar(0.0) )
        b1.push_back( std::pair<int, Scalar>( *(nsubMods[i].literals.begin()), nsubMods[i].score) );
      continue;
    }

    std::set<int>::const_iterator s_it (nsubMods[i].literals.begin()), s_end (nsubMods[i].literals.end());
    for (;s_it!=s_end;s_it++)
      //      for (int k = 0; k < nsubMods[i].literals.size(); k++ )
      lit2term[ *s_it ].insert( i );
  }
  //

  typedef HeapStoreEntryTb<Scalar>       HeapStoreEntry;
  typedef std::vector< HeapStoreEntry >  HeapEntryStore;
  typedef HeapInterfaceDb < HeapEntryB, HeapEntryStore > HeapInterface;

  /// stores relevant information about HE collapse
  HeapEntryStore _heapEntryStore;

  /// the heap-Interface structure
  HeapInterface _heapInterface(&_heapEntryStore);

  /// the heap for selecting the element with the least cost
  Utils::HeapT< HeapEntryB, HeapInterface > _heap(_heapInterface);

  _heapEntryStore.resize( nSegments, HeapStoreEntry(0, -1) ); 

  // init heap:
  for (int k=0;k<lit2term.size(); k ++)
  {
    if (lit2term[ k ].size() > 0)// literal k in more than 0 terms:
    {
      HeapStoreEntry &hse = _heapEntryStore[k];
      hse.setCost( lit2term[ k ].size() );
      hse.setLiteralNr( k );

      _heap.insert( HeapEntryB (k) );
    }
  }

  while ( !_heap.empty() )
  {
    HeapEntryB id = _heap.front();
    _heap.pop_front();
    int cSeg = id.pId();
    HeapStoreEntry& heapE = _heapEntryStore[ cSeg ];

    // process all terms the segment is present
    if ( lit2term[ cSeg ].size() == 0 ) // not in one term
      break; // finished

    // ?? : would make more sense to replace the terms by itself now instead of this approach - 
    //      but it delivers the same amount of variables since we have only positive literals !
    //if ( lit2term[ cSeg ].size() == 1 ) // in a single term
    //{
    //  continue;
    //}

    // special case: no new variables: if current exists only once in a binary, then the 2nd binary also only in this edge:
    if ( lit2term[ cSeg ].size() == 1  && nsubMods[ *((lit2term[ cSeg ]).begin()) ].literals.size() == 2)
    {
      int pos = *((lit2term[ cSeg ]).begin());
      // one of them is cseg?!
      assert( *(nsubMods[ pos ].literals.find(cSeg)) == cSeg );

      nsubMods[ *(lit2term[ cSeg ].begin()) ].literals.erase(cSeg);
      int var = *(nsubMods[ pos ].literals.begin());

      // binary non-submodular edge: var1, var2, score
      Scalar score = nsubMods[ pos ].score;

      // the non-sub edge:
      bl_11.push_back( edgePenalty<Scalar>(var, cSeg, score ) );
      assert (var<nVars);
      assert (cSeg<nVars);

      // make elem var unimportant
      HeapStoreEntry& hse = _heapEntryStore[ var ];
      hse.setCost(0);
      _heap.update( HeapEntryB( var ) );

      continue;
    }
    //

    // usual case, process uniformly for all variables
    std::set<int>::const_iterator s_it( lit2term[ cSeg ].begin() ), s_end (lit2term[ cSeg ].end() ); 
    // step one: sum the values/scores:
    Scalar sumScore(0);
    for (;s_it != s_end;s_it++)
      sumScore += nsubMods[*s_it].score;

    // introduce new variable:
    int newVar = nVars++;

    // the non-submodluar edge:
    bl_11.push_back( edgePenalty<Scalar>(newVar, cSeg, sumScore ) );
    assert (newVar<nVars);
    assert (cSeg<nVars);

    // submodular terms and reduced ones:
    for (s_it = lit2term[ cSeg ].begin(); s_it != s_end; s_it++)
    {
      std::set<int> &temp = nsubMods[ *s_it ].literals;
      Scalar score        = nsubMods[ *s_it ].score;
      assert( *(temp.find(cSeg)) == cSeg );

      temp.erase( cSeg );

      // non-submod edges:
      std::set<int>::const_iterator t_it( temp.begin() ), t_end ( temp.end() ); 
      literalTerm<Scalar> in;
      in.score = -score;
      in.literals.insert( newVar );
      // add to in all literals

      for (; t_it != t_end; t_it++)
        in.literals.insert( *t_it );
      subMods.push_back(in);

      ////// update other variables involved in term, if term becomes a single variable: (0
      if (temp.size() == 1) // term vanishes
      {
        int var = *(temp.begin());
        b1.push_back(std::pair<int, Scalar>( var, score) );

        // make elem var occur in one less term (which has disappered now)
        HeapStoreEntry& hse = _heapEntryStore[ var ];
        hse.setCost( hse.cost()-1 );
        lit2term[var].erase( *s_it );
        _heap.update( HeapEntryB( var ) );
      }
    }
  }
}


/// nsm already contains a vector of seg,label pairs and a score now convert those to positive literal terms
template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  convertToTerms ( higherOrders<Scalar> nsm, std::vector<literalTerm<Scalar> > &subMods, std::vector<literalTerm<Scalar> > &nsubMods )
{
  typename std::vector<literalTerm<Scalar> > in;
  literalTerm<Scalar> lTerm;
  lTerm.score = nsm.score;

  if ( fabs(lTerm.score) < 0.0001 ) return;

  recursiveConversion ( lTerm, nsm, 0, in );
  for(int i = 0; i<in.size(); i++ )
    if (in[i].score > 0 && in[i].literals.size() >1 )
      nsubMods.push_back(in[i]);
    else
      subMods.push_back(in[i]);
}


template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  recursiveConversion ( literalTerm<Scalar> lTerm, higherOrders<Scalar>& nsm, int position, std::vector<literalTerm<Scalar> > &vec_lit )
{
  //    std::vector< literalTerm<ScalarT> > vec_lit ;
  if (position == nsm.vars.size() || position > maxRecursionDepth) // else explosion in variables and in memory!
  {
    if ( lTerm.literals.size() > 0)
      vec_lit.push_back( lTerm );
    return;
  }
  if ( nsm.vars[position].second == 1 ) //0 )
  {
    lTerm.literals.insert( nsm.vars[position].first );
    recursiveConversion ( lTerm, nsm, position+1, vec_lit );
  }
  else
  {
    // (1-xi) * rest = rest - xi*rest

    //rest, shortened by one literal
    recursiveConversion ( lTerm, nsm, position+1, vec_lit );

    // - xi*rest
    lTerm.score = - lTerm.score;
    lTerm.literals.insert( nsm.vars[position].first );
    recursiveConversion ( lTerm, nsm, position+1, vec_lit );
  }
}

template<class Scalar>
void
  HigherOrderConverter<Scalar>::
  getSmallerIterator( 
  typename std::vector< higherOrders<Scalar> >::iterator& nsm_it, 
  typename std::vector< higherOrders<Scalar> >::iterator& sm_it, 
  typename std::vector< higherOrders<Scalar> >::iterator& nsm_end, 
  typename std::vector< higherOrders<Scalar> >::iterator& sm_end, 
  typename std::vector< higherOrders<Scalar> >::iterator*& it, 
  typename std::vector< higherOrders<Scalar> >::iterator*& end)
{
  if (nsm_it == nsm_end)
  {
    it = &sm_it;end = &sm_end;
    return;
  }
  if (sm_it == sm_end )
  {
    it = &nsm_it;end = &nsm_end;
    return;
  }
  if ( sm_it->vars[0].first > nsm_it->vars[0].first || (sm_it->vars[0].first == nsm_it->vars[0].first && sm_it->vars[0].second > nsm_it->vars[0].second))
  {
    it = &nsm_it;end = &nsm_end;
    return;
  }
  else
  {
    it = &sm_it;end = &sm_end;
    return;
  }
}

#endif
