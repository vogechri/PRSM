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
#ifndef __BinaryConversion_CPP_
#define __BinaryConversion_CPP_

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
BinaryConverter<Scalar>::
addOcclusionListOuterPixels( std::list< occPixelList< Scalar > >& _occListOuter )
  {
    occlusionListOuter.clear();
    occlusionListOuter.reserve( _occListOuter.size() );
    fillFromPixelList ( occlusionListOuter, _occListOuter ,0, 0 );
    std::sort( occlusionListOuter.begin(), occlusionListOuter.end() );    
  }

template<class Scalar>
void
BinaryConverter<Scalar>::
addOcclusionListOuterPixels( std::list< occPixelList< Scalar > >& _occListOuter, 
     std::list< occPixelList< Scalar > >& _occListOuter2)
  {
    // merge first, treat as normal lists:
    occlusionListOuter.reserve( _occListOuter.size() + _occListOuter2.size() );
    std::vector< occSegList<Scalar> > occListOuter;
    std::vector< occSegList<Scalar> > occListOuter2;

    fillFromPixelList ( occListOuter, _occListOuter , 0, 0 );
    fillFromPixelList ( occListOuter2, _occListOuter2, 0, 0 );
    std::sort( occListOuter.begin(), occListOuter.end() );
    std::sort( occListOuter2.begin(), occListOuter2.end() );

    // merge always, since 01 is already in the list 
    mergeListsOuter( occListOuter, occListOuter2, occlusionListOuter );
   };



template<class Scalar>
void
BinaryConverter<Scalar>::
fillFromPixelList ( std::vector< occSegList<Scalar> >&occlusionList, std::list< occPixelList< Scalar > >& _occList0, int ownlbl, int startId )
  {
    typename std::list< occPixelList< Scalar > >::iterator o_it (_occList0.begin()), o_end(_occList0.end());
    for ( ; o_it != o_end ; o_it++ )
    {
      std::vector < typename occSegList<Scalar>::SegLabelPair > _occluders;
      std::vector < PixSegPair >& occ_list = o_it->occluders;
      int ownid  = o_it->ownId;

      int _fullOcclusion = (o_it->fullOccluder>=0) ? 1:0;

      if ( !_fullOcclusion )
      {
        assert( !startId || ownid == occ_list[0].first);
        for ( int i =startId;i < occ_list.size(); i++ )
          _occluders.push_back( typename occSegList<Scalar>::SegLabelPair( occ_list[i].first, occ_list[i].second ) );
        std::sort( _occluders.begin(), _occluders.end(), sortOccluders<Scalar> );
      }
      occSegList<Scalar> newEntry ( -1, ownid, ownlbl, _occluders, _fullOcclusion );

      occlusionList.push_back( newEntry );
    }
  }


// the same for the outer occlusion list !
template<class Scalar>
void
BinaryConverter<Scalar>::
generateHigherOrderTerms_outer( std::vector<Scalar>& globalDataScores )
  {
    int nPixel = globalDataScores.size()/2;
    for (int i = 0; i < occlusionListOuter.size(); i++)
    {
      int sid   = occlusionListOuter[i].segId;
      int label = occlusionListOuter[i].label;

      Scalar penalty = occPenalty - globalDataScores[sid];

      if ( !occlusionListOuter[i].fullOcclusion && occlusionListOuter[i].occluders.size() >1 )
      {
        higherOrders<Scalar> ho(-penalty);
        ho.vars.reserve  ( occlusionListOuter[i].occluders.size() );
        for (int k=0; k < occlusionListOuter[i].occluders.size(); k++)
          ho.vars.push_back( typename higherOrders<Scalar>::SegLabelPair ( occlusionListOuter[i].occluders[k].first, 1-occlusionListOuter[i].occluders[k].second));// twist here to biuld the terms

        if (ho.isSubMod())
          submods.push_back(ho);
        else
          nonsubmods.push_back(ho);
      }
      else
      {
         if (occlusionListOuter[i].occluders.size() == 1)
          unaries.push_back (Unaries<Scalar> (occlusionListOuter[i].occluders[0].second, occlusionListOuter[i].occluders[0].first, penalty) );
      }
    }
  }


// generate the unaries to be added, split up into submodular elements and non-submods
template<class Scalar>
void
BinaryConverter<Scalar>::
generateHigherOrderTerms()
  {

    unaries.clear();
    submods.clear();
    nonsubmods.clear();

    for (int i = 0; i < occlusionList.size(); i++)
    {
      int sid = occlusionList[i].segId;

      int label = occlusionList[i].label;
      Scalar penalty;

      Scalar xtraPen =0;
      if (!label)
        penalty = occPenalty * freePix[sid + label*nSegments] - dataPen0[sid] + xtraPen;
      else
        penalty = occPenalty * freePix[sid + label*nSegments] - dataPen1[sid] + xtraPen;

      // now unary:
      unaries.push_back (Unaries<Scalar> (label, sid, penalty) );

      if ( !occlusionList[i].fullOcclusion )
      {
        higherOrders<Scalar> ho(-penalty);
        ho.vars.reserve  ( occlusionList[i].occluders.size()+1 );
        ho.vars.push_back( typename higherOrders<Scalar>::SegLabelPair ( occlusionList[i].segId, occlusionList[i].label ) );

        for (int k=0; k < occlusionList[i].occluders.size(); k++)
          ho.vars.push_back( typename higherOrders<Scalar>::SegLabelPair ( occlusionList[i].occluders[k].first, 1-occlusionList[i].occluders[k].second));// twist here to biuld the terms

        if (ho.isSubMod())
          submods.push_back(ho);
        else
          nonsubmods.push_back(ho);
      }
    }
  }


// merge 2 lists as two views combine, MUST check if no free Vars: ignore !
template<class Scalar>
void
BinaryConverter<Scalar>::
mergeLists ( std::vector< occSegList<Scalar> >& oListsa, std::vector< occSegList<Scalar> >& oListsb, 
             std::vector< occSegList<Scalar> >& mergedList ) //, std::vector<int>& freePixel )
  {

    mergedList.clear();

    typename std::vector<occSegList<Scalar> >::const_iterator o_it1 (oListsa.begin()), o_end1 (oListsa.end() );
    typename std::vector<occSegList<Scalar> >::const_iterator o_it2 (oListsb.begin()), o_end2 (oListsb.end() );

    // merge like a zipper into oLists0
    for (; o_it1 != o_end1 && o_it2 != o_end2;)
    {
      if (*o_it1 < *o_it2)// unique entry in first list
      {
        if ( freePix[o_it1->segId + o_it1->label*nSegments] >0 )
          mergedList.push_back( *o_it1 );    
        o_it1++;
        continue;
      }
      if (*o_it1 > *o_it2)// unique entry in first list
      {
        if ( freePix[o_it2->segId + o_it2->label*nSegments] >0 )
          mergedList.push_back( *o_it2 );
        o_it2++;
        continue;
      }
      if (*o_it1 == *o_it2) // merging both lists
      {
        if ( freePix[o_it1->segId + o_it1->label*nSegments] <=0 )
        {
          o_it1++;o_it2++;
          continue;
        }

        if (o_it1->fullOcclusion)
        {
          mergedList.push_back( *o_it1 );
          o_it1++;o_it2++;
          continue;
        }
        else if (o_it2->fullOcclusion)
        {
          mergedList.push_back( *o_it2 );
          o_it1++;o_it2++;
          continue;
        }
        /////////////////////////

        // assume these are sorted to begin with by segid
        // merge the occluders in the list: HOW? s set ? also checkk for fullOcclusion here
        std::vector < typename occSegList< Scalar>::SegLabelPair > mergedVec;
        mergedVec.reserve( (o_it1->occluders).size()+(o_it2->occluders).size() );
        typename std::vector < typename occSegList<Scalar>::SegLabelPair >::const_iterator l_it1 ((o_it1->occluders).begin()), l_end1 ((o_it1->occluders).end());
        typename std::vector < typename occSegList<Scalar>::SegLabelPair >::const_iterator l_it2 ((o_it2->occluders).begin()), l_end2 ((o_it2->occluders).end());
        // skip first entry (identical in both lists)
        // merge like a zipper into oLists0
        int fullOcclusion = 0;
        for (; l_it1 != l_end1 && l_it2 != l_end2;)
        {
          if (l_it1->first < l_it2->first)
          {
            mergedVec.push_back( *l_it1 );
            l_it1++;
            continue;
          }
          if (l_it1->first > l_it2->first)
          {
            mergedVec.push_back( *l_it2 );
            l_it2++;
            continue;
          }
          if (l_it1->first == l_it2->first)
          {
            if ( l_it1->second == l_it2->second) // twice keep only one
            {
              mergedVec.push_back( *l_it2 );
              l_it1++;l_it2++;
            }
            else // special case: fullOcclusion - in the list1 case add both (since both are one never ending here!)
            {
              fullOcclusion = 1;
              mergedVec.clear(); // not important any more
              break;
            }
          }
        }
        if (!fullOcclusion)
        {
          // insert the rest
          for (; l_it1 != l_end1;l_it1++)
            mergedVec.push_back( *l_it1 );
          // insert the rest
          for (; l_it2 != l_end2;l_it2++)
            mergedVec.push_back( *l_it2 );
        }

        mergedList.push_back( occSegList<Scalar>( o_it1->mpId, o_it1->segId, o_it1->label, mergedVec, fullOcclusion ) );
        o_it1++;o_it2++;
      }
    }

    // insert the rest:
    for (; o_it1 != o_end1; o_it1++)
      mergedList.push_back( *o_it1 );

    for (; o_it2 != o_end2; o_it2++)
      mergedList.push_back( *o_it2 );
  }


/// merge 2 lists as two views combine, MUST check if no free Vars: ignore !
template<class Scalar>
void
BinaryConverter<Scalar>::
mergeListsOuter ( std::vector< occSegList<Scalar> >& oListsa, std::vector< occSegList<Scalar> >& oListsb, 
                  std::vector< occSegList<Scalar> >& mergedList )
  {

    mergedList.clear();

    typename std::vector< occSegList<Scalar> >::const_iterator o_it1 (oListsa.begin()), o_end1 (oListsa.end() );
    typename std::vector< occSegList<Scalar> >::const_iterator o_it2 (oListsb.begin()), o_end2 (oListsb.end() );

    // merge like a zipper into oLists0
    for (; o_it1 != o_end1 && o_it2 != o_end2;)
    {
      if (*o_it1 < *o_it2)// unique entry in first list
      {
        mergedList.push_back( *o_it1 );    
        o_it1++;
        continue;
      }
      if (*o_it1 > *o_it2)// unique entry in first list
      {
        mergedList.push_back( *o_it2 );
        o_it2++;
        continue;
      }
      if (*o_it1 == *o_it2) // merging both lists
      {
        if (o_it1->fullOcclusion)
        {
          mergedList.push_back( *o_it1 );
          o_it1++;o_it2++;
          continue;
        }
        else if (o_it2->fullOcclusion)
        {
          mergedList.push_back( *o_it2 );
          o_it1++;o_it2++;
          continue;
        }
        /////////////////////////

        // assume these are sorted to begin with by segid
        // merge the occluders in the list: HOW? s set ? also checkk for fullOcclusion here
        typename std::vector < typename occSegList<Scalar>::SegLabelPair > mergedVec;
        mergedVec.reserve( (o_it1->occluders).size()+(o_it2->occluders).size() );
        typename std::vector < typename occSegList<Scalar>::SegLabelPair >::const_iterator l_it1 ((o_it1->occluders).begin()), l_end1 ((o_it1->occluders).end());
        typename std::vector < typename occSegList<Scalar>::SegLabelPair >::const_iterator l_it2 ((o_it2->occluders).begin()), l_end2 ((o_it2->occluders).end());
        // skip first entry (identical in both lists)
        // merge like a zipper into oLists0
        int fullOcclusion = 0;
        for (; l_it1 != l_end1 && l_it2 != l_end2;)
        {
          if (l_it1->first < l_it2->first)
          {
            mergedVec.push_back( *l_it1 );
            l_it1++;
            continue;
          }
          if (l_it1->first > l_it2->first)
          {
            mergedVec.push_back( *l_it2 );
            l_it2++;
            continue;
          }
          if (l_it1->first == l_it2->first)
          {
            if ( l_it1->second == l_it2->second) // twice keep only one
            {
              mergedVec.push_back( *l_it2 );
              l_it1++;l_it2++;
            }
            else // special case: fullOcclusion - in the list1 case add both (since both are one never ending here!)
            {
              fullOcclusion = 1;
              mergedVec.clear(); // not important any more
              break;
            }
          }
        }
        if (!fullOcclusion)
        {
          // insert the rest
          for (; l_it1 != l_end1;l_it1++)
            mergedVec.push_back( *l_it1 );
          // insert the rest
          for (; l_it2 != l_end2;l_it2++)
            mergedVec.push_back( *l_it2 );
        }

        mergedList.push_back( occSegList<Scalar>( o_it1->mpId, o_it1->segId, o_it1->label, mergedVec, fullOcclusion ) );
        o_it1++;o_it2++;
      }
    }

    // insert the rest:
    for (; o_it1 != o_end1; o_it1++)
      mergedList.push_back( *o_it1 );

    for (; o_it2 != o_end2; o_it2++)
      mergedList.push_back( *o_it2 );
  }

#endif