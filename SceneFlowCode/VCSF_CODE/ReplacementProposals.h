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
////////////////////////////////////////////////////////////////////////
// Generate a proposal set by mixing neighboring normals/rigid motions//
////////////////////////////////////////////////////////////////////////

#ifndef __Replacement_Proposals__h
#define __Replacement_Proposals__h

#include<set>
#include<map>
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "EvalEnergyStore.h"

using namespace std;
using namespace Math;


// output: in this form: 
//  int *prop2PlaneMap(NULL);
//  Scalar* expCenters (NULL);
// extend
// normals,    rotations,     translations
// in EvalEnergyStore

// locally to handle - later outside with all/new ones
//        EvalEnergyStore<Scalar> enStore;
//        enStore.prepareNormals( normals, nNormals, nStep );
//        enStore.setRotTra( rotations, nNormals );
// so EvalEnergyStore<Scalar> enStore; is an input 

// step1: identify spots for replacement == at solution boundaries 
// identify pairs of replacement proposals and their set of segments
// N replacement proposals and this proposal at all segments involved as proposal with id
// these can later be merged ?!
// need a map from pair of ints to set of segments 

// call: LocalReplacement lrp(esStore, segImg[0][0], nSegments[0][0], currestSolution[0][0], centers[0][0], edges[0][0]);
// prop2PlaneMap =  lrp.getPropMap();
// expCenters = lrp.getExpCenters();

template<typename Scalar> 
class LocalReplacement
{
  typedef std::pair<int,int> propPair;
  typedef std::set<int> segmentSet;
  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 15> P15;

public:

  LocalReplacement( EvalEnergyStore<Scalar>& _enstore, int* _segImg, int _nSegments, int* _currentSolution, Scalar* _segCenters, const mxArray* _edges, Scalar* K_ = NULL ) 
    : enStore(_enstore), nSegments(_nSegments), currentSolution(_currentSolution), segCenters( _segCenters ), edges(_edges)
  {
    clear();
    findReplacementIds(edges);
    genReplacementProps(K_);
  };

  ~LocalReplacement(){};

  void clear()
  { 
    propPairSet.clear();//mvpCenters.clear();
    prop2PlaneMap.clear();expCenters.clear();
  };

  /// idea is to gen new ones but also old props if located at new spots, iterate twice, redo the whole thing BUT only new mvps/locations
  void generateMore( Scalar* K_ = NULL )
  { 
    // push_back info so far..
    /*
    std::map<propPair, segmentSet >::iterator it(propPairSet.begin()), it_end(propPairSet.end())  ;
    
    for (;it != it_end;it++)
    {
      std::pair< std::map<propPair, segmentSet >::iterator,bool > ret;
      ret = oldPropPairSet.insert ( *it );
      if (!ret.second)
      {
        ((ret.first)->second).insert( it->second );
      }
    }
*/
    findReplacementIds(edges);
    genReplacementProps(K_);  
  };

  void findReplacementIds( const mxArray* edges )
  {
    // run along segments and neighs
    std::vector<P3>* norStored = enStore.manipulateNormals();
    std::vector<P3>* traStored = enStore.manipulateTra( );
    std::vector<M3>* rotStored = enStore.manipulateRot( );

    for (int sid = 0;sid < nSegments;sid++)
    {
      int pid = currentSolution[ sid ] ;// proposal in segment segi

      oldPropPairSet[pid].insert( sid ); // prposal pid on segment sid

      M3 r1 = (*rotStored)[ pid ];
      P3 t1 = (*traStored)[ pid ];
      P3 n1 = (*norStored)[ pid ];
      P15 mvpP1 = phys_mvp(n1,t1,r1);
      typename std::map<P15, int>::iterator it( mvp2Pid.find( mvpP1 ));
      if (it == mvp2Pid.end() )
        mvp2Pid[ mvpP1 ] = pid;
      else
        if (it->second != pid)
          printf("Hmm?\n");

      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges, sid) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges, sid) );

      for (int j =0; j < nIds ;j++)
      {
        int sid2 = (int) (edge [5*j ])-1;
        int pid2 = currentSolution[ sid2 ] ;// proposal in segment segi
        //          typename EvalEnergyStore<Scalar>::P3 centerID ( &centers[3 * id] );

        if ( pid != pid2 && sid < sid2 ) //&& ((centerSegi - centerID).sqrnorm() < maxDist*maxDist) )
        {
          //            valids[id] = 0;
          //            closeIds.push_back( id );
          //            conquering.push_back( id );

          propPair temp ( min(pid,pid2), max(pid,pid2) );//some order needed
          std::map<propPair, segmentSet >::iterator it = propPairSet.find( temp );

          propPairSet[ temp ].insert( sid  ); 
          propPairSet[ temp ].insert( sid2 ); 

          /*
          if (propPairSet.end() != it)
          {
          propPairSet[ temp ].insert( sid  ); 
          propPairSet[ temp ].insert( sid2 ); 
          }
          else // not new
          {
          propPairSet[ temp ].insert( sid  ); 
          propPairSet[ temp ].insert( sid2 );             
          }
          */
        }
      }
    }
  }


  void genReplacementProps( Scalar* K_ = NULL)
  {
    // can actually manipute this right here:
    //  std::vector<int> propIds(nSegments, 0);
    //  std::vector<int> valids (nSegments, 1);
    //  std::vector<int> result (nSegments, 0);

    M3 K;
    if (K_ != NULL)
      K=M3(K_);

    std::vector<P3>* norStored = enStore.manipulateNormals();
    std::vector<P3>* traStored = enStore.manipulateTra( );
    std::vector<M3>* rotStored = enStore.manipulateRot( );


    //norStored->resize(nSegments);
    //traStored->resize(nSegments);
    //rotStored->resize(nSegments);

    std::map<propPair, segmentSet >::iterator it(propPairSet.begin()), it_end(propPairSet.end());
    for(;it!=it_end;it++)
    {
      //
      propPair    pps = it->first;
      segmentSet  sgs = it->second;

      M3 r1 = (*rotStored)[ pps.first ];
      P3 t1 = (*traStored)[ pps.first ];
      P3 n1 = (*norStored)[ pps.first ];
      M3 r2 = (*rotStored)[ pps.second ];
      P3 t2 = (*traStored)[ pps.second ];
      P3 n2 = (*norStored)[ pps.second ];

      // 1. if already exists do not build new one -- could even do only in new locations by merging old and new sets below
      int pid1(-1), pid2(-1); 
      P15 mvpP1 = phys_mvp(n1,t2,r2);
      bool new1 (false), new2(false), new3(false), new4(false);
      typename std::map<P15, int>::iterator it( mvp2Pid.find( mvpP1 ));
      if (it == mvp2Pid.end() )
      {
        // simplistic yet can extend later / refit or so
        norStored->push_back( n1 );
        traStored->push_back( t2 );
        rotStored->push_back( r2 );
        pid1 = norStored->size()-1;
        mvp2Pid[mvpP1] = pid1;
        new1 = true;
      }
      else
        pid1 = it->second;

      P15 mvpP2 = phys_mvp(n2,t1,r1);
      it = mvp2Pid.find( mvpP2 );
      if (it == mvp2Pid.end() )
      {
        norStored->push_back( n2 );
        traStored->push_back( t1 );
        rotStored->push_back( r1 );
        pid2 = norStored->size()-1;
        mvp2Pid[mvpP2] = pid2;
        new2 = true;
      }
      else
        pid2 = it->second;
///
    int pid3(-1), pid4(-1);
    if (K_ != NULL)
    {
      /*
      P3 n3 = refit(r1, t1, n1, K, P3( &(segCenters[3*(*(sgs.begin())) ]) ), n2, t2, r2 );
      P3 n4 = refit(r2, t2, n2, K, P3( &(segCenters[3*(*(sgs.begin())) ]) ), n1, t1, r1 );

      norStored->push_back( n3 );
      traStored->push_back( t1 );
      rotStored->push_back( r1 );
      pid3 = norStored->size()-1;

      norStored->push_back( n4 );
      traStored->push_back( t2 );
      rotStored->push_back( r2 );
      pid4 = norStored->size()-1;
*/
      P3 t1out, t2out;
      refitC(r1, t1, n1, r2, t2, n2, K, P3( &(segCenters[3*(*(sgs.begin())) ]) ), P3( &(segCenters[3*(*(sgs.begin())) ]) ), t1out, t2out );

      P15 mvpP3 = phys_mvp(n1,t2out,r2);
      it = mvp2Pid.find( mvpP3 );
      if (it == mvp2Pid.end() )
      {
        norStored->push_back( n1 );
        traStored->push_back( t2out );
        rotStored->push_back( r2 );
        pid3 = norStored->size()-1;
        new3 = true;
      }
      else
        pid3 = it->second;

      P15 mvpP4 = phys_mvp(n2,t1out,r1);
      it = mvp2Pid.find( mvpP4 );
      if (it == mvp2Pid.end() )
      {
        norStored->push_back( n2 );
        traStored->push_back( t1out );
        rotStored->push_back( r1 );
        pid4 = norStored->size()-1;
        new4 = true;
      }
      else
        pid4 = it->second;
    }
///
      P3 mvpCenter(0,0,0);
      segmentSet::iterator sit(sgs.begin()), send(sgs.end());
      for(; sit!=send; sit++)
      {
        mvpCenter += P3 ( &(segCenters[3*(*sit) ]) );

        if ( oldPropPairSet[pid1].find( *sit ) == oldPropPairSet[pid1].end() )
        {
          // both proposals are located at the same segment
          prop2PlaneMap.push_back(pid1);
          expCenters.push_back( segCenters[3*(*sit) ] );
          expCenters.push_back( segCenters[3*(*sit)+1 ] );
          expCenters.push_back( segCenters[3*(*sit)+2 ] );

          oldPropPairSet[pid1].insert( (*sit) ); // prposal pid on segment sid
        }

        if ( oldPropPairSet[pid2].find( *sit ) == oldPropPairSet[pid2].end() )
        {
          prop2PlaneMap.push_back(pid2);
          expCenters.push_back( segCenters[3*(*sit) ] );
          expCenters.push_back( segCenters[3*(*sit)+1 ] );
          expCenters.push_back( segCenters[3*(*sit)+2 ] );

          oldPropPairSet[pid2].insert( (*sit) ); // prposal pid on segment sid
        }

        if (K_ != NULL)
        {
          if ( oldPropPairSet[pid3].find( *sit ) == oldPropPairSet[pid3].end() )
          {
            // both proposals are located at the same segment
            prop2PlaneMap.push_back(pid3);
            expCenters.push_back( segCenters[3*(*sit) ] );
            expCenters.push_back( segCenters[3*(*sit)+1 ] );
            expCenters.push_back( segCenters[3*(*sit)+2 ] );
            oldPropPairSet[pid3].insert( (*sit) ); // prposal pid on segment sid
          }
          if ( oldPropPairSet[pid4].find( *sit ) == oldPropPairSet[pid4].end() )
          {
            // both proposals are located at the same segment
            prop2PlaneMap.push_back(pid4);
            expCenters.push_back( segCenters[3*(*sit) ] );
            expCenters.push_back( segCenters[3*(*sit)+1 ] );
            expCenters.push_back( segCenters[3*(*sit)+2 ] );
            oldPropPairSet[pid4].insert( (*sit) ); // prposal pid on segment sid
          }
/*
          prop2PlaneMap.push_back(pid3);
          prop2PlaneMap.push_back(pid4);
          expCenters.push_back( segCenters[3*(*sit) ] );
          expCenters.push_back( segCenters[3*(*sit)+1 ] );
          expCenters.push_back( segCenters[3*(*sit)+2 ] );
          expCenters.push_back( segCenters[3*(*sit) ] );
          expCenters.push_back( segCenters[3*(*sit)+1 ] );
          expCenters.push_back( segCenters[3*(*sit)+2 ] );
          */
        }
      }

      mvpCenter /= sgs.size();

      if (new1)
      {
        mvpCenters.push_back( mvpCenter );
      }
      if (new2)
      {
        mvpCenters.push_back( mvpCenter );
      }
      if (K_ != NULL)
      {
        if (new3)
          mvpCenters.push_back( mvpCenter );
        if (new4)
          mvpCenters.push_back( mvpCenter );
      }
    }
  }

  int*    getPropMap    () {return (int*) (&( prop2PlaneMap[0] ));};
  Scalar* getExpCenters () {return (Scalar*) (&( expCenters[0] ));};
  int     getNProposals () {return prop2PlaneMap.size();};
  std::vector<P3>& getMvpCenters() {return mvpCenters;}

private:

  /// compute new normal -- s.t. old 2d motion remains, recomputes n_start as if R,t is used
  P3 refit(M3 R, P3 t, P3 N, const M3& K, P3 center_start, P3 n_start, P3 t_start, M3 R_start )
  {
    //     typename EvalEnergyStore<Scalar>::P3 center ( &centers[3 * id] );
    //     typename EvalEnergyStore<Scalar>::P3 n_start  = (*norStored)[currentSolution[id]];

    P3 kt       = K * t; // from new
    P3 c_3d     = center_start / (center_start|n_start);
    Scalar dStart = c_3d[2];
    if ( dStart != dStart )
      printf ("Fail generateNeighProposals_refit %f \n", dStart);

    P3 c_2d     = K * R * c_3d;  c_2d /= dStart;
    P3 c_2dG    = K * ( t_start + R_start *c_3d );c_2dG /= c_2dG[2];
    if ( c_2dG[2] != c_2dG[2] )
      printf ("Fail generateNeighProposals_refit %f \n", c_2dG[2]);

    P3 u = c_2d-c_2dG*c_2d[2];
    P3 v = kt - c_2dG * kt[2];// can be 0
//    Scalar s = -(u[0]*v[0] + u[1]*v[1]) / (v[0]*v[0]+v[1]*v[1]);
    Scalar s = -(u[0] + u[1]) / (v[0]+v[1]);

//    if ( fabs(v[0]*v[0]+v[1]*v[1]) < 0.0000001 )
    if ( fabs(v[0]+v[1]) < 0.0000001 )
      s=1;
    if ( s != s )
      printf ("Fail generateNeighProposals_refit %f \n", s);

    P3 N_scale = n_start * (s*dStart);

    // test: 
    P3 c_3dt     = center_start / (center_start|N_scale);
    P3 c_2dGt    = K * ( t + R *c_3d );c_2dGt /= c_2dGt[2];

    if ( fabs(s*dStart) < 0.0000001 )
      N_scale = n_start;

    return N_scale;
  }

  P15 phys_mvp(P3& n1,P3& t1,M3& r1)
  {
    return P15 ( n1[0], n1[1], n1[2], t1[0], t1[1], t1[2], r1(0,0), r1(0,1), r1(0,2), r1(1,0), r1(1,1), r1(1,2), r1(2,0), r1(2,1), r1(2,2) );
  }


  /// M3& r1out, M3& r2out, centers equal if normals differ still xchange
  void refitC(M3 R1, P3 t1, P3 n1, M3 R2, P3 t2, P3 n2, const M3& K, P3 center1, P3 center2, P3& t1out, P3& t2out )
  {
    P3 c_3d1    = center1 / (center1|n1);
    P3 c_3d2    = center2 / (center2|n2);     

    // Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) - Rt_lin(1:3,1:3,i) * centers(:,i) + centers(:,i);
   // back :     Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) + Rt_lin(1:3,1:3,i) * centers(:,i) - centers(:,i);

    // first convert to centered version:
    t1out = t1 - R1 * c_3d1 + c_3d1;
    t2out = t2 - R2 * c_3d2 + c_3d2;

    // xchange and convert back to origin centered version: R1 R2 swapped ? 
    P3 temp = t2out + R2 * c_3d1 - c_3d1;
    t1out   = t1out + R1 * c_3d2 - c_3d2;
    t2out = temp;

//    P3 temp = t2out + R2 * c_3d1 - c_3d1;
//    t2out   = t1out + R1 * c_3d2 - c_3d2;
//    t1out = temp;

//    t1out = t2 - R2 * (c_3d2-c_3d1) + (c_3d2 - c_3d1);
//    t2out = t1 - R1 * (c_3d1-c_3d2) + (c_3d1 - c_3d2);
  }


  std::map<propPair, segmentSet > propPairSet;

  /// contains the physical moving planes
  EvalEnergyStore<Scalar>& enStore;

  /// per esgment a proposal
  int* currentSolution;

  Scalar* segCenters;

  const mxArray* edges;

  int nSegments;

  // output: in this form: 
  std::vector<int>    prop2PlaneMap;
  std::vector<Scalar> expCenters;

  /// the mvp receive a true center -- for later
  std::vector<P3>     mvpCenters;

  /// the mvp receive a true center -- for later
//  std::vector<P3>     allMvpCenters;

  /// some stupid map (not even hash) from n,r,t to id, so plane2prop
  std::map<P15, int>  mvp2Pid;

  /// from propId to all covered segments so far, thus can check if new should be covered and where
  std::map<int, segmentSet > oldPropPairSet;
};

#endif //__Replacement_Proposals__h