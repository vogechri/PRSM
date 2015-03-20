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
//////////////////////////////////////////////////////
////////// Global store normals, rotation  ///////////
//////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANOR_STORE__
#define __ENERGY_ROTTRANOR_STORE__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#ifndef _NO_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

/// provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
template<typename Scalar> class EvalEnergyStore
{
public:

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::Mat3x3T<Scalar>    M3;

  EvalEnergyStore(): normals(), tra(), rot()
  {};

  ~EvalEnergyStore(){};
  
  const std::vector<P3>* getNormals( ) {return &normals;};
  const std::vector<P3>* getTra( )     {return &tra;};
  const std::vector<M3>* getRot( )     {return &rot;};

  std::vector<P3>* manipulateNormals( ) {return &normals;};
  std::vector<P3>* manipulateTra( )     {return &tra;};
  std::vector<M3>* manipulateRot( )     {return &rot;};

  void setRotTra( std::vector<P3>* tra_, std::vector<P3>* rot_ ) {tra = tra_;rot = rot_;}

  void prepareNormals( Scalar* normals_, int nElem, const int nStep=3 )
  {
    normals.clear();    
    normals.reserve(nElem);

    for (int i = 0; i < nElem ;i++)
      normals.push_back(P3(&(normals_[ nStep*i]) ));

  }

  void flipNormals( )
  {
    for (int i = 0; i < normals.size() ;i++)
      normals[i] = normals[i]*Scalar(-1);
  }

  /// from a array of translations  to a P3 vector
  void prepareTra( Scalar* dTra, int nElem, const int nStep=3 )
  {
    tra.clear();
    tra.reserve(nElem);

    for (int i=0; i < nElem; i++)
      tra.push_back(P3( &(dTra[nStep*i]) ));
  }

  void prepareRotTra( Scalar* dRt, int nElem )
  {
    rot.clear();
    tra.clear();
    rot.reserve(nElem);
    tra.reserve(nElem);

    for (int i=0; i < nElem; i++)
    {
      P3 ri_x = P3 (&( dRt[ 6*i   ] ));
      P3 ti_x = P3 (&( dRt[ 6*i +3] ));

      Scalar sinA = ri_x.norm();
      sinA = min(sinA, 1.);sinA = max(sinA, 0.000000001);
      
      P3 rVec       = ri_x * (1./sinA);
      Scalar alpha  = asin(sinA);
      Scalar cosA   = cos(alpha);

      Scalar cosA_1 = Scalar(1) - cosA;

      M3 A1( cosA,   -ri_x[2],  ri_x[1],
             ri_x[2], cosA,    -ri_x[0],
             -ri_x[1], ri_x[0],  cosA );

      M3 A2( cosA_1 * rVec[0]*rVec[0], cosA_1 * rVec[0]*rVec[1] ,cosA_1 * rVec[0]*rVec[2],
             cosA_1 * rVec[0]*rVec[1], cosA_1 * rVec[1]*rVec[1], cosA_1 * rVec[2]*rVec[1],
             cosA_1 * rVec[0]*rVec[2], cosA_1 * rVec[2]*rVec[1], cosA_1 * rVec[2]*rVec[1] );

      rot.push_back( A1+A2 );
      tra.push_back( ti_x  );
    }
  }

  /// from 3 parametric rotations to full matrix 
  void prepareRot( Scalar* dRot, int nRot )
  {
    rot.clear();
    rot.reserve(nRot);

    for (int i=0; i < nRot; i++)
    {
      P3 ri_x = P3 ( &( dRot[ 3*i ] ) );

      Scalar sinA = ri_x.norm();
      sinA = min(sinA, 1.);sinA = max(sinA, 0.000000001);
      
      P3 rVec = ri_x * (1./sinA);
      Scalar alpha = asin(sinA);
      Scalar cosA = cos(alpha);

      Scalar cosA_1 = Scalar(1) - cosA;

      M3 A1( cosA,   -ri_x[2],  ri_x[1],
             ri_x[2], cosA,    -ri_x[0],
             -ri_x[1], ri_x[0],  cosA );

      M3 A2( cosA_1 * rVec[0]*rVec[0], cosA_1 * rVec[0]*rVec[1] ,cosA_1 * rVec[0]*rVec[2],
             cosA_1 * rVec[0]*rVec[1], cosA_1 * rVec[1]*rVec[1], cosA_1 * rVec[2]*rVec[1],
             cosA_1 * rVec[0]*rVec[2], cosA_1 * rVec[2]*rVec[1], cosA_1 * rVec[2]*rVec[1] );

      rot.push_back( A1+A2 );
    }
  }

  /// from 4x4 array to rotation matrix, translation vectors
  void setRotTra( Scalar* dRot, int nElem )
  {
    rot.clear();
    tra.clear();
    rot.reserve(nElem);
    tra.reserve(nElem);

    for (int i=0; i < nElem; i++)
    {
      const Scalar*  Rt_i  = &( dRot[16*i   ] );
      tra.push_back( P3( &(Rt_i[12]) ) );
      rot.push_back(  M3(Rt_i[0], Rt_i[4], Rt_i[8], Rt_i[1], Rt_i[5], Rt_i[9], Rt_i[2], Rt_i[6], Rt_i[10]) );
    }
  }


private:

  /// normals
  std::vector<P3> normals;
  /// translations
  std::vector<P3> tra;
  /// rotations
  std::vector<M3> rot;
};
#endif
