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
////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#ifndef __DATA_Containers__h
#define __DATA_Containers__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;

template<typename Scalar> struct Binary
{
  typedef Math::VectorT<Scalar, 4>  P4;
  Binary(int _segI, int _segJ, P4 _edge): segJ(_segJ), segI(_segI), edge(_edge) {};
  Binary(int _segI, int _segJ): segJ(_segJ), segI(_segI) {};
  Binary(): segJ(-1), segI(-1) {};
  ////////
  P4 edge;
  int segI;
  int segJ;
};

/*! these are pairs to be evaluated. as different moving planes which map onto each other
must contain: variable ids and segment/pixel ids - i have seg2var - right
*/
struct MotionPairs
{
  MotionPairs(int _segI, int _varId, int _segJ0, int _varIdJ0, int _segJ1, int _varIdJ1)
    : segI(_segI), segJ0(_segJ0), segJ1(_segJ1), varId(_varId), varIdJ0(_varIdJ0), varIdJ1(_varIdJ1) {};

  MotionPairs(): varId(-1), segI(-1) {};

  ////////
  /// segment or pixel in first view
  int segI;
  /// segment or pixel in 2nd view if i is assigned 1
  int segJ0;
  /// segment or pixel in 2nd view if i is assigned 1
  int segJ1;
  /// the variable of the segment/pixel (can be -1 == not a variable)
  int varId;
  /// the variable of the segment/pixel (can be -1 == not a variable)
//  int varid2;

  /// varid of match if segI is assigned a 0
  int varIdJ0;

  /// varid of match if segI is assigned a 1
  int varIdJ1;

  /// the assignment of the match (if varid1 = 1 or 0 ) [proposal or current solution]
  //  int assign;
  // no second assignment need - both cases have to be tested if varids!=-1 or just current solution else
};

template<typename Scalar> struct Unary
{
  typedef Math::VectorT<Scalar, 2>  P2;
  Unary(int _segI, P2 _node): var(_segI), node(_node) {};
  Unary(int _segI): var(_segI) {};
  Unary(): var(-1), node(0.,0.) {};
  ////////
  P2 node;
  int var;
};

template<typename Scalar> struct dataElem
{
    typedef Math::Mat3x3T<Scalar>     M3;
    typedef Math::VectorT<Scalar, 4>  P4;
    typedef Math::VectorT<Scalar, 3>  P3;
    typedef Math::VectorT<Scalar, 2>  P2;
    typedef Math::VectorT<int, 4>  P4i;
    typedef Math::VectorT<int, 2>  P2i;

    dataElem( ): firstIdFW(0), firstIdBW(0), varsFW(0), varsBW(0), proposal(-1), w(0),h(0), xtraPen(0.1) {};

    void clear()
    { 
      seg2segFW.clear(); seg2segBW.clear(); freeBW.clear(); freeFW.clear();
      dataFW.clear(); dataBW.clear(); freedataBW.clear(); freedataFW.clear();    
      //depthsFW.clear(); depthsBW.clear(); 
      seg2varFW.clear();seg2varBW.clear();
    };


    inline Scalar getDataG_FW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxFW[0] || qy<boxFW[1] || qx>=boxFW[2] || qy>=boxFW[3] ) 
        return 0;
      else
        return dataFW[qy-boxFW[1]+(qx-boxFW[0])*(boxFW[3]-boxFW[1])];
//        return qx-boxFW[0]+(qy-boxFW[1])*(boxFW[2]-boxFW[0]);
    }
    /// needs this for pixel stuff
    inline Scalar getDataG_BW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxBW[0] || qy<boxBW[1] || qx>=boxBW[2] || qy>=boxBW[3] ) 
        return 0;
      else
        return dataBW[qy-boxBW[1]+(qx-boxBW[0])*(boxBW[3]-boxBW[1])];
//        return qx-boxBW[0]+(qy-boxBW[1])*(boxBW[2]-boxBW[0]);
    }
    inline Scalar getFreeDataG_FW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxFW[0] || qy<boxFW[1] || qx>=boxFW[2] || qy>=boxFW[3] ) 
        return 0;
      else
        return freedataFW[qy-boxFW[1]+(qx-boxFW[0])*(boxFW[3]-boxFW[1])];
//        return qx-boxFW[0]+(qy-boxFW[1])*(boxFW[2]-boxFW[0]);
    }
    /// needs this for pixel stuff
    inline Scalar getFreeDataG_BW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxBW[0] || qy<boxBW[1] || qx>=boxBW[2] || qy>=boxBW[3] ) 
        return 0;
      else
        return freedataBW[qy-boxBW[1]+(qx-boxBW[0])*(boxBW[3]-boxBW[1])];
//        return qx-boxBW[0]+(qy-boxBW[1])*(boxBW[2]-boxBW[0]);
    }
    /// needs this 
    inline int getBoxIdFW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxFW[0] || qy<boxFW[1] || qx>=boxFW[2] || qy>=boxFW[3] ) 
        return -1;
      else
        return qy-boxFW[1]+(qx-boxFW[0])*(boxFW[3]-boxFW[1]);
//        return qx-boxFW[0]+(qy-boxFW[1])*(boxFW[2]-boxFW[0]);
    }
    /// needs this for pixel stuff
    inline int getBoxIdBW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxBW[0] || qy<boxBW[1] || qx>=boxBW[2] || qy>=boxBW[3] ) 
        return -1;
      else
        return qy-boxBW[1]+(qx-boxBW[0])*(boxBW[3]-boxBW[1]);
//        return qx-boxBW[0]+(qy-boxBW[1])*(boxBW[2]-boxBW[0]);
    }
    /// needs this 
    inline int getLocalIdFW( int globalID ) const
    {
      int qx = globalID / h; 
      int qy = globalID % h;

      if ( qx<boxFW[0] || qy<boxFW[1] || qx>=boxFW[2] || qy>=boxFW[3] ) 
        return -1;
      else
        return seg2varFW[qy-boxFW[1]+(qx-boxFW[0])*(boxFW[3]-boxFW[1])];
    }
    /// needs this for pixel stuff
    inline int getLocalIdBW( int globalID ) const 
    {
      int qx = globalID /h;//% w; 
      int qy = globalID %h;/// w;

      if ( qx<boxBW[0] || qy<boxBW[1] || qx>=boxBW[2] || qy>=boxBW[3] ) 
        return -1;
      else
        return seg2varBW[qy-boxBW[1]+(qx-boxBW[0])*(boxBW[3]-boxBW[1])];
    }

    /// from global pixel ids
    P3 getImageCoordG( int globalId ) const
    {
      int qx = globalId /h; 
      int qy = globalId %h;
      return P3( qx+1., qy+1., 1. );
    }
    /// from local pixel ids
    P3 getImageCoordFW( int boxId ) const 
    {
//      int localWidth  = (boxFW[2]-boxFW[0]);
//      int qx = boxId % localWidth;
//      int qy = boxId / localWidth;
      int localHeight = (boxFW[3]-boxFW[1]);
      int qx = boxId /localHeight;
      int qy = boxId %localHeight;
      return P3( boxFW[0]+qx+1, boxFW[1]+qy+1, 1. );
    }
    P3 getImageCoordBW( int boxId ) const 
    {
//      int localWidth = (boxBW[2]-boxBW[0]);
//      int qx = boxId % localWidth;
//      int qy = boxId / localWidth;
      int localHeight = (boxBW[3]-boxBW[1]);
      int qx = boxId /localHeight;
      int qy = boxId %localHeight;
      return P3( boxBW[0]+qx+1, boxBW[1]+qy+1, 1. );
    }

    /*
    void consolidateSeg2VarBWMapping( dataElem<Scalar>& _rhs )
    {
      assert( _rhs.seg2varBW.size() == seg2varBW.size());

      int nElems=0;
      for ( int i=0;i<seg2varBW.size(); i++ )
      {
        if (seg2varBW[i] > -1 || _rhs.seg2varBW[i] > -1 )
        {
          seg2varBW[i]      = nElems++;
          _rhs.seg2varBW[i] = nElems;
        }
      }
      varsBW = nElems;
      _rhs.varsBW = nElems;
    }
    */

    // fix oob by new oob penalty: (could also use xtra oob prediction)
    void updateOobPenaltyFW ( const std::vector<Scalar>& oobPenaltiesFW, const std::vector<int>& areasFW )
    {
      for (int i=0;i<seg2segFW.size();i++)
      {
        int id = seg2segFW[i].first;
        dataFW[i] = freedataFW[i] + (areasFW[id]-freeFW[i]) * (xtraPen+oobPenaltiesFW[id]);// non free variables are oob
      }
    }
    void updateOobPenaltyBW ( const std::vector<Scalar>& oobPenaltiesBW, const std::vector<int>& areasBW )
    {
      for (int i=0;i<seg2segBW.size();i++)
      {
        int id = seg2segBW[i].first;
        dataBW[i] = freedataBW[i] + (areasBW[id]-freeBW[i]) * (xtraPen+oobPenaltiesBW[id]);
      }
    }

    Scalar xtraPen;

    /// forward mapping  
    std::vector< std::pair<int,int> > seg2segFW;
    /// backward mapping  
    std::vector< std::pair<int,int> > seg2segBW;

    std::vector<int> freeFW;
    std::vector<int> freeBW;

    std::vector<Scalar> freedataFW;
    std::vector<Scalar> freedataBW;

    std::vector<Scalar> dataFW;
    std::vector<Scalar> dataBW;

    /// which normal/rotation is stored?
    int proposal;

    /// start segment
//    int seg;

    /// number of variables forward
    int varsFW;
    /// number of variables backward
    int varsBW;

    /// to prevent merging of edges
    int firstIdFW;
    int firstIdBW;
    // not working, H*pix must be dehomogenized
    /// computes the depth in the other view: hvNom|pix = vn' * H*pix
    // P3 hvNom
    
    /// needed for global pixel to box id function
    P4i boxFW;
    P4i boxBW;

    /// width of the image - needed to transfer global pixel id to local box id, and also to seg2seg == pixel2pixel id
    int w;
    /// height of the image
    int h;

    M3  Hom;
    M3 iHom; // from pixel to pixel?
    P3  Nom;// normals in view1, act on pixel in pix coords?
    P3 iNom;// normals in view2

    std::vector<int> seg2varFW;//-1: undef / could also be constructed from scratch each time
    std::vector<int> seg2varBW;//-1: undef / could also be constructed from scratch each time
  };

  template<typename Scalar> struct PatchCenter
  {
    typedef Math::Mat3x3T<Scalar>     M3;
    typedef Math::VectorT<Scalar, 3>  P3;
    typedef Math::VectorT<Scalar, 2>  P2;

    PatchCenter() : px(-1), py (-1), center (P3(0,0,0)), segId (-1)
    {};

    PatchCenter(int px_, int py_, P3 center_, int segId_ ) : px(px_), py (py_), center (center_), segId (segId_)
    {};
  
    bool operator<(const PatchCenter<Scalar>& _rhs) const 
    {
      return (px < _rhs.px) || ( (px == _rhs.px) && (py < _rhs.py) );
    }

    bool operator>(const PatchCenter<Scalar>& _rhs) const 
    {
      return (px > _rhs.px) || ( (px == _rhs.px) && (py > _rhs.py) );
    }
     
    bool operator==(const PatchCenter<Scalar>& _rhs) const 
    {
      return ( (px == _rhs.px) && (py == _rhs.py) );
    }
    //  
    int px; 
    int py; 
    P3 center;
    int segId; // could be removed
  };

#endif