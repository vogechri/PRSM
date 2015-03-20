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

#ifndef _SEGMENT_HEAP_HELPER_H
#define _SEGMENT_HEAP_HELPER_H

#include "Math/VectorT.h"
#include <limits>

#ifdef max
#undef max
#define REDIFINE_MAX
#endif


/*!
* \class HeapEntryD
* \brief The entry type to decimate an segment image
*
* \author C.Vogel
*/
class HeapEntryD
{

public:

  HeapEntryD(int pId) : _pId(pId)
  {};

  HeapEntryD(const HeapEntryD &he) : _pId(he._pId)
  {};

  HeapEntryD() : _pId(-1)
  {};

  ~HeapEntryD()
  {};

  int pId() const
  {
    return _pId;
  }

  void setP( int id ) 
  {
    _pId = id;
  }

  HeapEntryD& operator=(const HeapEntryD &he)
  {
    _pId = he._pId;
    return *this;
  }

private:

  // index in the store 
  int _pId;
};

// -----------------------------------------------------------------------

/*!
* \class HeapStoreEntryT
* \brief The type stored for each potential segment to collapse.
*
* \author C.Vogel
*/
template <typename ScalarT>
class HeapStoreEntryT
{
  typedef ScalarT Scalar;

public:

  HeapStoreEntryT(Scalar cost, int startSeg) 
    : _startSeg(startSeg), _cost(cost), _heapId (-1)
  {};

  HeapStoreEntryT() 
    : _startSeg(-1), _cost(std::numeric_limits<Scalar>::max()), _heapId (-1)
  {};

  ~HeapStoreEntryT()
  {};

  Scalar cost() 
  {
    return _cost;
  }

  int startSeg()
  {
    return _startSeg;
  }

  int heapId() 
  {
    return _heapId;
  }

  void setCost( Scalar cost ) 
  {
    _cost = cost;
  }

  void setStartSeg( int startSeg ) 
  {
    _startSeg= startSeg;
  }

  void setHeapId( int id ) 
  {
    _heapId = id;
  }


  HeapStoreEntryT& operator=(const HeapStoreEntryT &he)
  {
    _cost = he._cost;
    _startSeg = he._startSeg;
    _heapId = he._heapId;
    return *this;
  }

  void update( HeapStoreEntryT &he)
  {
    _cost = he._cost;
    _startSeg = he._startSeg;
  }

private:

  Scalar _cost;
  int _startSeg;
  int _heapId;
};


//------------------------------------------------------------------------
/*!
* \class HeapInterfaceD
* \brief The Interface class to handle the sorting in the heap.
*
* \author C.Vogel
*/
template <class HeapEntryT, class HeapStoreT>
struct HeapInterfaceD
{

  HeapInterfaceD(HeapStoreT *entryStore) 
  {
    _entryStore = entryStore;
  };

  ~HeapInterfaceD() {};

  /// Comparison of two HeapEntry's: strict less
  bool less(const HeapEntryT& _e1, const HeapEntryT& _e2)
  {
    return ((*_entryStore)[_e1.pId()].cost() < (*_entryStore)[_e2.pId()].cost());
  }

  /// Comparison of two HeapEntry's: strict greater
  bool greater( const HeapEntryT& _e1, const HeapEntryT& _e2)
  {
    return ((*_entryStore)[_e1.pId()].cost() > (*_entryStore)[_e2.pId()].cost());
  }

  /// Get the heap position of HeapEntry _e
  int get_heap_position(const HeapEntryT& _e)
  {
    return (*_entryStore)[_e.pId()].heapId();
  }

  /// Set the heap position of HeapEntry _e
  void set_heap_position(const HeapEntryT& _e, const int _i)
  {
    (*_entryStore)[_e.pId()].setHeapId(_i);
  }


private:
  HeapStoreT *_entryStore;
};
////////////////////////////////////////////////////////////////////////////

#ifdef REDIFINE_MAX
#undef REDIFINE_MAX
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

/////////////////////////////////////////////////
/////////////////////////////////////////////////
#endif // _SEGMENT_HEAP_HELPER_H
