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

#ifndef _BINARY_HEAP_HELPER_H
#define _BINARY_HEAP_HELPER_H

#include "Math/VectorT.h"
#include <limits>

#ifdef max
#undef max
#define REDIFINE_MAX
#endif


/*!
* \class HeapEntryB
* \brief The entry type to decimate an segment image
*
* \author C.Vogel
*/
class HeapEntryB
{

public:

  HeapEntryB(long pId) : _pId(pId)
  {};

  HeapEntryB(const HeapEntryB &he) : _pId(he._pId)
  {};

  HeapEntryB() : _pId(-1)
  {};

  ~HeapEntryB()
  {};

  long pId() const
  {
    return _pId;
  }

  void setP( long id ) 
  {
    _pId = id;
  }

  HeapEntryB& operator=(const HeapEntryB &he)
  {
    _pId = he._pId;
    return *this;
  }

private:

  // index in the store 
  long _pId;
};

// -----------------------------------------------------------------------

/*!
* \class HeapStoreEntryTb
* \brief The type stored for each potential variable to be substituted.
*
* \author C.Vogel
*/
template <typename ScalarT>
class HeapStoreEntryTb
{
  typedef ScalarT Scalar;

public:

  HeapStoreEntryTb(Scalar cost, long vid) 
    : _vid(vid), _cost(cost), _heapId (-1)
  {};

  HeapStoreEntryTb() 
    : _vid(-1), _cost(std::numeric_limits<Scalar>::min()), _heapId (-1)
  {};

  ~HeapStoreEntryTb()
  {};

  Scalar cost() 
  {
    return _cost;
  }

  long literalNr()
  {
    return _vid;
  }

  long heapId() 
  {
    return _heapId;
  }

  void setCost( Scalar cost ) 
  {
    _cost = cost;
  }

  void setLiteralNr( long vid ) 
  {
    _vid = vid;
  }

  void setHeapId( long id ) 
  {
    _heapId = id;
  }


  HeapStoreEntryTb& operator=(const HeapStoreEntryTb &he)
  {
    _cost = he._cost;
    _vid  = he._vid;
    _heapId = he._heapId;
    return *this;
  }

  void update( HeapStoreEntryTb &he)
  {
    _cost = he._cost;
//    _heapId = he._heapId; // always like this ! activating is dead wrong (set to -1->bang)
    _vid = he._vid;
  }

private:

  Scalar _cost;
  long   _vid;
  long   _heapId;
};


//------------------------------------------------------------------------
/*!
* \class HeapInterfaceDb
* \brief The Interface class to handle the sorting in the heap.
*
* \author C.Vogel
*/
template <class HeapEntryT, class HeapStoreT>
struct HeapInterfaceDb
{

  HeapInterfaceDb(HeapStoreT *entryStore) 
  {
    _entryStore = entryStore;
  };

  ~HeapInterfaceDb() {};

  /// Comparison of two HeapEntry's: strict less
  bool less(const HeapEntryT& _e1, const HeapEntryT& _e2)
  {
    return ((*_entryStore)[_e1.pId()].cost() > (*_entryStore)[_e2.pId()].cost());
  }

  /// Comparison of two HeapEntry's: strict greater
  bool greater( const HeapEntryT& _e1, const HeapEntryT& _e2)
  {
    return ((*_entryStore)[_e1.pId()].cost() < (*_entryStore)[_e2.pId()].cost());
  }

  /// Get the heap position of HeapEntry _e
  long get_heap_position(const HeapEntryT& _e)
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
#endif // _BINARY_HEAP_HELPER_H