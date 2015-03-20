#ifndef _SEGGROW_HEAP_HELPER_H
#define _SEGGROW_HEAP_HELPER_H

#include "Math/VectorT.h"
#include <limits>

#ifdef max
#undef max
#define REDIFINE_MAX
#endif


/*!
* \class HeapEntryG
* \brief The entry type to decimate an segment image
*
* \author C.Vogel
*/
class HeapEntryG
{

public:

  HeapEntryG(long pId) : _pId(pId)
  {};

  HeapEntryG(const HeapEntryG &he) : _pId(he._pId)
  {};

  HeapEntryG() : _pId(-1)
  {};

  ~HeapEntryG()
  {};

  long pId() const
  {
    return _pId;
  }

  void setP( long id ) 
  {
    _pId = id;
  }

  HeapEntryG& operator=(const HeapEntryG &he)
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
class HeapStoreEntryTg
{
  typedef ScalarT Scalar;

public:

  HeapStoreEntryTg(Scalar cost, long vid, long pid, long sidold) 
    : _vid(vid), _cost(cost), _pid(pid), _heapId (-1), _sidold( sidold )
  {};

  HeapStoreEntryTg() 
    : _vid(-1), _cost(std::numeric_limits<Scalar>::max()), _heapId (-1), _sidold( -1 )
  {};

  ~HeapStoreEntryTg()
  {};

  Scalar cost() 
  {
    return _cost;
  }

  long segNr()
  {
    return _vid;
  }

  long segNr_left()
  {
    return _sidold;
  }

  long heapId() 
  {
    return _heapId;
  }

  long pixNr( ) 
  {
    return _pid;
  }

  void setCost( Scalar cost ) 
  {
    _cost = cost;
  }

  void setSegNr( long vid ) 
  {
    _vid = vid;
  }

  void setSegNrLeft( long vid ) 
  {
    _sidold = vid;
  }

  void setPixNr( long pid ) 
  {
    _pid = pid;
  }

  void setHeapId( long id ) 
  {
    _heapId = id;
  }


  HeapStoreEntryTg& operator=(const HeapStoreEntryTg &he)
  {
    _cost = he._cost;
    _vid  = he._vid;
    _pid  = he._pid;
    _heapId = he._heapId;
    return *this;
  }

  void update( HeapStoreEntryTg &he)
  {
    _cost = he._cost;
//    _heapId = he._heapId; // always like this ! activating is dead wrong (set to -1->bang)
    _vid = he._vid;
  }

private:

  Scalar _cost;
  long   _vid;
  long   _pid;
  long   _sidold;
  long   _heapId;
};


//------------------------------------------------------------------------
/*!
* \class HeapInterfaceDg
* \brief The Interface class to handle the sorting in the heap.
*
* \author C.Vogel
*/
template <class HeapEntryT, class HeapStoreT>
struct HeapInterfaceDg
{

  HeapInterfaceDg(HeapStoreT *entryStore) 
  {
    _entryStore = entryStore;
  };

  ~HeapInterfaceDg() {};

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
#endif // _SEGGROW_HEAP_HELPER_H