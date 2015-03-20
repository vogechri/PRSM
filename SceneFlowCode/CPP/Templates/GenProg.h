#ifndef _GEN_PROG_H
#define _GEN_PROG_H

#include <limits> // for inf and nan checks


#ifdef max
#undef max
#define REDEFINE_MAX
#endif

#ifdef min
#undef min
#define REDEFINE_MIN
#endif

//== NAMESPACES ===============================================================

namespace GenProg  {

using namespace std;
//== IMPLEMENTATION ===========================================================


template<typename T>
inline bool isAnyInf(T value)
{
  return !((value >= (numeric_limits<T>::min())) && 
           (value <= (numeric_limits<T>::max())));
}

template<typename T>
inline bool isNan(T value)
{
  return (value != value);
};

template<typename T>
inline bool isInf(T value)
{
  return ((numeric_limits<T>::has_infinity) &&
          (value == numeric_limits<T>::infinity()));
}

/// This type maps \c true or \c false to different types.
template <bool b> class Bool2Type { enum { my_bool = b }; };

/// This class generates different types from different \c int 's.
template <int i>  class Int2Type  { enum { my_int = i }; };

/// Handy typedef for Bool2Type<true> classes
//typedef Bool2Type<true> True;

/// Handy typedef for Bool2Type<false> classes
//typedef Bool2Type<false> False;

/// Handy type to use in function head for representing a void type
/// see e.g. Functor.h
struct EmptyType {};

/// Handy type to use in function implementation for representing an iterator
/// without behaviour see e.g. MaskTester.h
struct EmptyIterator 
{
public:

  typedef EmptyType      value_type;

  typedef value_type&    reference;
  typedef value_type*    pointer;

  EmptyIterator(){};
//  EmptyIterator(EmptyIterator e_it){};
  template <typename T> EmptyIterator(T t){};

  ~EmptyIterator(){};

  /// Standard dereferencing operator.
  reference operator*()  const { return EmptyType(); }
  
  /// Standard pointer operator.
  pointer   operator->() const { return NULL; }

  bool operator==(const EmptyIterator) const {return true;}
  bool operator!=(const EmptyIterator) const {return false;}

  EmptyIterator& operator++(){return *this;};
  EmptyIterator& operator--(){return *this;};

  template <typename T>
  EmptyIterator& operator+=(T t){return *this;};

  template <typename T>
  EmptyIterator& operator-=(T t){return *this;};
};


//-----------------------------------------------------------------------------

/// not needed until now
#ifdef _my_false 

//--- Template "if" w/ partial specialization ---------------------------------
#if _PARTIAL_SPECIALIZATION


template <bool condition, class Then, class Else>
struct IF { typedef Then Result; };

/** Template \c IF w/ partial specialization
\code
typedef IF<bool, Then, Else>::Result  ResultType;
\endcode    
*/
template <class Then, class Else>
struct IF<false, Then, Else> { typedef Else Result; };





//--- Template "if" w/o partial specialization --------------------------------
#else


struct SelectThen 
{
  template <class Then, class Else> struct Select {
    typedef Then Result;
  };
};

struct SelectElse
{
  template <class Then, class Else> struct Select {
    typedef Else Result;
  };
};

template <bool condition> struct ChooseSelector {
  typedef SelectThen Result;
};

template <> struct ChooseSelector<false> {
  typedef SelectElse Result;
};


/** Template \c IF w/o partial specialization. Use it like
\code
typedef IF<bool, Then, Else>::Result  ResultType;
\endcode    
*/

template <bool condition, class Then, class Else>
class IF 
{ 
  typedef typename ChooseSelector<condition>::Result  Selector;
public:
  typedef typename Selector::template Select<Then, Else>::Result  Result;
};

#endif
//=============================================================================
#endif  // false
} // namespace GenProg

#ifdef REDEFINE_MAX
#undef REDEFINE_MAX
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifdef REDEFINE_MIN
#undef REDEFINE_MIN
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#endif

//=============================================================================
#endif // GEN_PROG
//=============================================================================