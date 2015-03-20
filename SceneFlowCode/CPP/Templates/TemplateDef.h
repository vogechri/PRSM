#ifndef _TEMPLATE_DEF_H
#define _TEMPLATE_DEF_H

// --------------------------------------------------------- MS Visual C++ ----
// Compiler _MSC_VER
// .NET 2002 1300 
// .NET 2003 1310
#if defined(_MSC_VER)
//#  define CC_MSVC
#  if (_MSC_VER == 1300)
#    define OUT_OF_CLASS_TEMPLATE       0
#    define PARTIAL_SPECIALIZATION      0
#    define INCLUDE_TEMPLATES           1
#  elif (_MSC_VER > 1300)
#    define OUT_OF_CLASS_TEMPLATE       1
#    define PARTIAL_SPECIALIZATION      1
#    define INCLUDE_TEMPLATES           1
#  else
#    error "Version 7 (.NET 2002) or higher of the MS VC++ is required!"
#  endif
/*
#  if defined(__cplusplus) && !defined(_CPPRTTI)
#    error "Enable Runtime Type Information (Compiler Option /GR)!"
#  endif
#  if !defined(_USE_MATH_DEFINES)
#    error "You have to define _USE_MATH_DEFINES in the compiler settings!"
#  endif
*/
#else
#    define OUT_OF_CLASS_TEMPLATE       1
#    define PARTIAL_SPECIALIZATION      1
#    define INCLUDE_TEMPLATES           1
#endif

//=============================================================================
#endif // ifdef _TEMPLATE_DEF_H
//=============================================================================
