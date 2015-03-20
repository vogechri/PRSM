/*----------------------------------------------------------------------*/
/* modul      : basetype.h                                              */
/* description: interface                                               */
/*                                                                      */
/* author     : Siemens AG, Bensheim                                    */
/* date       : unknown                                                 */
/* changed by : Alireza Esmailpour                                      */
/* date       : 04.11.1994                                              */
/*                                                                      */
/* changed by : Matthias Block                                          */
/* date       : 21.12.1999                                              */
/* change     : Zusaetzliche Konstante "NO_MEMORY" fuer die             */
/*              Speicherverwaltung der Previews eingefuegt              */
/*----------------------------------------------------------------------*/


#ifndef _BASETYPE_H_
#define _BASETYPE_H_


/*----------------------------------------------------------------------*/
/* common defines used by ordinary C-programmers                        */
/*----------------------------------------------------------------------*/
#ifndef SUCCESS
#define SUCCESS         (0)
#endif

#ifndef FAILURE
#define FAILURE         (-1)
#endif

#ifndef FOREVER
#define FOREVER         for(;;)
#endif

#ifndef NIL
#define NIL             (-1)
#endif


#ifndef NO_MEMORY
#define NO_MEMORY       (-2)
#endif

/*----------------------------------------------------------------------*/
/* defines                                                              */
/*----------------------------------------------------------------------*/
#define SHL(j,c) ((j)<<(c))
#define SHR(j,c) ((j)>>(c))

#define ADDR(var) (Pointer) &(var)
#define OR ||
#define AND &&
#define XOR ^
#define MOD %
#define DIV /
#define NOT(arg)  (! (arg))
#define NEG(arg)  (~ (arg))

/*----------------------------------------------------------------------*/
/* datatypes                                                            */
/*----------------------------------------------------------------------*/
/*
typedef char*           String;
*/
typedef unsigned char   uByte;
/*
typedef signed char     sByte;
*/
typedef unsigned short  uWord;
typedef short           sWord;
typedef unsigned long   uDoubleWord;
typedef long            sDoubleWord;
typedef double          ExtendedReal;

typedef uByte*          Pointer;
typedef Pointer*        Handle;


/*----------------------------------------------------------------------*/
/* Constants                                                            */
/*----------------------------------------------------------------------*/
#define kMaxuByte         0xFF
#define kMaxsByte         0x7F
#define kMaxuWord         0xFFFF
#define kMaxsWord         0x7FFF
#define kMaxuDoubleWord   0xFFFFFFFFUL
#define kMaxsDoubleWord   0x7FFFFFFFL
#define kMaxReal          3.37E+38
#define kMinReal          8.43E-37
#define kEpsReal          1.19209290E-07F
#define kMaxExtendedReal  1.797693E+308
#define kMinExtendedReal  2.225074E-308
#define kEpsExtendedReal  2.2204460492503131E-16


#endif /* _BASETYPE_H_ */
