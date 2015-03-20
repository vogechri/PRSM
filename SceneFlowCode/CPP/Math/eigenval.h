/*----------------------------------------------------------------------*/
/* modul      : EIGENVAL.H                                              */
/* description: interface                                               */
/*                                                                      */
/* author     : Peter J. Neugebauer                                     */
/* date       : 15.06.1993                                              */
/*----------------------------------------------------------------------*/

#ifndef __EIGENVAL_H
#define __EIGENVAL_H


/*----------------------------------------------------------------------*/
/* includes                                                             */
/*----------------------------------------------------------------------*/
#include "basetype.h"


/*----------------------------------------------------------------------*/
/* function prototypes                                                  */
/*----------------------------------------------------------------------*/
void  JacobiEigenvector(ExtendedReal *smat, sWord n, ExtendedReal *eigenvec);
sWord LowestEigenval   (ExtendedReal *smat, sWord n, ExtendedReal *eigenvec,
 		        ExtendedReal *lowesteigenvec,
		        ExtendedReal *lowesteigenval);


#endif /* __EIGENVAL_H */
