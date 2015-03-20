/*----------------------------------------------------------------------*/
/* modul      : EIGENVAL.H                                              */
/* description: interface                                               */
/*                                                                      */
/*----------------------------------------------------------------------*/

#ifndef __EIGENVALUE_H
#define __EIGENVALUE_H


/*----------------------------------------------------------------------*/
/* includes                                                             */
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* function prototypes                                                  */
/*----------------------------------------------------------------------*/
void  JacobiEigenvector(double *smat, int n, double *eigenvec);
int LowestEigenval   (double *smat, int n, double *eigenvec,
		 									double *lowesteigenvec, double *lowesteigenval);

#endif /* __EIGENVALUE_H */
