/*----------------------------------------------------------------------*/
/* modul      : EIGENVAL.C						*/
/*                                                                      */
/* description: Some utility-functions for computing eigenvalues and 	*/
/*		eigenvectors.						*/
/*                                                                      */
/* conventions: As a common convention in C a matrix is considered as a */
/*		set of row vectors. MATRIX[i][j] specifies the i-th row */
/*		and the j-th column, starting with row and column	*/
/*		number 0.						*/
/*		    If for indices i<j all elements are 0, then the	*/
/*		matrix is called an upper right triangular matrix. The	*/
/*		values in the lower left part then are generally not	*/
/*		considered by the routines of this modul. In case of a	*/
/*		symmetric matrix the elements in the lower left part	*/
/*		can be easily obtained from the values in the upper	*/
/*		right part. So the routines do not consider the values	*/
/*		in the lower left part of a given symmetric matrix.	*/
/*                                                                      */
/* references : Schwarz, Rutishauser, Stiefel:                          */
/*              Numerik symmetrischer Matrizen.                         */
/*              B.G.Teubner-Verlag, Stuttgart 1968                      */
/*                                                                      */
/*              Eduard Stiefel:                                         */
/*              Einfuehrung in die numerische Mathematik.               */
/*              B.G.Teubner-Verlag, Stuttgart 1969                      */
/*                                                                      */
/*              I.Bronstein:                                            */
/*              Taschenbuch der Mathematik (21./22.Auflage)             */
/*              Verlag Harri Deutsch, Thun und Frankfurt(Main) 1981     */
/*                                                                      */
/* author     : P. J. Neugebauer                                        */
/* date       : 13.06.91                                                */
/* updated    : 15.12.93						*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* includes                                                             */
/*----------------------------------------------------------------------*/
#include <math.h>
#include "eigenval.h"


/*----------------------------------------------------------------------*/
/* implementation                                                       */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* name       : sqr							*/
/*                                                                      */
/* parameters : ExtendedReal x			the argument		*/
/*                                                                      */
/* description: Computes x*x						*/
/*                                                                      */
/* returnvalue: the product x*x						*/
/*----------------------------------------------------------------------*/

static ExtendedReal sqr(ExtendedReal x)
{
  return x*x;
}


/*----------------------------------------------------------------------*/
/* name       : JacobiRotation                                          */
/*                                                                      */
/* parameters : ExtendedReal *smat	symmetrical matrix smat[n][n]	*/
/*		sWord n			order n of the matrix		*/
/*		ExtendedReal cos_phi	cos of the rotation angle	*/
/*		ExtendedReal sin_phi	sin of the rotation angle	*/
/*		sWord ind1, ind2	rotation index pair (ind1<ind2)	*/
/*                                                                      */
/* description: Performs an elementary jacobi rotation for a		*/
/*		symmetrical matrix.					*/
/*                                                                      */
/* returnvalue: void                                                    */
/*                                                                      */
/* author     : P. J. Neugebauer                                        */
/* date       : 13.06.1991                                              */
/*----------------------------------------------------------------------*/

static void JacobiRotation(ExtendedReal *smat, sWord n, ExtendedReal cos_phi,
                           ExtendedReal sin_phi, sWord ind1, sWord ind2)
{
  ExtendedReal g, h;
  ExtendedReal sqr_c, sqr_s;
  ExtendedReal prd_cs, prd2cs;
  sWord j;
  ExtendedReal *act1, *act2, *act3;


  /* act1 = &smat[ind1][ind1] */
  /* act2 = &smat[ind1][ind2] */
  /* act3 = &smat[ind2][ind2] */
  act1 = smat + (sDoubleWord)(n+1) * ind1;
  act2 = smat + (sDoubleWord)n * ind1 + ind2;
  act3 = smat + (sDoubleWord)(n+1) * ind2;

  h = (sqr_c = sqr(cos_phi)) * *act1
    - (prd2cs = 2.0 * (prd_cs = cos_phi * sin_phi)) * *act2
    + (sqr_s = sqr(sin_phi)) * *act3;
  g = sqr_s * *act1 + prd2cs * *act2 + sqr_c * *act3;
  *act2 = prd_cs * (*act1 - *act3) + (sqr_c - sqr_s) * *act2;
  *act1 = h;
  *act3 = g;

  act1 = smat;
  for (j=0; j<ind1; j++, act1 += n) {
    h          = cos_phi * act1[ind1] - sin_phi * act1[ind2];
    act1[ind2] = sin_phi * act1[ind1] + cos_phi * act1[ind2];
    act1[ind1] = h;
  }

  act2 = act1 + n;
  for (j=ind1 + 1; j<ind2; j++, act2 += n) {
    h          = cos_phi * act1[j] - sin_phi * act2[ind2];
    act2[ind2] = sin_phi * act1[j] + cos_phi * act2[ind2];
    act1[j] = h;
  }

  for (j=ind2 + 1; j<n; j++) {
    h       = cos_phi * act1[j] - sin_phi * act2[j];
    act2[j] = sin_phi * act1[j] + cos_phi * act2[j];
    act1[j] = h;
  }
}


/*----------------------------------------------------------------------*/
/* name       : JacobiEigenvector					*/
/*                                                                      */
/* parameters : ExtendedReal *smat	symmetric matrix smat[n][n]	*/
/*		sWord n			order n of the matrix		*/
/*		ExtendedReal *eigenvec	matrix[n][n] with eigenvectors	*/
/*                                                                      */
/* description: Computes all eigenvalues and eigenvectors of a		*/
/*		symmetric matrix with the special cyclic jacobi		*/
/*		algorithmen.						*/
/*		After the computation the matrix smat contains the	*/
/*		eigenvalues on the diagonal. The matrix eigenvec	*/
/*		contains the eigenvectors as column-vectors.		*/
/* note       : The cyclic jacobi algorithm works even with matrices    */
/*              which are not regular. It is assured that the returned  */
/*              eigenvectors are perpendicular to each other and that   */
/*              they are all of length 1.                               */
/*                                                                      */
/* returnvalue: void							*/
/*                                                                      */
/* author     : P. J. Neugebauer                                        */
/* date       : 13.06.1991                                              */
/*----------------------------------------------------------------------*/

void JacobiEigenvector(ExtendedReal *smat, sWord n, ExtendedReal *eigenvec) 
{
  sWord n_minus1 = n - 1;
  sWord i, j, k;                      /* loop and index-variables       */
  ExtendedReal ss;                    /* compare value for termination  */
  ExtendedReal theta, t, h;
  ExtendedReal cos_phi, sin_phi;
  ExtendedReal *act1, *act2, *act3;


  /* initialize the matrix with the initial eigenvectors */
  act1 = eigenvec;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++)
      *act1++ = 0.0;
    *act1++ = 1.0;
    for (j=i+1; j<n; j++)
      *act1++ = 0.0;
  }

  FOREVER {
    /* compute and test the termination criterium */
    act1 = smat;
    ss = 0.0;
    for (i=0; i<n_minus1; i++, act1 += n)
      for (j=i+1; j<n; j++)
        ss += sqr(act1[j]);
    if (ss < 1.0E3 * kEpsExtendedReal / 2.0)
      break;

    /* run cyclic through all matrixelements in the upper right part */
    for (i=0; i<n_minus1; i++) {
      for (j=i+1; j<n; j++) {
        act1 = smat + (sDoubleWord)(n+1) * i;
        act2 = smat + (sDoubleWord)n * i + j;
        act3 = smat + (sDoubleWord)(n+1) * j;

        /* rotate matrix */
        if (*act2 != 0.0) {
          theta = 0.5 * (*act3 - *act1) / *act2;
          if (theta > 0.0)
            t = 1.0 / (theta + sqrt(1.0 + sqr(theta)));
          else if (theta < 0.0)
            t = 1.0 / (theta - sqrt(1.0 + sqr(theta)));
          else
            t = 1.0;
          cos_phi = 1.0 / sqrt(1.0 + sqr(t));
          sin_phi = cos_phi * t;
          JacobiRotation(smat, n, cos_phi, sin_phi, i, j);

          /* accumulate eigenvector matrix */
          act1 = eigenvec;
          for (k=0; k<n; k++, act1 += n) {
            h       = cos_phi * act1[i] - sin_phi * act1[j];
            act1[j] = sin_phi * act1[i] + cos_phi * act1[j];
            act1[i] = h;
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/* name       : LowestEigenval						*/
/*                                                                      */
/* parameters : ExtendedReal *smat      symmetric matrix smat[n][n]     */
/*              sWord n                 order n of the matrix           */
/*              ExtendedReal *eigenvec, matrix[n][n] with eigenvectors  */
/*		       *lowesteigenvec, eigenvector[n] belonging to the	*/
/*					lowest eigenvalue		*/
/*		       *lowesteigenval  lowest eigenvalue		*/
/*                                                                      */
/* description: Computes all eigenvalues and eigenvectors of a          */
/*              symmetric matrix with the special cyclic jacobi         */
/*              algorithmen.                                            */
/*              After the computation the matrix smat contains the      */
/*              eigenvalues on the diagonal. The matrix eigenvec        */
/*              contains the eigenvectors as column-vectors.            */
/*		In addition, the vector lowesteigenvec[n] is a copy	*/
/*		of the eigenvector belonging to the lowest eigenvalue.	*/
/*		The lowest eigenvalue will be returned in the argument	*/
/*		lowesteigenval.						*/
/* note       : The cyclic jacobi algorithm works even with matrices	*/
/*		which are not regular. It is assured that the returned	*/
/*		eigenvectors are perpendicular to each other and that	*/
/*		they are all of length 1.				*/
/*                                                                      */
/* returnvalue: SUCCESS: all well done					*/
/*		FAILURE: the input matrix is not regular		*/
/*                                                                      */
/* author     : P. J. Neugebauer                                        */
/* date       : 15.12.93                                                */
/*----------------------------------------------------------------------*/

sWord LowestEigenval(ExtendedReal *smat, sWord n,  ExtendedReal *eigenvec,
                     ExtendedReal *lowesteigenvec, ExtendedReal *lowesteigenval)
{
  sWord i,j, index;			/* loop and index variables	*/
  sWord nzero;				/* number of eigenvalues which	*/
					/* are zero			*/
  ExtendedReal *act;			/* row-pointer			*/
  ExtendedReal epsilon;			/* x < epsilon => x is zero	*/


  /* break on abnormal cases */
  if (n < 1)
    return SUCCESS;

  /* compute eigenvectors and eigenvalues */
  JacobiEigenvector(smat, n, eigenvec);

  /* search minimum and count number of eigenvalues which are 0 */
  epsilon = sqrt(1.0E3 * kEpsExtendedReal / 2.0);
  nzero = 0;
  act = smat;
  *lowesteigenval = kMaxExtendedReal;
  for (i=0; i<n; i++, act += n) {
    if (act[i] < *lowesteigenval) {
      index = i;
      *lowesteigenval = act[i];
    }
    if (act[i] < epsilon)
      nzero++;
  }

  /* reduce numerical problems */
  if (*lowesteigenval < 0)
    *lowesteigenval = 0.0;

  /* copy eigenvector to the output vector */
  act = eigenvec;
  for (j=0; j<n; j++, act += n)
    lowesteigenvec[j] = act[index];

  /* is the matrix regular? */
  if (nzero <= 1)
    return SUCCESS;
  else
    return FAILURE;
}
