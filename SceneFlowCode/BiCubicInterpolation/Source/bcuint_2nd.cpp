#include "bcucof.cpp"
#ifndef bcuint_2nd_FCT
#define bcuint_2nd_FCT
 
#define NRANSI
#include "nrutil.h"

// y f(x), y1 derivative in x direction y2 derivative in y, y12 derivative in x and y direction 
// all at the 4 corners of the block : (x1l, x2l), (x1l, x2l)
// 
void bcuint_2nd(double y[], double y1[], double y2[], double y12[], double x1l,
	double x1u, double x2l, double x2u, double x1, double x2, double *ansy,
	double *ansy1, double *ansy2, double *ansy11, double *ansy22, double *ansy12, double **c)
{
	int i;
	double t,u,d1,d2;

	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2, c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {

    // matrix notation indices are weird
    // p(x,y) = sum_I sum_j (a_ij * x^i + y^j)
    //        = t^3 * (a_30 + a_31 * u + a_32 * u^2 + a_33 * u^3+ t^2 * (...) + ...
	  *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];

    // 2nd derivatices:
    *ansy11 = u*(*ansy11)+6.0*c[4][i]*t+2.0*c[3][i];
    *ansy22 = t*(*ansy22)+6.0*c[i][4]*u+2.0*c[i][3];
    if (i > 1)
      *ansy12 = u*(*ansy12) + ((double)(i-1))*(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
	*ansy1 /= d1;
	*ansy2 /= d2;
	*ansy11 /= d1*d1;
	*ansy22 /= d2*d2;
	*ansy12 /= d1*d2;

//	mxFree(c0);  mxFree(c1); mxFree(c2); mxFree(c3);

//	mxFree(c);

//	free_matrix(c,1,4,1,4);
};
#undef NRANSI

#endif
