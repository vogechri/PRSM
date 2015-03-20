#include "bcucof.cpp"
#ifndef bcuint_FCT
#define bcuint_FCT
 
#define NRANSI
#include "nrutil.h"

// y f(x), y1 derivative in x direction y2 derivative in y, y12 derivative in x and y direction 
// all at the 4 corners of the block : (x1l, x2l), (x1l, x2l)
// 
void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
	double x1u, double x2l, double x2u, double x1, double x2, double *ansy,
	double *ansy1, double *ansy2, double **c)
{
	int i;
	double t,u,d1,d2;//,**c;
//	double *c0, *c1,*c2,*c3;

//	double *c;
//  c  = (double*)  mxMalloc( 16*sizeof(double) );

	//	double c[4][4];
//	c=matrix(1,4,1,4);
//	c  = (double**)  mxMalloc( 4*sizeof(double*) );
//  c0 = (double*)   mxMalloc( 4*sizeof(double) );
//  c1 = (double*)   mxMalloc( 4*sizeof(double) );
//  c2 = (double*)   mxMalloc( 4*sizeof(double) );
//  c3 = (double*)   mxMalloc( 4*sizeof(double) );
//  c[0] = c0; c[1] = c1;  c[2] = c2;  c[3] = c3;


	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2, c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
//  for (i=3;i>=0;i--) {
//    *ansy=t*(*ansy)+((c[i*4+3]*u+c[i*4+2])*u+c[i*4+1])*u+c[i*4+0];
//    *ansy2=t*(*ansy2)+(3.0*c[i*4+3]*u+2.0*c[i*4+2])*u+c[i*4+1];
//    *ansy1=u*(*ansy1)+(3.0*c[3*4+i]*t+2.0*c[2*4+i])*t+c[1*4+i];
//	   *ansy=t*(*ansy)+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
//	   *ansy2=t*(*ansy2)+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
//	   *ansy1=u*(*ansy1)+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];

    // matrix notation indices are weird
	  *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];

    // p(x,y) = sum_I sum_j (a_ij * x^i + y^j)
    //        = t^3 * (a_30 + a_31 * u + a_32 * u^2 + a_33 * u^3+ t^2 * (...) + ...

  }
	*ansy1 /= d1;
	*ansy2 /= d2;

//	mxFree(c0);  mxFree(c1); mxFree(c2); mxFree(c3);

//	mxFree(c);

//	free_matrix(c,1,4,1,4);
};

void bcuint_noDeriv(double y[], double y1[], double y2[], double y12[], double x1l,
	double x1u, double x2l, double x2u, double x1, double x2, double *ansy, double **c)
{
	int i;
	double t,u,d1,d2;//,**c;
//	double *c0, *c1,*c2,*c3;

//	double *c;
//  c  = (double*)  mxMalloc( 16*sizeof(double) );

	//	double c[4][4];
//	c=matrix(1,4,1,4);
//	c  = (double**)  mxMalloc( 4*sizeof(double*) );
//  c0 = (double*)   mxMalloc( 4*sizeof(double) );
//  c1 = (double*)   mxMalloc( 4*sizeof(double) );
//  c2 = (double*)   mxMalloc( 4*sizeof(double) );
//  c3 = (double*)   mxMalloc( 4*sizeof(double) );
//  c[0] = c0; c[1] = c1;  c[2] = c2;  c[3] = c3;


	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2, c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=0.0;
	for (i=4;i>=1;i--) {
//  for (i=3;i>=0;i--) {
//    *ansy=t*(*ansy)+((c[i*4+3]*u+c[i*4+2])*u+c[i*4+1])*u+c[i*4+0];
//    *ansy2=t*(*ansy2)+(3.0*c[i*4+3]*u+2.0*c[i*4+2])*u+c[i*4+1];
//    *ansy1=u*(*ansy1)+(3.0*c[3*4+i]*t+2.0*c[2*4+i])*t+c[1*4+i];
//	   *ansy=t*(*ansy)+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
//	   *ansy2=t*(*ansy2)+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
//	   *ansy1=u*(*ansy1)+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];

    // matrix notation indices are weird
	  *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    // p(x,y) = sum_I sum_j (a_ij * x^i + y^j)
    //        = t^3 * (a_30 + a_31 * u + a_32 * u^2 + a_33 * u^3+ t^2 * (...) + ...
  }

//	mxFree(c0);  mxFree(c1); mxFree(c2); mxFree(c3);
//	mxFree(c);

//	free_matrix(c,1,4,1,4);
};

void bcuint(float y[], float y1[], float y2[], float y12[], float x1l,
	float x1u, float x2l, float x2u, float x1, float x2, float *ansy,
	float *ansy1, float *ansy2, float **c)
{
	int i;
	float t,u,d1,d2;

	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2, c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0f;
	for (i=4;i>=1;i--) {

    // matrix notation indices are weird
	  *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];

  }
	*ansy1 /= d1;
	*ansy2 /= d2;
};

void bcuint_noDeriv(float y[], float y1[], float y2[], float y12[], float x1l,
	float x1u, float x2l, float x2u, float x1, float x2, float *ansy, float **c)
{
	int i;
	float t,u,d1,d2;

  d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2, c);
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=0.0f;
	for (i=4;i>=1;i--) {

    // matrix notation indices are weird
	  *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
  }
};
#undef NRANSI

#endif
