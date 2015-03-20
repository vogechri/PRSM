#ifndef _BICUBIC_IP_H
#define _BICUBIC_IP_H

#include <math.h>
#include "bcuint.cpp"
#include<vector>
#include <assert.h>
//#include "bcuint_2nd.cpp"
//#include "mex.h"

//#define _NO_OPENMP

#ifndef _NO_OPENMP
#include <omp.h>
#endif

#define NRANSI
#include "nrutil.h"
#undef NRANSI

#undef DebugOut

template <typename S_>
class BiCubicIp
{
  public:

    typedef S_ Scalar;

  BiCubicIp(int _N, int _M, Scalar*& _Img ) :Idx(NULL), Idy(NULL), Idxy(NULL), N(_N), M(_M), I(_Img)
  {generateMem();};

  ~BiCubicIp()
  {
    freeMem();    
    deleteMatrices();
  };

/// compute interpolation and store in the respective containers
void interpolate( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF, Scalar *&outputdFX, Scalar *&outputdFY )
{
//  double **cc;
//  cc=matrix(1,4,1,4);


#pragma omp parallel for schedule (static)
for(int id=0;id<elements;id++)
{

#ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif

  int tNum(0); 
#ifndef _NO_OPENMP
  tNum = omp_get_thread_num( ); // the current thread, defines where to store the temp per thread solution
#endif

    int intx, inty;
    Scalar x,y;

    // c++ style indices 
    Scalar idx = X[id] -Scalar(1.0);
    Scalar idy = Y[id] -Scalar(1.0);

    if ( (idx >= Scalar(N-1)) )
      idx = idx - Scalar(0.000001);

    if ( (idy >= Scalar(M-1)) )
      idy = idy - Scalar(0.000001);

    if ( (idx < Scalar(0)) )
      idx = idx + Scalar(0.000001);

    if ( (idy < Scalar(0)) )
      idy = idy + Scalar(0.000001);

    x = floor( idx );
    y = floor( idy );

    intx = int(x);
    inty = int(y);

    if ( (intx >= N-1) || (intx < 0) || (inty >= M-1) || (inty < 0) )
    {
      if (intx < 0)   intx = 0;
      if (intx > N-1) intx = N-1;
      if (inty < 0)   inty = 0;
      if (inty > M-1) inty = M-1;
      int myid = intx*M+inty;
      outputF[id]  = I[myid];
      outputdFX[id] = Scalar(0.0);
      outputdFY[id] = Scalar(0.0);
      continue;
    }

#ifdef DebugOut
    mexPrintf("n=%d m=%d : %d %d  %f - %f\n", n, m, intx, inty, X[id], Y[id]);
    int myid = intx*M+inty+1;
    Scalar a  = I[myid];
    Scalar b = Idx[myid];
    Scalar c = Idy[myid];
    Scalar d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty+1, a, b, c, d);
     myid = intx*M+inty+1+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty+1, a, b, c, d);
    myid = intx*M+inty+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty, a, b, c, d);
     myid = intx*M+inty;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty, a, b, c, d);
#endif

    // first get t and u and the 
//    t = idx-x;
//    u = idy-y;

    Scalar f[5], y1[5], y2[5], y12[5];

    int k = 2, add = 1;
    // start at lower left and go counterclockwise
    int idCC = intx*M + inty + 1;

    k = 4;
    f[k]   =    I[idCC];
    y1[k]  =  Idx[idCC];
    y2[k]  =  Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 3;
    idCC   = intx*M + inty + 1 + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 2;
    idCC   = intx*M + inty + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 1;
    idCC   = intx*M + inty;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

/*
    {
    Scalar d1,d2,**c;
    c=matrix(1,4,1,4);

    d1=idx-x;
    d2=idy-y;
    bcucof(f,y1,y2,y12,d1,d2, c);
    for (i=1;i<=4;i++)
      for (j=1;j<=4;j++)
        mexPrintf(" %f ", c[i][j]);
    mexPrintf("\n");
    }
    */

//    bcuint_noDeriv(f, y1, y2, y12, x, x+Scalar(1.0), y, y+Scalar(1.0), idx, idy, &(outputF[id]), cc[tNum] );
    bcuint (f, y1, y2, y12, x, x+Scalar(1.0), y, y+Scalar(1.0), idx, idy, &(outputF[id]), &(outputdFX[id]), &(outputdFY[id]), cc[tNum] );

    //    free_matrix(cc,1,4,1,4);

}
//for (int i=0;i<mThreads;i++)
//  free_matrix(cc[i],1,4,1,4);

}

/// compute interpolation and store in the respective containers
void interpolate_2nd( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF, Scalar *&outputdFX, Scalar *&outputdFY, Scalar *&outputdFXX, Scalar *&outputdFYY, Scalar *&outputdFXY )
{
  #pragma omp parallel for schedule (static)
for(int id=0;id<elements;id++)
{

#ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif

  int tNum(0); 
#ifndef _NO_OPENMP
  tNum = omp_get_thread_num( ); // the current thread, defines where to store the temp per thread solution
#endif

    int intx, inty;
    Scalar x,y;

    // c++ style indices 
    Scalar idx = X[id] -Scalar(1.0);
    Scalar idy = Y[id] -Scalar(1.0);

//    int id  = m*N+n;

    Scalar t,u;

    if ( (idx >= Scalar (N-1)) )
      idx = idx - 0.000001;

    if ( (idy >= Scalar (M-1)) )
      idy = idy - 0.000001;

    if ( (idx < Scalar(0)) )
      idx = idx + 0.000001;

    if ( (idy < Scalar(0)) )
      idy = idy + 0.000001;

    x = floor( idx );
    y = floor( idy );

    intx = int(x);
    inty = int(y);

    if ( (intx >= N-1) || (intx < 0) || (inty >= M-1) || (inty < 0) )
    {
      if (intx < 0)   intx = 0;
      if (intx > N-1) intx = N-1;
      if (inty < 0)   inty = 0;
      if (inty > M-1) inty = M-1;
      int myid = intx*M+inty;
      outputF[id]  = I[myid];
      outputdFX[id] = 0.0;
      outputdFY[id] = 0.0;
      outputdFXX[id] = 0.0;
      outputdFYY[id] = 0.0;
      outputdFXY[id] = 0.0;
      continue;
    }

#ifdef DebugOut
    mexPrintf("n=%d m=%d : %d %d  %f - %f\n", n, m, intx, inty, X[id], Y[id]);
    int myid = intx*M+inty+1;
    double a  = I[myid];
    double b = Idx[myid];
    double c = Idy[myid];
    double d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty+1, a, b, c, d);
     myid = intx*M+inty+1+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty+1, a, b, c, d);
    myid = intx*M+inty+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty, a, b, c, d);
     myid = intx*M+inty;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty, a, b, c, d);
#endif

    // first get t and u and the 
//    t = idx-x;
//    u = idy-y;

    double f[5], y1[5], y2[5], y12[5];

    int k = 2, add = 1;
    // start at lower left and go counterclockwise
    int idCC = intx*M + inty + 1;

    k = 4;
    f[k]   =    I[idCC];
    y1[k]  =  Idx[idCC];
    y2[k]  =  Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 3;
    idCC   = intx*M + inty + 1 + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 2;
    idCC   = intx*M + inty + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 1;
    idCC   = intx*M + inty;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    bcuint_2nd(f, y1, y2, y12, x, x+1.0, y, y+1.0, idx, idy, &(outputF[id]), &(outputdFX[id]), &(outputdFY[id]), &(outputdFXX[id]), &(outputdFYY[id]), &(outputdFXY[id]), cc[tNum] );

  }
}
/// compute interpolation and store in the respective containers
void interpolate( int elements, Scalar*& X, Scalar *&Y, Scalar *&outputF )
{

#pragma omp parallel for schedule (static)
for(int id=0;id<elements;id++)
{
//  Scalar **cc;
//  cc=matrix(1,4,1,4);

#ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif
  int tNum(0); 
#ifndef _NO_OPENMP
    tNum = omp_get_thread_num( ); // the current thread, defines where to store the temp per thread solutuion
#endif
    int intx, inty;
    Scalar x,y;
//    Scalar t,u;

    // c++ style indices 
    Scalar idx = X[id] -Scalar(1.0);
    Scalar idy = Y[id] -Scalar(1.0);

    if ( (idx >= Scalar(N-1)) )
      idx = idx - Scalar(0.000001);

    if ( (idy >= Scalar(M-1)) )
      idy = idy - Scalar(0.000001);

    if ( (idx < Scalar(0)) )
      idx = idx + Scalar(0.000001);

    if ( (idy < Scalar(0)) )
      idy = idy + Scalar(0.000001);

    x = floor( idx );
    y = floor( idy );

    intx = int(x);
    inty = int(y);

    if ( (intx >= N-1) || (intx < 0) || (inty >= M-1) || (inty < 0) )
    {
      if (intx < 0)   intx = 0;
      if (intx > N-1) intx = N-1;
      if (inty < 0)   inty = 0;
      if (inty > M-1) inty = M-1;
      int myid = intx*M+inty;
      outputF[id]  = I[myid];
      continue;
    }

#ifdef DebugOut
    mexPrintf("n=%d m=%d : %d %d  %f - %f\n", n, m, intx, inty, X[id], Y[id]);
    int myid = intx*M+inty+1;
    Scalar a  = I[myid];
    Scalar b = Idx[myid];
    Scalar c = Idy[myid];
    Scalar d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty+1, a, b, c, d);
     myid = intx*M+inty+1+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty+1, a, b, c, d);
    myid = intx*M+inty+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty, a, b, c, d);
     myid = intx*M+inty;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty, a, b, c, d);
#endif

    // first get t and u and the 
//    t = idx-x;
//    u = idy-y;

    Scalar f[5], y1[5], y2[5], y12[5];

    int k = 2, add = 1;
    // start at lower left and go counterclockwise
    int idCC = intx*M + inty + 1;

    k = 4;
    f[k]   =    I[idCC];
    y1[k]  =  Idx[idCC];
    y2[k]  =  Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 3;
    idCC   = intx*M + inty + 1 + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 2;
    idCC   = intx*M + inty + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 1;
    idCC   = intx*M + inty;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

/*
    {
    Scalar d1,d2,**c;
    c=matrix(1,4,1,4);

    d1=idx-x;
    d2=idy-y;
    bcucof(f,y1,y2,y12,d1,d2, c);
    for (i=1;i<=4;i++)
      for (j=1;j<=4;j++)
        mexPrintf(" %f ", c[i][j]);
    mexPrintf("\n");
    }
    */
    bcuint_noDeriv(f, y1, y2, y12, x, x+Scalar(1.0), y, y+Scalar(1.0), idx, idy, &(outputF[id]), cc[tNum] );
}

}


/// compute interpolation and store in the respective containers
void interpolate( std::vector<Scalar>& X, std::vector<Scalar>& Y, std::vector<Scalar>& outputF, std::vector<Scalar>& outputdFX, std::vector<Scalar>& outputdFY )
{
  int elements = X.size();
  assert (outputF.size() == X.size() );

//#pragma omp parallel for schedule (static)
for(int id=0;id<elements;id++)
{

#ifdef DebugOut
  int numPrcs = omp_get_num_procs();
  int th_id = omp_get_thread_num();
  printf("row %d threadid is : %d\n", n, th_id);
#endif

  int tNum(0); 
#ifndef _NO_OPENMP
  tNum = omp_get_thread_num( ); // the current thread, defines where to store the temp per thread solution
#endif

    int intx, inty;
    Scalar x,y;

    // c++ style indices 
    Scalar idx = X[id] -Scalar(1.0);
    Scalar idy = Y[id] -Scalar(1.0);

    if ( (idx >= Scalar(N-1)) )
      idx = idx - Scalar(0.000001);

    if ( (idy >= Scalar(M-1)) )
      idy = idy - Scalar(0.000001);

    if ( (idx < Scalar(0)) )
      idx = idx + Scalar(0.000001);

    if ( (idy < Scalar(0)) )
      idy = idy + Scalar(0.000001);

    x = floor( idx );
    y = floor( idy );

    intx = int(x);
    inty = int(y);

    if ( (intx >= N-1) || (intx < 0) || (inty >= M-1) || (inty < 0) )
    {
      if (intx < 0)   intx = 0;
      if (intx > N-1) intx = N-1;
      if (inty < 0)   inty = 0;
      if (inty > M-1) inty = M-1;
      int myid = intx*M+inty;
      outputF[id]  = I[myid];
      outputdFX[id] = Scalar(0.0);
      outputdFY[id] = Scalar(0.0);
      continue;
    }

#ifdef DebugOut
    mexPrintf("n=%d m=%d : %d %d  %f - %f\n", n, m, intx, inty, X[id], Y[id]);
    int myid = intx*M+inty+1;
    Scalar a  = I[myid];
    Scalar b = Idx[myid];
    Scalar c = Idy[myid];
    Scalar d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty+1, a, b, c, d);
     myid = intx*M+inty+1+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty+1, a, b, c, d);
    myid = intx*M+inty+M;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx+1, inty, a, b, c, d);
     myid = intx*M+inty;
     a  = I[myid];
     b = Idx[myid];
     c = Idy[myid];
     d = Idxy[myid];
    mexPrintf("id=(%d %d) : %f %f %f %f\n", intx, inty, a, b, c, d);
#endif

    Scalar f[5], y1[5], y2[5], y12[5];

    int k = 2, add = 1;
    // start at lower left and go counterclockwise
    int idCC = intx*M + inty + 1;

    k = 4;
    f[k]   =    I[idCC];
    y1[k]  =  Idx[idCC];
    y2[k]  =  Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 3;
    idCC   = intx*M + inty + 1 + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 2;
    idCC   = intx*M + inty + M;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    k = 1;
    idCC   = intx*M + inty;
    f[k]   = I[idCC];
    y1[k]  = Idx[idCC];
    y2[k]  = Idy[idCC];
    y12[k] = Idxy[idCC];

    bcuint (f, y1, y2, y12, x, x+Scalar(1.0), y, y+Scalar(1.0), idx, idy, &(outputF[id]), &(outputdFX[id]), &(outputdFY[id]), cc[tNum] );
}
}

void approximateImageGradients()
{
  freeMem();
  Idx  = (Scalar*) malloc( N*M*sizeof(Scalar) );
  Idy  = (Scalar*) malloc( N*M*sizeof(Scalar) );
  Idxy = (Scalar*) malloc( N*M*sizeof(Scalar) );

/////////////////////////////////////////////////////////
// Start OF APPROXIMATING IMAGE GRADIENTS AS PRE-STEP (those are used further below for bi-cubic stuff)
/////////////////////////////////////////////////////////
#pragma omp for //num_threads(omp_get_num_procs())
for(int n=0;n<N;n++)
{
  for(int m=0;m<M;m++)
  {
    int id = n*M+m; // m is the row index, that is determines y, id +1 is row below current one
//    int id_test  = m*N+n;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f\n", n, m, I[id]);
#endif


  if ( (m>0) && (m<M-1) )
    Idy[id]  = (I[id+1]-I[id-1])/Scalar(2.0);
  if ( (n>0) && (n<N-1) )
    Idx[id]  = (I[id+M]-I[id-M])/Scalar(2.0);

  if ( (m>0) && (m<M-1) && (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M+1]-I[id+M-1] - I[id-M+1] + I[id-M-1]) /Scalar(4.0);

  }
}

#ifdef DebugOut
//  mexPrintf("First Row\n");
#endif

// along first row ______
//                 ...
for(int n=0;n<N;n++)
{
  int m = 0;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id+1], I[id]);
#endif

  Idy[id]  = (I[id+1]-I[id]);//1.0;

  if ( (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M+1]-I[id+M] - I[id-M+1] + I[id-M]) /Scalar(2.0);
  else
  {
    if ( (n>0) ) // corner upper right
      Idxy[id] = (I[id+1]-I[id] - I[id-M+1] + I[id-M]) /Scalar(2.0);
    else // n== 0 // corner upper left
      Idxy[id] = (I[id+M+1]-I[id+M] - I[id+1] + I[id]) /Scalar(2.0);
  }
}

#ifdef DebugOut
//  mexPrintf("Last Row\n");
#endif
//              ...
// along last row ______
for(int n=0;n<N;n++)
{
  int m = M-1;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, M-1, I[id], I[id-1]);
#endif

  Idy[id]  = (I[id]-I[id-1]);///1.0;

  if ( (n>0) && (n<N-1) )
    Idxy[id] = (I[id+M]-I[id+M-1] - I[id-M] + I[id-M-1]) /Scalar(2.0);
  else
  {
    if ( (n>0) )// corner lower right
      Idxy[id] = (I[id]-I[id-1] - I[id-M] + I[id-M-1]) /Scalar(2.0);
    else// corner lower left
      Idxy[id] = (I[id+M]-I[id+M-1] - I[id] + I[id-1]) /Scalar(2.0);
  }
}

#ifdef DebugOut
//  mexPrintf("First Column\n");
#endif
// along first column |...
for(int m=0;m<M;m++)
{
  int n = 0;
  int id = n*M+m;

#ifdef DebugOut
//  mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id+M], I[id]);
#endif

  Idx[id]  = (I[id+M]-I[id]);//1.0;

  if ( (m>0) && (m<M-1) )
    Idxy[id] = (I[id+M+1]-I[id+M-1] - I[id+1] + I[id-1]) /Scalar(2.0);
  else
  {
    if ( (m>0) )// corner lower left
      Idxy[id] = (I[id+M]-I[id+M-1] - I[id] + I[id-1]) /Scalar(2.0);
    else// corner upper left
      Idxy[id] = (I[id+M+1]-I[id+M] - I[id+1] + I[id]) /Scalar(2.0);
  }
}

#ifdef DebugOut
 // mexPrintf("Last Column\n");
#endif
// along last column ...|
for(int m=0;m<M;m++)
{
  int n = N-1;
  int id = n*M+m;

#ifdef DebugOut
 // mexPrintf("n=%d m=%d = %f - %f\n", n, m, I[id], I[id-M]);
#endif

  Idx[id]  = (I[id]-I[id-M]);//1.0;

  if ( (m>0) && (m<M-1) )
    Idxy[id] = (I[id+1]-I[id-1] - I[id-M+1] + I[id-M-1]) /Scalar(2.0);
  else
  {
    if ( (m>0) ) // corner lower right
      Idxy[id] = (I[id]-I[id-1] - I[id-M] + I[id-M-1]) /Scalar(2.0);
    else// corner upper right
      Idxy[id] = (I[id+1]-I[id] - I[id-M+1] + I[id-M]) /Scalar(2.0);
  }
}
/////////////////////////////////////////////////////////
// END OF APPROXIMATING THE IMMAGE GRADIENTS AS PRE-STEP
/////////////////////////////////////////////////////////
}

private:

  void createMatrices();
  void deleteMatrices();

  void generateMem()
  {
    int mThreads(1);

#ifndef _NO_OPENMP
    mThreads = omp_get_max_threads();
#endif

    cc.resize(mThreads, NULL);
    createMatrices();
    //  std::vector<Scalar**> cc(mThreads, NULL);
    //  for (int i=0;i<mThreads;i++)
    //    cc[i] = matrix(1,4,1,4);
  }

  void freeMem()
  {
    if(Idx !=NULL)
    {
      free (Idx); Idx = NULL;
    }
    if(Idy !=NULL)
    {
      free (Idy); Idy = NULL;
    }
    if(Idxy !=NULL)
    {
      free (Idxy); Idxy = NULL;
    }
  };

  int N;
  int M;

  std::vector<Scalar**> cc;

  Scalar *Idx, *Idy, *Idxy;
  const Scalar* I;
};


template <typename S_>
void BiCubicIp<S_>::createMatrices()
{
  int mThreads(1);

#ifndef _NO_OPENMP
  mThreads = omp_get_max_threads();
#endif

  // dull but else mem leak !!!!
//  std::vector<Scalar**> cc(mThreads, NULL);
  for (int i=0;i<mThreads;i++)
    cc[i] = matrix(1,4,1,4);
}

template <>
void BiCubicIp<float>::createMatrices()
{
  int mThreads(1);

#ifndef _NO_OPENMP
  mThreads = omp_get_max_threads();
#endif

  // dull but else mem leak !!!!
//  std::vector<Scalar**> cc(mThreads, NULL);
  for (int i=0;i<mThreads;i++)
    cc[i] = fmatrix(1,4,1,4);
}

template <typename S_>
void BiCubicIp<S_>::deleteMatrices()
{
  int mThreads(1);

#ifndef _NO_OPENMP
  mThreads = omp_get_max_threads();
#endif

  for (int i=0;i<mThreads;i++)
    free_matrix(cc[i],1,4,1,4);
};

template <>
void BiCubicIp<float>::deleteMatrices()
{
  int mThreads(1);

#ifndef _NO_OPENMP
  mThreads = omp_get_max_threads();
#endif

  for (int i=0;i<mThreads;i++)
    free_fmatrix(cc[i],1,4,1,4);
};


#endif