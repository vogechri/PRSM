/*
Copyright (c) 2013, Christoph Vogel, ETH Zurich

The code may be used free of charge for non-commercial and
educational purposes, the only requirement is that this text is
preserved within the derivative work. For any other purpose you
must contact the authors for permission. This code may not be
redistributed without written permission from the authors.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES 
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE 
FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY 
DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, 
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, 
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <assert.h>
#include <vector>
#include <set>
#include <algorithm>
#include <omp.h>
#include "mex.h"
#include "lapack.h"
#include "MatMatl2.h"
#include "GridElement.h"

/// matlab calling
void growGrid ( int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{

  typedef double Scalar;
  Scalar *output1;

  Scalar *Matrices = (Scalar*)   mxGetPr(prhs[0]);//
  Scalar *Vectors  = (Scalar*)   mxGetPr(prhs[1]);//
  Scalar *Scalars  = (Scalar*)   mxGetPr(prhs[2]);//
  Scalar  nSegs    = (Scalar ) (*mxGetPr(prhs[3]));//

  // the segment neighbors
  const mxArray* Segments  = (prhs[4]);
  int elements  = (int)mxGetNumberOfElements(prhs[2]);
  int nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

  Scalar  nGoal    = (Scalar ) (*mxGetPr(prhs[5]));

  Scalar *pps      = NULL;
  Scalar  nErr     = 100;
  if ( nrhs > 6)
  {
    pps      = (Scalar*)   mxGetPr(prhs[6]);// pixel per segment 
    nErr     = (Scalar ) (*mxGetPr(prhs[7]));
  }

  ///////////////////////
  int nsubs=mxGetNumberOfDimensions(prhs[0]);
  const mwSize* dims = mxGetDimensions(prhs[0]);
  int dx = dims[0];
  int dy = dims[1];
  int dz = dims[2];

  printf("Elements: %d, %d\n", elements, dz);

  plhs[0] = mxCreateDoubleMatrix( elements, 1, mxREAL);
  output1  = mxGetPr(plhs[0]);

  if (dx==6)
  {
    std::vector< Matmatvecl2_d6* > mat_container( elements, NULL );
    if ( pps  == NULL )
      for (int i=0; i< elements;i++)
    	  mat_container[i] = new Matmatvecl2_d6( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]) );
    else
      for (int i=0; i< elements;i++)
      	mat_container[i] = new Matmatvecl2_d6( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]), &(pps[i]) );

    grid myGrid( nSegments, Segments, mat_container, nGoal );
    myGrid.copy_seg ( output1 );
    for (int i=0; i< elements;i++)
    	delete( mat_container[i]); 
  }

  if (dx==3)
  {
    std::vector< Matmatvecl2_d3* > mat_container( elements, NULL );
    if ( pps  == NULL )
      for (int i=0; i< elements;i++)
    	  mat_container[i] = new Matmatvecl2_d3( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]) );
    else
      for (int i=0; i< elements;i++)
      	mat_container[i] = new Matmatvecl2_d3( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]), &(pps[i]) );

    grid myGrid( nSegments, Segments, mat_container, nGoal );
    myGrid.copy_seg ( output1 );
    for (int i=0; i< elements;i++)
  	  delete( mat_container[i]); 
  }

  if (dx==9)
  {
    std::vector< Matmatvecl2_d9* > mat_container( elements, NULL );
    if ( pps  == NULL )
      for (int i=0; i< elements;i++)
    	  mat_container[i] = new Matmatvecl2_d9( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]) );
    else
      for (int i=0; i< elements;i++)
      	mat_container[i] = new Matmatvecl2_d9( &Matrices[dx*dy*i], &Vectors[dy*i], &(Scalars[i]), &(pps[i]) );

    grid myGrid( nSegments, Segments, mat_container, nGoal );
    myGrid.copy_seg ( output1 );
    for (int i=0; i< elements;i++)
    	delete( mat_container[i]); 
  }
}