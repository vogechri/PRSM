///////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#ifndef __ENERGY_MOTIONCHECKER__
#define __ENERGY_MOTIONCHECKER__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

//#define _Debug_out_

#include <map>
#include <vector>

using namespace std;
using namespace Math;

//#define __doubleCheck__

//////////////////////////////////////////////////////////////////////////////////////

/// provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
template<typename Scalar> class MotionCheck
{
public:

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 3> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::Mat3x3T<Scalar>    M3; 

  MotionCheck( Scalar* p2d_, int nSegments_, int nProposals_ )
    : iK(1.,0.,0.,0.,1.,0.,0.,0.,1.), epsilon(1.0), penalty(50.0), p2d(p2d_), nSegments(nSegments_), nProposals(nProposals_)
  {epsilon2 = epsilon*epsilon;};

  ~MotionCheck(){};

  void setEpsilon (Scalar epsilon_)           {epsilon = epsilon_;epsilon2 = epsilon*epsilon;}

  void setPenalty (Scalar penalty_)           {penalty = penalty_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  void setK (Scalar* K_)
  {
    // both ways are identical !!!!
    iK=M3(K_);
    iK.invert();
  }

//  const std::vector<P3>& getTra()      { return tra; }
//  const std::vector<M3>& getRot()      { return rot; }
//  const std::vector<P3>& getCenters()  { return centers; }

  /// need also centers in 2d (3d)
  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_ )
  {
    if ( p2d == NULL )
      return;

    scores.clear();
    scores.resize( nSegments*nProposals ,0);

    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

#pragma omp parallel for schedule (static)
    for (int k = 0;k< nProposals;k++)
    {
      // get rotation, translation, normal
      const P3& norm = (*normals)[k];
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      for (int i = 0; i < nSegments ;i++)
      {
        // get center:
        const P3 c2d( &centers_[ 3*i ] );
        P3 c3d = c2d * (1. / (norm | c2d));
        P3 add = trak + c3d -rotk * c3d;

        int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
        int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));
        
        for (int j=0;j<nIds;j++)
        {
          P3 p_3d   = P3(&p2d[3*ids[j]]);
          Scalar dp = 1. / (norm | p_3d);
          p_3d      *= dp;
          
//          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
          P3 p_3d_t1 = rotk * p_3d + add;
          P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
          Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
          Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

          if (distance > epsilon || dp < Scalar(1))
            scores[ k*nSegments + i] += penalty;

        }
      }
    }
  }


  void computeGeometryConsistencyPerPixel(int width, int height, int trial, int aoix, int aoiX, int aoiy, int aoiY, 
                                          std::vector<int>& loc2glob, int nPixel, int* segImg, Scalar* centers )
  {
    if ( p2d == NULL )
      return;

    scores.clear();
    scores.resize( 2*nPixel ,0);

    const P3& normT = (*normals)[trial];
    const M3& rotkT  = (*rot)[trial];
    const P3& trakT  = (*tra)[trial];
    const P3  c2dT( &centers[ 3*trial ] );
    P3 c3dT = c2dT * (1. / (normT | c2dT));
    P3 addT = trakT + c3dT -rotkT * c3dT;

#pragma omp parallel for schedule (static)
    for (int j=aoiy;j<aoiY;j++)
    {
      int py = j;
      int pos = j*width;
      int i = pos+aoix;
      for (int ii=aoix;ii<aoiX;ii++,i++)
      {
        int id_pix1_local = loc2glob[ii-aoix + (aoiX-aoix)*(j-aoiy)];
        if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

        int px = ii;
        int segId = segImg[ i ];
        int id_pix_global = i;
//        P3 centerA ( &centers[3*segId] );

#ifdef _DEBUG
      if(py==278 && px==1059 )
        int breakhere = 1;
#endif

        P3 p_2d = iK * P3( px, py, 1. );

        // get rotation, translation, normal
        const P3& norm = (*normals)[segId];
        const M3& rotk  = (*rot)[segId];
        const P3& trak  = (*tra)[segId];

        // get center:
        const P3 c2d( &centers[ 3*segId ] );
        P3 c3d = c2d * (1. / (norm | c2d));
        P3 add = trak + c3d -rotk * c3d;
      
        //      P3 p_3d   = P3(&p2d[3*ids[j]]);
        Scalar dp = 1. / (norm | p_2d);
        P3 p_3d   = p_2d*dp;
          
        //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
        P3 p_3d_t1 = rotk * p_3d + add;
        P3 motion  = p_3d_t1 - p_3d;
        
#ifndef __doubleCheck__
        Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
        Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

        if (distance > epsilon || dp < Scalar(1))
          scores[ 2*id_pix1_local ] = penalty;
        
        ///////////////////////////////////////////////////////////////////

        dp = 1. / (normT | p_2d);
        p_3d    = p_2d*dp;
          
        //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
        p_3d_t1 = rotkT * p_3d + addT;
        motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
        distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
        distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

        if (distance > epsilon || dp < Scalar(1) )
          scores[ 2*id_pix1_local+1 ] = penalty;

      }
    }
  }

  /// need also centers in 2d (3d)
  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, const std::vector<int> ids_)
  {
    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
	const int skip = 8; // eval every skip'th pixel
    // also assume that the max ids is smaller then the amount of normals
    assert(ids_.size() == nSegments);

    if ( p2d == NULL )
      return;

    scores.clear();
    scores.resize( nSegments, 0);

#pragma omp parallel for schedule (static)
    for (int i = 0; i < nSegments ;i++)
    {
      int k = ids_[i];

#ifdef _DEBUG
      if(k==1689 )
        int breakhere = 1;
#endif

      // get rotation, translation, normal
      const P3& norm = (*normals)[k];
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      // get center:
      const P3 c2d( &centers_[ 3*i ] );
      P3 c3d = c2d * (1. / (norm | c2d));
      P3 add = trak + c3d -rotk * c3d;

      int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
      int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

      for (int j=0;j<nIds;j+=skip)
      {
        P3 p_3d = iK * P3( ids[j]/height+1, ids[j]%height+1, 1. );// a lot faster !!??
//        P3 p_3d = P3(&p2d[3*ids[j]]);

        Scalar dp = 1. / (norm | p_3d);
        p_3d      *= dp;

        //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
        P3 p_3d_t1 = rotk * p_3d + add;
        P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
//        Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), p_3d_t1[2]));
        Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), p_3d_t1[2]*p_3d_t1[2]));
#else
        Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), std::min( p_3d_t1[2]*p_3d_t1[2], (p_3d[2]-motion[2] )*(p_3d[2]-motion[2] ) )));
#endif

        if (distance2 > epsilon2 || dp < Scalar(1))
          scores[i] += Scalar(skip)*penalty;

      }
    }
  }


  std::vector<Scalar>& getScores () { return scores; }

  void setWidthHeight(int width_, int height_) {width = width_; height = height_;}

private:

  //////////////////////////////////////////////////

  /// the centers of the pathces in 2d
  //  std::vector<P3> centers;
//  std::vector<P3> centersA;
//  std::vector<P3> centersB;

  M3 iK;

  const std::vector<P3>* normals;

  /// translations
  const std::vector<P3>* tra;
  /// rotations
  const std::vector<M3>* rot;

  ///  a penalty if geometry violates constraints 0 else
  std::vector<Scalar>    scores;

  /// motion 3d distance ratio to be smaller than epsilon
  Scalar epsilon;

  /// epsilon^2
  Scalar epsilon2;

  /// penalty on segment if it is larger 
  Scalar penalty;

  /// 2d points corresponding to the pixel
  Scalar* p2d;

  int    nSegments;

  int    nProposals;

  int width;
  int height;
};
#endif
