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
  typedef Math::VectorT<int, 4>    P4i;
  typedef Math::VectorT<Scalar, 3> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::Mat3x3T<Scalar>    M3; 

  MotionCheck( Scalar* p2d_, int nSegments_, int nProposals_ )
    : iK(1.,0.,0.,0.,1.,0.,0.,0.,1.), epsilon(1.0), K(1.,0.,0.,0.,1.,0.,0.,0.,1.),
      penalty(10.0), p2d(p2d_), nSegments(nSegments_), nProposals(nProposals_)
  {epsilon2 = epsilon*epsilon;};

  ~MotionCheck(){};

  void setEpsilon (Scalar epsilon_)           {epsilon = epsilon_;epsilon2 = epsilon*epsilon;}

  void setPenalty (Scalar penalty_)           {penalty = penalty_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  void setK (Scalar* K_)
  {
    // both ways are identical !!!!
    K =M3(K_); 
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
#ifdef _use_patchCenters_
        const P3 c2d( &centers_[ 3*i ] );
        P3 c3d = c2d * (1. / (norm | c2d));
        P3 add = trak + c3d -rotk * c3d;
#else
      P3 add(0.,0.,0.);
#endif

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
//    if ( p2d == NULL )  return;

    scores.clear();
    scores.resize( 2*nPixel ,0);

    const P3& normT = (*normals)[trial];
    const M3& rotkT  = (*rot)[trial];
    const P3& trakT  = (*tra)[trial];
#ifdef _use_patchCenters_
    const P3  c2dT( &centers[ 3*trial ] );
    P3 c3dT = c2dT * (1. / (normT | c2dT));
    P3 addT = trakT + c3dT -rotkT * c3dT;
#else
      P3 addT(0.,0.,0.);
#endif

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
#ifdef _use_patchCenters_
        const P3 c2d( &centers[ 3*segId ] );
        P3 c3d = c2d * (1. / (norm | c2d));
        P3 add = trak + c3d -rotk * c3d;
#else
        P3 add (0.,0.,0.);
#endif
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

  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, int trial)
  {
    const std::vector<int> ids_((int) (mxGetM(Segments) * mxGetN(Segments)), trial);
    computeGeometryConsistency( Segments, centers_,  ids_);
  }

  /// need also centers in 2d (3d)
  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, const std::vector<int> ids_)
  {
    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
	const int skip = 8; // eval every skip'th pixel
	
    // also assume that the max ids is smaller then the amount of normals
    assert(ids_.size() == nSegments);

//    if ( p2d == NULL )      return;

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
#ifdef _use_patchCenters_
      const P3 c2d( &centers_[ 3*i ] );
      P3 c3d = c2d * (1. / (norm | c2d));
      P3 add = trak + c3d -rotk * c3d;
#else
      P3 add(0.,0.,0.);
#endif
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

  void computeGeometryConsistencyLocal(const mxArray*& Segments, const std::vector<int> &ids_, const std::vector<int> &seg2var, P3 mc, Scalar dm, std::vector<Scalar>& _scores)
  {
    int _nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
    Scalar dm2 = dm*dm;
    P3  addDispDiff = K * mc;
    // also assume that the max ids is smaller then the amount of normals
    assert(ids_.size() == nSegments);

//    if ( p2d == NULL )
//      return;

    _scores.clear();
    _scores.resize( nSegments, 0);

//#pragma omp parallel for schedule (static)
    for (int i = 0; i < _nSegments ;i++)
    {
      int localId = seg2var[i];
      if( localId <0 ) continue;
      int k = ids_[i];

      // get rotation, translation, normal
      const P3& norm  = (*normals)[k];// new
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      // get center: // here no cneter stuff:
//      const P3 c2d( &centers_[ 3*i ] );
//      P3 c3d = c2d * (1. / (norm | c2d));
//      P3 add = trak + c3d -rotk * c3d;
      P3 add = trak;

      int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
      int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

      const int step = 8;

      for (int j=0;j<nIds;j+=step)
      {
//        P3 p_3d( ids[j]/height+1, ids[j]%height+1, 1. );// a lot faster !!??
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

        // disp max compi:
        P3 p_img  = K * p_3d_t1;
        P3 p_img2 = p_img + addDispDiff;
        p_img /= p_img[2];p_img2 /= p_img2[2];
        Scalar dist = (p_img - p_img2).sqrnorm();
        //Scalar dist = (p_img - P3( ids[j]/height+1, ids[j]%height+1, 1. )).sqrnorm();

        if (distance2 > epsilon2 || dp < Scalar(1) || dist > dm2)
          _scores[i] += step*penalty; // per segment
//          scores[localId] += 4*penalty;
      }
    }
  }

  void computeGeometryConsistencyLocal(int* segImg, int _nSegments, const std::vector<int> &ids_, const std::vector<int> &seg2var, P3 mc, Scalar dm, std::vector<Scalar>& _scores)
  {
    Scalar dm2 = dm*dm;
    P3  addDispDiff = K * mc;
    // also assume that the max ids is smaller then the amount of normals
    assert(ids_.size() == nSegments);

    _scores.clear();
    _scores.resize( nSegments, 0);
//    std::vector<Scalar> hits(nSegments, 0);
    const int step  = 4;
    const int step2 = step*step;

    for(int i=step/2; i< width;i+=step)
      for(int j=step/2; j< height;j+=step)
      {
      int pix = j+i*height;
      int segId   = segImg[pix];
      int localId = seg2var[segId];
      if( localId <0 ) continue;
      int k = ids_[segId];

      // get rotation, translation, normal
      const P3& norm  = (*normals)[k];// new
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      P3 add = trak;
      P3 p_3d = iK * P3( i+1, j+1, 1. );// a lot faster !!??

      Scalar dp = 1. / (norm | p_3d);
      p_3d      *= dp;

      P3 p_3d_t1 = rotk * p_3d + add;
      P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
        Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), p_3d_t1[2]*p_3d_t1[2]));
#else
        Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), std::min( p_3d_t1[2]*p_3d_t1[2], (p_3d[2]-motion[2] )*(p_3d[2]-motion[2] ) )));
#endif

        // disp max compi:
        P3 p_img  = K * p_3d_t1;
        P3 p_img2 = p_img + addDispDiff;
        p_img /= p_img[2];p_img2 /= p_img2[2];
        Scalar dist = (p_img - p_img2).sqrnorm();
        //Scalar dist = (p_img - P3( ids[j]/height+1, ids[j]%height+1, 1. )).sqrnorm();

        if (distance2 > epsilon2 || dp < Scalar(1) || dist > dm2)
        {
          _scores[segId] += step2*penalty; // per segment
//          hits[segId] += 1.0;
        }
      }
//      for (int i=0;i<_scores.size();i++)
//      {
//      _scores[segId] *= 
//      }
  }

  // ignoring the rotation areounf a center - here only around 0
  void computeGeometryConsistencyPerPixelLocal( int width, int height, int trial, P4i aoi, 
                                                const std::vector<int>& loc2glob, int* segImg, const std::vector<P3>* normalz, 
                                                std::vector<P2>& scoresUn, P3 mc, Scalar dm )
  {
    const P3  normT  = (*normalz)[trial];
    const M3& rotkT  = (*rot)[trial];
    const P3& trakT  = (*tra)[trial];
    const Scalar dm2 = dm*dm;
    const P3 addDispDiff = K * mc;

    int aoix = aoi[0];
    int aoiX = aoi[2];
    int aoiy = aoi[1];
    int aoiY = aoi[3];

//#pragma omp parallel for schedule (static)
    for (int ii=aoix;ii<aoiX;ii++)
    {
      int pos = ii*height;
      for (int j=aoiy;j<aoiY;j++)
      {
        int py = j;
//        int i = pos + j;

        int id_pix1_local = loc2glob[j-aoiy + (aoiY-aoiy)*(ii-aoix)];
        if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

        int px = ii;
        int segId = segImg[ pos + j ];
//        int id_pix_global = pos + j;

        P3 p_2d = iK * P3( px+1, py+1, 1. );

        // get rotation, translation, normal
        const P3& norm  = (*normalz)[segId];
        const M3& rotk  =  (*rot)[segId];
        const P3& trak  =  (*tra)[segId];
      
        Scalar dp = 1. / (norm | p_2d);
        P3 p_3d   = p_2d*dp;

        P3 p_3d_t1 = rotk * p_3d + trak;
        P3 motion  = p_3d_t1 - p_3d;

        // disp max compi 2nd timestep
        P3 p_img  = K * p_3d_t1;
        P3 p_img2 = p_img + addDispDiff;
        p_img /= p_img[2];p_img2 /= p_img2[2];
        Scalar dist = (p_img - p_img2).sqrnorm();        

#ifndef __doubleCheck__
        Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
        Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

        P2 pen(0,0);

        if (distance > epsilon || dp < Scalar(1) || dist > dm2)
          pen[0] = penalty;
        
        ///////////////////////////////////////////////////////////////////

        dp      = 1. / (normT | p_2d);
        p_3d    = p_2d*dp;

        p_3d_t1 = rotkT * p_3d + trakT;
        motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
        distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
//        distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
        distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

        // disp max compi:
        p_img  = K * p_3d_t1;
        p_img2 = p_img + addDispDiff;
        p_img /= p_img[2];p_img2 /= p_img2[2];
        dist = (p_img - p_img2).sqrnorm();

        if (distance > epsilon || dp < Scalar(1) || dist > dm2)
          pen[1] = penalty;
          
        scoresUn[ id_pix1_local ] = pen;
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
  M3  K;

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
