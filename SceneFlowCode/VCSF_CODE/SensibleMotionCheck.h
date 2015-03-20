/////////////////////////////////////////////////////////
/// put high score on non-sense motions == super fast ///
/////////////////////////////////////////////////////////
#ifndef __ENERGY_MOTIONCHECKER__
#define __ENERGY_MOTIONCHECKER__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include "mex.h" // include after VectorT and Mat3x3T
#include <vector>

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

  /// used for view consistent case
  void computeGeometryConsistencyLocal(int* segImg, int _nSegments, const std::vector<int> &ids_, 
                                       const std::vector<int> &seg2var, P3 mc, 
                                       Scalar dm, std::vector<Scalar>& _scores);

  // ignoring the rotation around a center - here only around 0
  /// used for view consistent case per pixel
  void computeGeometryConsistencyPerPixelLocal( int width, int height, int trial, P4i aoi, 
                                                const std::vector<int>& loc2glob, int* segImg, const std::vector<P3>* normalz, 
                                                std::vector<P2>& scoresUn, P3 mc, Scalar dm );

  /// used for standard algorithm
  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, int trial)
  {
    const std::vector<int> ids_((int) (mxGetM(Segments) * mxGetN(Segments)), trial);
    computeGeometryConsistency( Segments, centers_,  ids_);
  }

  /// need also centers in 2d (3d)
  void computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, const std::vector<int> ids_);

  /// used for standard algorithm, per pixel
  void computeGeometryConsistencyPerPixel(int width, int height, int trial, int aoix, int aoiX, int aoiy, int aoiY, 
                                          std::vector<int>& loc2glob, int nPixel, int* segImg, Scalar* centers );

  std::vector<Scalar>& getScores () { return scores; }

  void setWidthHeight(int width_, int height_) {width = width_; height = height_;}

private:

  //////////////////////////////////////////////////

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

#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_MOTIONCHECKER__CPP)
#include "SensibleMotionCheck.cpp"
#endif

#endif
