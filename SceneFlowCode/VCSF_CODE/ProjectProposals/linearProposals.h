////////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#ifndef __LINEAR_PROPOSALS__
#define __LINEAR_PROPOSALS__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

// for growing the maps
#include "Templates/HeapT.h"
#include "Segment_HeapHelper.h"

#ifndef _NO_OPENMP
#include <omp.h>
#endif

//#include "mex.h" // include after VectorT and Mat3x3T

#include <map>
#include <vector>

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

/// provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
template<typename Scalar> class linearMotionProposals
{
public:

  typedef HeapStoreEntryT<double>  HeapStoreEntry;
  typedef vector< HeapStoreEntry > HeapEntryStore;

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 3> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::Mat3x3T<Scalar>    M3; 

  // Scalar* p2d_,
  linearMotionProposals(  int nSegments_, int nProposals_, int w_, int h_)
    : iK(1.,0.,0.,0.,1.,0.,0.,0.,1.), nSegments(nSegments_), nProposals(nProposals_),
    width(w_), height(h_), 
    _heapEntryStore(), _heapInterface(&_heapEntryStore), _heap(_heapInterface),//, p2d(p2d_)
    centeredView(true)
  {};

  ~linearMotionProposals(){};

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  void useCentred( bool onOff) {centeredView = onOff;};


  void setK (Scalar* K_)
  {
    K = M3(K_);
    iK=M3(K_);
    iK.invert();
  }

  // projected and given in the form of new segments:
  std::vector<P3>& getNewN()        { return propNormals; };
  std::vector<P3>& getNewTra()      { return propTranslations; };
  std::vector<M3>& getNewRot()      { return propRotations; };
  std::vector<int> getPMap()        { return proposalMap;}; 

  // projected only, reasonable if center (projected spot is given)
  std::vector<P2>& getNewCenters ()  { return newCenters; };
  std::vector<P3>& getProjN()        { return newNormals; };
  std::vector<P3>& getProjTra()      { return newTranslations; };
  std::vector<M3>& getProjRot()      { return newRotations; };

  /// need also centers in 2d (3d), ids_ : current solution: segment i received proposal ids_[i]
  void computeNew3dmovingplanes(int nSegments, Scalar* centers_, const std::vector<int> ids_ )
  {
    newCenters.clear();
    newCenters.resize( (int) nSegments, P2(0,0));
    newNormals.clear();
    newNormals.resize( (int) nSegments, P3(0,0,1));
    newTranslations.clear();
    newTranslations.resize( (int) nSegments, P3(0,0,0));
    newRotations.clear();
    newRotations.resize( (int) nSegments, M3(1,0,0, 0,1,0, 0,0,1));

//    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

#pragma omp parallel for schedule (static)
    for (int i = 0; i < nSegments ;i++)
    {
      int k = ids_[i];

      // get rotation, translation, normal
      const P3& norm = (*normals)[k];
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      // get center:
      const P3 c2d( &centers_[ 3*i ] );
      P3 c3d = c2d * (1. / (norm | c2d));
      P3 add = trak + c3d - rotk * c3d;

      if ( !centeredView )
        add = trak;

      P3 p_3d_1    = c2d + P3(1.,0,0);
      Scalar dp_1  = 1. / (norm | p_3d_1);
      p_3d_1      *= dp_1;
      P3 p_3d_1_t1 = rotk * p_3d_1 + add;

      P3 p_3d_2    = c2d + P3(0,1.,0);
      Scalar dp_2  = 1. / (norm | p_3d_2);
      p_3d_2      *= dp_2;
      P3 p_3d_2_t1 = rotk * p_3d_2 + add;

      ///////////////////// 
      P3  newC3d = (c3d+trak);
      if ( !centeredView )
        newC3d = rotk*c3d + trak;

      P3  projC  = K * newC3d;
      newCenters[i] =  P2(projC[0]/projC[2], projC[1]/projC[2]);
      newTranslations[i] = trak;
      newRotations[i] = rotk;

      M3 mat ( newC3d[0], newC3d[1], newC3d[2], p_3d_1_t1[0]-newC3d[0], p_3d_1_t1[1]-newC3d[1], p_3d_1_t1[2]-newC3d[2], p_3d_2_t1[0]-newC3d[0], p_3d_2_t1[1]-newC3d[1], p_3d_2_t1[2]-newC3d[2] );
      mat.invert();
//      P3 newN = mat* P3(1,0,0);
      P3 newN = mat(0);
      Scalar test1 = newN|c3d;
      Scalar test2 = newN|p_3d_1;
      Scalar test3 = newN|p_3d_2;

      newNormals[i] = newN;
    }
  }

  /// need also centers in 2d (3d), ids_ : current solution: segment i received proposal ids_[i]
  void noMotionMovingPlanes(int nSegments, Scalar* centers_, const std::vector<int> ids_ )
  {
    newCenters.clear();
    newCenters.resize( nSegments, P2(0,0));
    newNormals.clear();
    newNormals.resize( nSegments, P3(0,0,1));
    newTranslations.clear();
    newTranslations.resize( nSegments, P3(0,0,0));
    newRotations.clear();
    newRotations.resize( nSegments, M3(1,0,0, 0,1,0, 0,0,1));

//    nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));

#pragma omp parallel for schedule (static)
    for (int i = 0; i < nSegments ;i++)
    {
      int k = ids_[i];

      // get rotation, translation, normal
      const P3& norm = (*normals)[k];
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      // get center:
      const P3 c2d( &centers_[ 3*i ] );

      P3  projC  = K * c2d;
      newCenters[i] =  P2(projC[0]/projC[2], projC[1]/projC[2]);
      newTranslations[i] = trak;
      newRotations[i] = rotk;

      newNormals[i] = norm;
    }
  }

  void adjustNormals()
  {
    for(int i=0;i<newNormals.size();i++)
      newNormals[i] = newNormals[i]*Scalar(-1);
  }

  // edges yield the neighbor segments 
  /*! idea : first newCenters are projected into closest segment: store closest distance and id per segment
   * there will be empty segments to begin with - not important
   * use Q: pick closest, update neighbors (or put into Q)
   * finally eeach segment has an id == proposal
  */
  void generateProposals(const mxArray* edges_, int* segImg, int nSegments, Scalar* centers_t1)
  {
    std::vector<int> bestId(nSegments, -1); // -1 : none yet
    std::vector<int> fixed( nSegments,  0); // 0 free, 1: in Q, else fixed
    std::vector<Scalar> distances( nSegments, width*height*width );
    for( int i =0; i < newCenters.size(); i++ )
    {
      P2 propCenter = newCenters[i];
      P2 pc2d = P2( max(0.,min(width-1., propCenter[0])), max(0.,min(height-1., propCenter[1]) ) );

      int pixPos = (int)(pc2d[0]+0.5) * height + (int)(pc2d[1]+0.5);
      int initSeg = segImg[pixPos];

      P3 segCenter = K * P3(&centers_t1[3*initSeg]);
      P2 seg2DCenter = P2( segCenter[0]/segCenter[2], segCenter[1]/segCenter[2] );

      Scalar distance = ( propCenter - seg2DCenter ).sqrnorm();
      if ( initSeg<bestId.size() && ( bestId[initSeg] == -1 || distances[initSeg] > distance ))
      {
        bestId[initSeg] = i;
        distances[initSeg] = distance;
      }
    }

    // now put into Q:
    _heap.clear();
    _heapEntryStore.clear();
//      _heapInterface.clear();
    // cost w*h (larger than maximal) , startSeg : -1 == not touched yet, -2: fixed already 
    _heapEntryStore.resize( nSegments, HeapStoreEntry(width*height*width, -1) );

    // prep intialization
    for (int k=0; k< bestId.size(); k++)
    {
      if (bestId[k]<0) continue;
      HeapStoreEntry& hse = _heapEntryStore[k];
      hse.setStartSeg(bestId[k]); 
      hse.setCost( distances[k] );
      _heap.insert( HeapEntryD (k) );
      fixed[k]=1; // in Q
    }

      int numRed(0);
      while ( !_heap.empty() )
      {
        numRed++;

        HeapEntryD id = _heap.front();
        _heap.pop_front();
        int cSeg = id.pId();
        HeapStoreEntry& heapE = _heapEntryStore[ cSeg ];

        fixed[cSeg] = -1; // fixed
        // fix segment , add neighbors if smaller

        int sSeg = heapE.startSeg();
        bestId[cSeg] = sSeg;
        P2 propCenter = newCenters[sSeg];

        // loop neighs:
        Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, cSeg) );
        int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, cSeg) )/5;
        for (int j =0; j < nIds ;j++)
        {
          int nId = (int) (edge [5*j ])-1;
          if ( fixed[nId] == -1) continue; // already fix

          HeapStoreEntry& hse = _heapEntryStore[nId];
          int bestSeg = hse.startSeg();

          P3 segCenter = K * P3(&centers_t1[3*nId]);
          P2 seg2DCenter = P2( segCenter[0]/segCenter[2], segCenter[1]/segCenter[2] );  
          Scalar distance = ( propCenter - seg2DCenter ).sqrnorm();

          // new
          if ( fixed[nId] == 0)
          {
            hse.setCost( distance );
            hse.setStartSeg( sSeg );
            _heap.insert( HeapEntryD (nId) );
            fixed[nId] = 1;
          }
          else // update if lower only
          {
            if ( distance< hse.cost() )
            {
              hse.setCost( distance );
              hse.setStartSeg( sSeg );
              _heap.update( HeapEntryD( nId ) );
            }
          }
        }
      } // while heap not empty

    proposalMap = bestId;

    propNormals.clear();
    propNormals.resize( nSegments, P3(0,0,1));
    propTranslations.clear();
    propTranslations.resize( nSegments, P3(0,0,0));
    propRotations.clear();
    propRotations.resize( nSegments, M3(1,0,0, 0,1,0, 0,0,1));

    for(int i=0;i<nSegments;i++)
    {
      propNormals[i]      = newNormals[proposalMap[i]]*Scalar(-1);
      propTranslations[i] = newTranslations[proposalMap[i]];
      propRotations[i]    = newRotations[proposalMap[i]];
    }

  }

private:

  //////////////////////////////////////////////////

  /// the centers of the pathces in 2d
  //  std::vector<P3> centers;

  M3 K;
  M3 iK;

  const std::vector<P3>* normals;
  /// translations
  const std::vector<P3>* tra;
  /// rotations
  const std::vector<M3>* rot;

  std::vector<P2>    newCenters;

  std::vector<P3>    newNormals;

  std::vector<M3>    newRotations;

  std::vector<P3>    newTranslations;

  std::vector<P3>    propNormals;

  std::vector<M3>    propRotations;

  std::vector<P3>    propTranslations;

  /// 2d points corresponding to the pixel
//  Scalar* p2d;

  int    nSegments;

  int    nProposals;

  int width, height;

  std::vector<int> proposalMap; 

  /// whether or not use rotation around center or origin
  bool centeredView;

  ////////////////////
    /// stores relevant information about HE collapse
  HeapEntryStore _heapEntryStore;

  /// the heap-Interface structure
  HeapInterfaceD < HeapEntryD, HeapEntryStore > _heapInterface;

  /// the heap for selecting the element with the least cost
  Utils::HeapT< HeapEntryD, HeapInterfaceD < HeapEntryD, HeapEntryStore > > _heap;
};
#endif //__LINEAR_PROPOSALS__
