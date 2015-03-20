////////////////////////////////////////////////////////////////
// handle occlusions for an image pair: z-Buffer with memory //
// generate lists of (potential) occluders                   //
////////////////////////////////////////////////////////////////

#ifndef __OcclusionMappingBuffer__h
#define __OcclusionMappingBuffer__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#include <list>
#include <set>
#include <limits>
#include <math.h>
#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include "QPBO_Converter.h" // convert occusion configuration to binary terms

using namespace std;
using namespace Math;

typedef std::pair <int,int>  PixSegPair;

/// the list of the occPixelList's has to be sorted by the ownIds to be able to merge lists
template<typename Scalar>
struct occPixelList
{
  typedef std::pair <int,int>  PixSegPair;
  occPixelList( ) : ownId(-1), segId(-1), occluders(), fullOccluder(-1) //, outerOccluder(-1)
  {};

  occPixelList( int _pid, int _segid, std::vector < PixSegPair > _occluders, int _fullOccluder )
    : ownId(_pid), segId(_segid), occluders(_occluders), fullOccluder(_fullOccluder), occludes(-1)
  {};

  ~occPixelList() {};
  ////////////////////////////////////////////

  /// stored are pairs: pixel id, assign a 0 or 1, interpretation 0: occludes if id=0, encode x_i in poly: if x_i=1 not occluded
  std::vector < PixSegPair > occluders;

  int ownId;        // id of the occluded pixel GLOBAL ? NO!! - LOCAL !

  int occludes;

  int segId;        // seg assignment to pixel

  /// the pixel is always occluded - either from outside or from a fre pixel - no matter what is assigned ot the pixel
  int fullOccluder; // needed to encode a general occlusion (pixel id not really needed)
  // replace data term by occ pen or add unary poly with (occ-data)
};

template<typename Scalar>
struct oMapBufferStorage
{
public:

  oMapBufferStorage() : pixPos(-1), depth(-1.0), segId(-1), px(-1), py(-1) {};

  // with segid interesting for the lookup of new id assignments (a 1 to a pixel)
  oMapBufferStorage( int _pixPos, int _px, int _py, Scalar _depth, int _segId )
    : pixPos(_pixPos), px(_px), py(_py), depth(_depth), segId(_segId)
  {};

  ~oMapBufferStorage() {};
  ///////////////////////////////////////
  int px;
  int py;
  int pixPos;
  Scalar depth;
  int segId;
};

template<typename Scalar>
struct oMapBufferLookup
{
public:

  oMapBufferLookup( ) : vecPos(-1), pixPos(-1) {};// oob
  oMapBufferLookup( int _pixPos, int _vecPos) :  pixPos(_pixPos), vecPos(_vecPos)
  {};

  ~oMapBufferLookup() {};
  ///////////////////////////////////////
  int pixPos;
  int vecPos;
};
//////////////////////////////////////////////////////////////////////////

/*
call like this:
initNewSeg( int segid1, int _startX, int _endX, int _startY, int _endY )

generateOcclusionLists( int segId, int aoix, int aoiX, int aoiy, int aoiY )

for all pId in the image: with pid being the global id
//getOcclusionLists( occPixelList< Scalar >& pixList0, occPixelList< Scalar >& pixList1, int pId, int expSeg, int aoix, int aoiX, int aoiy, int aoiY )

finishNewSeg( int segid1, int startX, int endX, int startY, int endY, std::vector<int> changes )
*/

///////////////////////////////////////
template<typename Scalar> class OcclusionMapBuffer
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
//  typedef Math::VectorT<Scalar, 5>  P5;

  typedef std::pair <int,int>       Seg_ij;


  OcclusionMapBuffer(int _maxPatchSize, bool doOcclusionStuff = true)
    : segImg(NULL), w(0), h(0), occPenalty(0), p2d(NULL), newIteration(0), 
    startX(0), endX(0), startY(0), endY(0), maxPatchSize(_maxPatchSize), locIds(NULL), genOccData(doOcclusionStuff)
  {}

  OcclusionMapBuffer(int _maxPatchSize, int w_, int h_, bool doOcclusionStuff = true)
    : segImg(NULL), w(w_), h(h_), occPenalty(0), p2d(NULL), newIteration(0), 
    startX(0), endX(0), startY(0), endY(0), maxPatchSize(_maxPatchSize), locIds(NULL), genOccData(doOcclusionStuff)
  {}

  ~OcclusionMapBuffer(){};

  void setP2d (Scalar* p2d_) { p2d = p2d_;}

    /// return the list of occlusions for segments assigned a 0
  std::list< occPixelList< Scalar > >& getOccList0()
  {return freeOcclusionsSeg0;}

  /// return the list of occlusions for segments assigned a 1
  std::list< occPixelList< Scalar > >& getOccList1()
  {return freeOcclusionsSeg1;}

  /// return the list of outer occlusions
  std::list< occPixelList< Scalar > >& getOuterList()
  {return outerOcclusions;}

  void setLocalIds (std::vector<int>& _locIds) {locIds = &_locIds;}

  /// the lists of the occPixelList's has to be sorted by the ownIds to be able to merge lists
  void generateOcclusionLists( int segId, int aoix, int aoiX, int aoiy, int aoiY );

  // 
  void init();

  /// convert to directly usable stuff
  void setHomvNom( Scalar* homs_, Scalar* vNoms_ );

  void setSegImg( int _w, int _h, int* _segImg, int _nSegments );

  void SetOccPenalty (Scalar _occPenalty) {occPenalty = _occPenalty;};
//  void setOccSegImg ( int* _occSegImg ) { occSegImg = _occSegImg;};

  ////////////////////////////////////////////////////////////////////////////
  // functions which are called several times
  ////////////////////////////////////////////////////////////////////////////

  /// clears the part out of the 
  void initNewSeg( int segid1, int _startX, int _endX, int _startY, int _endY );

  // update : the other pixel have an id in the mrf stored here HERE updates must be 0: no update 1:changed
  /// clears some stuff and updates according to segment map, assume updates encodes a 0 for no change and 1 else
//  void finishNewSeg( int segid1, int startX, int endX, int startY, int endY, std::vector<int>& changes )
  void finishNewSeg( int segid1, std::vector<int>& loc2glob, std::vector<int>& changes );

  inline int isFreeOld( const oMapBufferStorage<Scalar>& info, int expSeg, int aoix, int aoiX, int aoiy, int aoiY )
  {
//    if (info.segId != expSeg && info.px >= aoix && info.px < aoiX && info.py >= aoiy && info.py < aoiY)
//      return (*locIds)[(aoiX-aoix)*(info.py-aoiy)+(info.px-aoix)]; 
    int testY = info.pixPos/w;
    int testX = info.pixPos%w;

    if (info.segId != expSeg && (testX >= aoix) && (testX < aoiX) && (testY >= aoiy) && (testY < aoiY))
      return (*locIds)[(aoiX-aoix)*(testY-aoiy)+(testX-aoix)]; 
    else
      return -1;
  }

  inline int isFree( const oMapBufferStorage<Scalar>& info, int aoix, int aoiX, int aoiy, int aoiY )
  {
//    if (info.segId != expSeg && info.px >= aoix && info.px < aoiX && info.py >= aoiy && info.py < aoiY)
//      return (*locIds)[(aoiX-aoix)*(info.py-aoiy)+(info.px-aoix)]; 
    int testY = info.pixPos/w;
    int testX = info.pixPos%w;

    if ( (testX >= aoix) && (testX < aoiX) && (testY >= aoiy) && (testY < aoiY) )
      return (*locIds)[(aoiX-aoix)*(testY-aoiy)+(testX-aoix)]; 
    else
      return -1;
  }

  // projPix: the pixel is projected onto that pixel, info the list of pixel assigned a 0, falling on the same pixel 
bool constructList0( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY);

  // projPix: the pixel is projected onto that pixel, info the list of pixel assigned a 0, falling on the same pixel 
  bool constructList1( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY);

  /// assumes that pId is already the id of a pixel closest to the camera of all fixed pixel and occluded by one (or more) free pixel ! 
  bool generateExtraOcclusions ( int pId, occPixelList< Scalar >& pixList, int expSeg, int aoix, int aoiX, int aoiy, int aoiY );

  /// returns a list for each pixel asked, expSeg is to check if a pixel within the update area aoi is fixed
  int getOcclusionLists( occPixelList< Scalar >& pixList0, occPixelList< Scalar >& pixList1, int pId, int expSeg, int aoix, int aoiX, int aoiy, int aoiY, int locId );

  void getCurrentOcclusions( Scalar* oMap );

  int getNPixel() { return (endY-startY) * (endX-startX);}; 

  void getFullIdx(std::vector< Scalar >& idx)
  {
    idx.resize(w*h);
    std::copy( idX.begin() ,idX.end(), idx.begin() );
  }

  void getFullIdy(std::vector< Scalar >& idy)
  {
    idy.resize(w*h);
    std::copy( idY.begin(), idY.end(), idy.begin() );
  }

  void getIdx(std::vector< Scalar >& idx)
  {
    idx.resize((endY-startY) * (endX-startX));
//    std::copy( idX.begin() ,idX.end(), idx.begin() );
    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++,i++)
        idx[i] = idX[i];
    }
  }

  void getIdy(std::vector< Scalar >& idy)
  {
    idy.resize((endY-startY) * (endX-startX));
//    std::copy( idY.begin(), idY.end(), idy.begin() );

    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++,i++)
        idy[i] = idY[i];
    }
  }

  void getSegIdx(std::vector< Scalar >& idx)
  {
    idx.resize((endY-startY) * (endX-startX));
    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++,i++)
        idx[i] = idXCurSeg[i];
    }
  }

  void getSegIdy(std::vector< Scalar >& idy)
  {
    idy.resize((endY-startY) * (endX-startX));
    for (int j=startY;j<endY;j++)
    {
      int pos = j*w;int i = pos+startX;
      for (int ii=startX;ii<endX;ii++,i++)
        idy[i] = idYCurSeg[i];
    }
  }


  /// input is a vector with global indices which we want to look up
  void getIdx(std::vector< Scalar >& idx, const std::vector< int >& globFree)
  {
    idx.resize(globFree.size());
    for (int j=0;j<globFree.size();j++)
        idx[j] = idX[globFree[j]];
  }

  void getIdy(std::vector< Scalar >& idy, const std::vector< int >& globFree)
  {
    idy.resize(globFree.size());
    for (int j=0;j<globFree.size();j++)
      idy[j] = idY[globFree[j]];
  }

  /// input is a vector with global indices which we want to look up
  void getSegIdx(std::vector< Scalar >& idx, const std::vector< int >& globFree)
  {
    idx.resize(globFree.size());
    for (int j=0;j<globFree.size();j++)
      idx[j] = idXCurSeg[ globFree[j] ];
  }

  void getSegIdy(std::vector< Scalar >& idy, const std::vector< int >& globFree)
  {
    idy.resize(globFree.size());
    for (int j=0;j<globFree.size();j++)
      idy[j] = idYCurSeg[ globFree[j] ];
  }

private:

  ///////////////////////////////////////////////

  /// from the pixel ina the reference image to the pixel it is projected onto
  std::vector<oMapBufferLookup<Scalar> > lookup0;
  std::vector<oMapBufferLookup<Scalar> > lookup1;

  /// information about the pixel in the data image: depth for occlusions and id of the original pixel
  std::vector< std::vector< oMapBufferStorage<Scalar> > > store0;
  /// same as above but instead of the id also a iteration number 
  std::vector< std::vector< oMapBufferStorage<Scalar> > > store1;

  /// the x-coordinate of the warped pixels
  std::vector< Scalar > idX;
  /// the y-coordinate of the warped pixels
  std::vector< Scalar > idY;

  /// the x-coordinate of the warped pixels by the current segment
  std::vector< Scalar > idXCurSeg;
  /// the y-coordinate of the warped pixels by the current segment
  std::vector< Scalar > idYCurSeg;


  int startX; 
  int endX; 
  int startY;
  int endY;

  /// helper variable to identify correct occlusions in the current setting
  int newIteration;

  /// width
  int w;
  /// height
  int h;

  /// the amount of segmetns in the images
  int nSegments;

  /// a vector mapping the ids in the patch to local ids for the energy 
  std::vector<int>* locIds;

  Scalar *p2d;

  /// the segment ids for later test w.r.t. occlusions
  int* segImg;

  std::vector<M3> homs; // homographies
  std::vector<P3> vNoms;// view normals, to compute the depth at that pixel quickly

  Scalar occPenalty;

  /// stores pixel ids which are occluded by free pixel - this might lead to potential multiple occlusions which uncover this pixel
  std::set<int> outsideOccluded;

  int maxPatchSize;

  bool genOccData;

  /// storing the occlusions of the non free pixel
  std::list< occPixelList< Scalar > > outerOcclusions;
  /// storing the oclcusions of the free pixel in the case of assignment of 0
  std::list< occPixelList< Scalar > > freeOcclusionsSeg0;
  /// storing the oclcusions of the free pixel in the case of assignment of 0
  std::list< occPixelList< Scalar > > freeOcclusionsSeg1;
};


/// to store the current data scores - to be updated when the solution changes
template<typename Scalar>
struct dataScoresBuffer
{
  typedef std::pair <int,int>  PixSegPair;
  dataScoresBuffer( int _w, int _h, std::vector< Scalar >& _initDataScore ) : w(_w), h(_h)
  {
    initDataScore( _initDataScore );
  };

  ~dataScoresBuffer() {};

  /// initialize the data score 
  void initDataScore( std::vector<Scalar>& _dataScore )
  {
    assert (_dataScore.size() == w*h);
    dataScore.resize( _dataScore.size() );
    std::copy( _dataScore.begin(), _dataScore.end(), dataScore.begin() );
  }

  /// update the data score, based on a score map and the changes: 0 nothing changed 1 changed segment map
  void updateDataScore( std::vector<Scalar>& newScores, std::vector<int>& loc2Glob, std::vector<int>& changes )
  {
    int nPixel = loc2Glob.size();
    for (int i=0;i<loc2Glob.size();i++)
    {
      if ( changes[i]<=0 ) continue;
        dataScore[loc2Glob[i]] = newScores[nPixel+i];
    }
  }

  ////////////////////////////////////////////

  ///width
  int w;
  ///height
  int h;
  /// the data score of the current solution, without considering occlusions
  std::vector< Scalar > dataScore;
};

#if defined(INCLUDE_TEMPLATES) && !defined(__ENERGY_ROTTRANOR__cpp)
#include "OcclusionMappingBuffer.cpp"
#endif

#endif // __OcclusionMapping__h
