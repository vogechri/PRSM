////////////////////////////////////////////////////////////////
// Compute a matrix of scores assigning normal i to segment j //
////////////////////////////////////////////////////////////////

#ifndef __OcclusionExpansionBuffer__h
#define __OcclusionExpansionBuffer__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

using namespace std;
using namespace Math;

template<typename Scalar>
struct oMapBufferLookup
{
public:

  oMapBufferLookup( ) : vecPos(-1), pixPos(-1) {};// oob
  oMapBufferLookup( int _pixPos, int _vecPos) :  pixPos(_pixPos), vecPos(_vecPos)
  {};

  ~oMapBufferLookup() {};
  ///////////////////////////////////////
  int pixPos;//where is it stored/project to
  int vecPos;//position in vector is ..
};
//////////////////////////////////////////////////////////////////////////

/*!
* Note that the buffer is agnostic - it has no knpwledge about homographies or geometry/motion involved
* The implementation is very stupid and very far from good but sufficient here.
*/
template<typename Scalar> class OcclusionExpansionBuffer
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
  typedef Math::VectorT<int, 4>  P4i;

  typedef std::pair <int,int>       Seg_ij;

  OcclusionExpansionBuffer(int w_, int h_, int* _segImg)
    : segImg(_segImg), w(w_), h(h_)
  {}

  ~OcclusionExpansionBuffer(){};

  /// in all homographies, lookup proposal from segment, might need pixel to segment mapping
  void initOccBuffer( std::vector<M3>& homs )
  {
    init( w, h, segImg );
    // precompute the current situation
    int lastSeg = -1;
    M3 H;P3 N;
    {
      for (int i=0;i<w*h;i++)
      {
        int segI = segImg[i];
        if (lastSeg != segI)
        {
          lastSeg = segI;
          assert( lastSeg < homs.size() );
          H = homs [segI];
        }
        /// get the new position
        P3 pixpos = H * P3( int(i/h)+1., int(i%h) +1., 1. );pixpos /= pixpos[2];

        int px ( floor(pixpos[0]-0.5) );
        int py ( floor(pixpos[1]-0.5) );

        if (px>=0 && px < w && py>=0 && py < h)
        {
          int pixId = px*h+py;
          store0[pixId].push_back( i );
          lookup0[i] = oMapBufferLookup<Scalar>( pixId, store0[pixId].size()-1 );
        }
        else
          lookup0[i] = oMapBufferLookup<Scalar>( -1, -1 );
      }
    }
  }

  typename std::vector< int >& get_Occluders( int pixId )
  {
    if ( lookup0[pixId].pixPos > -1 )
      return store0[lookup0[pixId].pixPos];
    else 
      return dummy;
  }
  /// returns a list of all pixel falling onto pixId itself
  typename std::vector< int >& get_Partners( int pixId )
  {
    return store0[ pixId ];
  }

  /// update giving i: global id of pixel, hom its 'new' homography assigned
  bool update( int i, const M3& Hom )
  {
    // first find and remove old entry:
    if (lookup0[i].pixPos >= 0)
    {
      // the list of pixel falling onto the smae as this one:
      std::vector< int >& tmp = store0[lookup0[i].pixPos];
      if (tmp.size() > 1)
      {
        int vp = lookup0[i].vecPos;// position in vec
        tmp[vp] = tmp.back();      // duplicate
        lookup0[i].vecPos = -1;    // erase
        lookup0[i].pixPos = -1;
        lookup0[tmp[vp]].vecPos = vp;// change lookup of other pixel
        tmp.pop_back();            // erase
      }
      else
        tmp.clear();
    }
    //
    // now add entry again: project forward
    P3 pixpos = Hom * P3( int(i/h) +1., int(i%h) +1., 1);pixpos /= pixpos[2];

    // find the pixel occupied:
    int px = floor(pixpos[0]-0.5);
    int py = floor(pixpos[1]-0.5);

    // put into buffer
    if (px>=0 && px < w && py>=0 && py < h)
    {
//    int pixId = px+py*w;
      int pixId = px*h+py;
      store0[pixId].push_back( i );
      lookup0[i] = oMapBufferLookup<Scalar>( pixId, store0[pixId].size()-1 );
    }
    else
      lookup0[i] =oMapBufferLookup<Scalar>(); // redundant

    return true;
  }

private:

  void init( int _w, int _h, int* _segImg )
  { 
    h=_h;
    w=_w;
    //    nSegments = _nSegments;
    segImg    = _segImg;

    store0.clear();
    store0.resize(w*h);
    for (int i=0;i<w*h;i++)
    {
      store0[i].reserve(10);// now that is ugly - this is costing something
    }

    lookup0.clear();lookup0.resize(w*h);
  };

  /////
  int w;
  int h;

  int* segImg;

  /// from the pixel ina the reference image to the pixel it is projected onto
  std::vector<oMapBufferLookup<Scalar> > lookup0;
  
  std::vector< std::vector< int > > store0;

  std::vector< int > dummy;
};

#undef _OccPct_
#endif // __OcclusionMapping__h
