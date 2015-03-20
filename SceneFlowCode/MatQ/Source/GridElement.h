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


#include <map>
#include <vector>
#include <algorithm>
//#include "lapack.h"
#include "Templates/HeapT.h"
#include "Grid_HeapHelper.h"

typedef std::pair<int, double> neighScore;

// got a container of matrices, now a grid element == 
template<typename T= double>
struct grid_element
{ 
  typedef std::map< int, T > nset;

  grid_element() : bestId(-1), bestScore(100000000.0) {};

  ~grid_element() {};

  void setId (int id ) {ownids.clear();ownids.insert( id );};

  std::set<int> & getOwnIds() {return ownids;};

  bool isNeigh( int neighId) { return (ns.end() != ns.find (neighId) ); };
  T getBestScore() {return bestScore;};
  int getBestID()  {return bestId;};

  /// merge two elements
  std::map<int,T>& merge( const grid_element& ge )
  {
    ns.insert(ge.ns.begin(), ge.ns.end());
    ownids.insert(ge.ownids.begin(), ge.ownids.end());
    return ns;
  }

  void remove( int id )
  {
  	typename nset::iterator it = ns.find( id );
    if (it != ns.end()) ns.erase( id );
  }

  void resetBestScore () {bestScore = 100000000.0;};

  void updateBestScore()
  {
	  typename nset::iterator it = ns.begin();
    for ( ;it != ns.end(); it++)
      if ( bestScore > it->second )
      {
        bestScore = it->second;
        bestId    = it->first;
      }
  }

  /// does also update if entry exists
  void setScore( int neighId, T score )
  {
  	typename nset::iterator it = ns.find( neighId );

    if ( it != ns.end() ) // already in map -- unclear if not in map then returns what ?????
    {
      it->second = score;
      if ( it->first == bestId && bestScore < score) // score gets worse
      {
        bestScore = score;
        bestId    = neighId;
        // check whole list 
        updateBestScore();// must check all to find best
      }
      else // not same as current best - or the same and better
        if ( bestScore > score ) // we have a new best or the best is already better so do nothing
        {
          bestId = it->first;
          bestScore = score; 
        }
    }
    else // new insert
    {
      ns.insert ( neighScore( neighId, score ) );
      if ( score < bestScore )
      {
        bestScore = score;
        bestId    = neighId;
      }
    }
  }

  /// remove from neighlist -- update best score stuff
  bool removeNeigh( int neighId )
  {
    if (!isNeigh( neighId )) return false;
  	typename nset::iterator it = ns.find( neighId ); 
    if ( it->first == bestId )
    {
      ns.erase( it );
      resetBestScore ();
      updateBestScore();
    }
    else
      ns.erase( it );
    return true;
  }

  /// 
  T getNeighScore( int neighId )//, T score )
  {
  	typename nset::iterator it = ns.find( neighId ); 
    return (it == ns.end()) ? -1 : it->second;
  }

    std::map< int, T> ns;
    std::set<int> ownids;

    int bestId;
    T bestScore;
};


class grid
{
public:
    typedef HeapStoreEntryT<double>  HeapStoreEntry;
    typedef std::vector< HeapStoreEntry > HeapEntryStore;

    template<int N>
  grid(int nSegments_, const mxArray* Segments_, std::vector< Matmatvecl2<double,N>* >& mat_container, int finalSegs_ = 1000 )
    : Segments(Segments_), //mat_container( mat_container_  ),
      _heapEntryStore(), _heapInterface(&_heapEntryStore), _heap(_heapInterface), finalSegs(finalSegs_)
  {
    // it can be nElements == N*nSegments these are different 'layers'
    // thus edges are to be between nElements % nSegments to all 0..K-1 * nSegments + (nElements % nSegments) with K = nElements / nSegments 
    int nElements = mat_container.size();
    elements.resize(nElements);
    int nSegments = nSegments_;

    int nLayers = nElements/nSegments;
    assert( nElements%nSegments ==0 );
    assert( nLayers > 0 );

    lastScore = -1;
    nElements = mat_container.size();
    //// init
    for (int i=0;i <nElements;i++)
    {
      int sid = i%nSegments;// seg id
      int lid = i/nSegments;// layer id

      int* ids = (int*) mxGetPr(mxGetCell(Segments, sid));
      int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, sid));

      elements[i].setId ( i );

      for (int layer = 0; layer < lid; layer++ )
      {
         int id2 = layer * nSegments + sid;
         double score = mat_container[i]->solve( *mat_container[ id2 ] );
         elements[i]  .setScore( id2, score );
         elements[id2].setScore( i, score );
      }

      // compute stupid score: then put in both grid cells
      for (int j=0;j<nIds;j++)
      {
        int id2 = ids[j];
        id2 = id2 + lid * nSegments;

        if (id2<=i) continue;
        double score = mat_container[i]->solve( *mat_container[ id2 ] );
        elements[i]  .setScore( id2, score );
        elements[id2].setScore( i, score );
      }
    }
    // each 'element' has its score and neigh and neighs in grid

    //// now build heap structure
    _heap.clear();
    _heapEntryStore.clear();
    // cost w*h (larger than maximal) , startSeg : -1 == not touched yet, -2: fixed already 
    _heapEntryStore.resize( nElements, HeapStoreEntry(-1., -1) ); 

    std::vector<char> valid(nElements, 1); // all valid - later more and more get removed

    for (int i=0;i <nElements;i++)
    {
      double score = elements[i].getBestScore();
      HeapStoreEntry& hse = _heapEntryStore[i];
      hse.setCost( score );
      hse.setId( i );
      _heap.insert( HeapEntryD (i) );
      // later instead: // _heap.update( HeapEntryD( nId) );
    }

    //// pop best element, get all neighs, merge 
    // in heap by segment id 
    // pop -> merge list
    /// -> joint neighhood but needs to recompute for all stuff 
    while ( !_heap.empty() && nElements > finalSegs )
    {
        HeapEntryD id = _heap.front();
        _heap.pop_front();
        int cSeg = id.pId();

        HeapStoreEntry& heapE = _heapEntryStore[ cSeg ];

        int id1 = heapE.Id();
        int id2 = elements[ heapE.Id()].getBestID();

        // wait -- if 
        if ( !valid[ id1 ] || !valid[ id2 ]  )
        {
#ifdef _DEBUG
          int test = valid[ id1 ];
          int notGood =0;
          if (test)
            notGood =1;
          grid_element<double>& e1 = elements[id1];
          grid_element<double>& e2 = elements[id2];
#endif
          continue;
        }
        // valid collapse: 
        nElements--;

        // merge , pick lower for this, invalidate other:
        int tmp(-1);
        if (id1 > id2){ tmp = id1;id1 = id2;id2 = tmp;};
        valid[ id2 ] =0;

        grid_element<double>& e1 = elements[id1];
        grid_element<double>& e2 = elements[id2];

        double bestScore = e1.getBestScore();


        if ( bestScore >  1000000.) break; // leave while

        lastScore = bestScore;

        e1.remove( id2 );
        e2.remove( id1 );
        std::map<int,double>& nList = e1.merge( e2 );
        e1.resetBestScore ();

        // also ads in the matstructure
        (*mat_container[id1]) += *mat_container[id2];// update own score, solve then joint - singles - constructor update own score

        // gets updated on the fly: needs new upheap stuff also 
        std::map<int,double>::iterator  it(nList.begin()), it_end(nList.end()); 
        for (;it != it_end;it++)
        {
          int idN = it->first;
          double scoreNew = mat_container[id1]->solve( *mat_container[it->first] );
          it->second = scoreNew; // current element -- needs an update as well - no insert if 
          elements[ idN ].removeNeigh( id2 ); // gone
          elements[ idN ].setScore( id1, scoreNew );

          double score = elements[idN].getBestScore();
          HeapStoreEntry& hse = _heapEntryStore[idN];
          hse.setCost( score );
          _heap.update( HeapEntryD( idN ) );
        }
        // now the element itself, updated internally already (& reference map )
        e1.updateBestScore();
        double score = e1.getBestScore();
        _heapEntryStore[id1].setCost( score );
        if (tmp < 0 ) // no id switch
          _heap.insert( HeapEntryD( id1 ) );
        else
          _heap.update( HeapEntryD( id1 ) );
    };

    clusters.clear();
    for ( int i=0; i < elements.size(); i++)
      if ( valid[ i ] )
        clusters.push_back(  elements[i].getOwnIds() );
    printf("Last Score: %.2f\n", lastScore);
  }

  ~grid() {};

  /// stores id -> joint id
  void copy_seg(double* out_array)
  {
    for (int i=0;i<clusters.size();i++)
    {
      std::set<int>::iterator it(clusters[i].begin()), it_end( clusters[i].end());
      for (;it != it_end;it++)
        out_array[*it] = i;
    }
  }


private:

//  std::vector< Matmatvecl2_d6* >& mat_container;
  std::vector< grid_element<double> > elements;
  const mxArray* Segments;

  /// the amount of segmetns in the images
  int nElements;

    /// stores relevant information about HE collapse
  HeapEntryStore _heapEntryStore;

  /// the heap-Interface structure
  HeapInterfaceD < HeapEntryD, HeapEntryStore > _heapInterface;

  /// the heap for selecting the element with the least cost
  Utils::HeapT< HeapEntryD, HeapInterfaceD < HeapEntryD, HeapEntryStore > > _heap;

  std::vector< std::set<int> > clusters;

  double lastScore;

  int finalSegs;
};