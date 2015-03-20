////////////////////////////////////////////////////////////////
// Container for graph cut by truncation //
////////////////////////////////////////////////////////////////

#ifndef __QPBO_Cutting__h
#define __QPBO_Cutting__h

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include <map>
using namespace std;
using namespace Math;

#include "DataContainers.h"
#include "../QPBO-v1.3.src/QPBO.h"

//#define _truncation_test_

//same name same game:
template<typename Scalar> 
class GraphCutContainer
{

//  typedef Graph<int,int,int> GraphType;
  typedef QPBO<int> GraphType;
  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 4>  P4;
  typedef Math::VectorT<Scalar, 3>  P3;
  typedef Math::VectorT<Scalar, 2>  P2;
  typedef Math::VectorT<int, 4>  P4i;
  typedef Math::VectorT<int, 2>  P2i;

  typedef std::pair<int,int> mapkey;
  typedef std::map< mapkey, P4 > edgeMap;


public:

  GraphCutContainer() : nNodes(0), nEdges(0), g(NULL), non_subs(0), num_edges(0), num_nodes(0), mapping(NULL) {};//, lsa0Score(0) {};

  ~GraphCutContainer() {deleteMem();};

  void deleteMem()
  {
    if (g!=NULL) 
      delete g;
    g=NULL;
  }

  void create(int _nNodes, int _nEdges)
  {
    nNodes = _nNodes;
    nEdges = _nEdges;
    deleteMem();
  	g = new GraphType( nNodes, nEdges);
  }

  void setup(int nVars )
  {
    num_nodes = nVars;
  	g->AddNode(nVars);
    unaries.clear();
    unaries.assign(nVars, P2(0.,0.));
  }

  void Reset() {reset();};

  void reset()
  {
    num_nodes = 0;
    non_subs  = 0;
    num_edges = 0;
    g->Reset();
    freeMapping();
  };

  // finializes the submodular graph
  void finalize()
  {
    for (int i=0;i<unaries.size(); i++)
    {
      Scalar minUn = std::min(unaries[i][0], unaries[i][1]);
      g->AddUnaryTerm( i, unaries[i][0]-minUn, unaries[i][1]-minUn );
    }
  };

  int solve()
  {
    g->Solve();
		g->ComputeWeakPersistencies();
    int flow = g -> ComputeTwiceEnergy();
    return flow/2;
  }

  int resolve(bool doProbe, bool doImprove)
  {
    assert(g->GetNodeNum() == num_nodes);
     int maxProbe =3;
     int non_sol =0;
     int impIts =0;
     int allImpIts = 0;

     for (int j = 0; j < g->GetNodeNum();j++)
     {
        int x = g->GetLabel(j);
        if (x<0)
          non_sol++;
     }
     if (non_sol==0) return non_sol;

     // fix
     if (doProbe || doImprove)
     for (int k=0; k< g->GetNodeNum();k++)
     {
        int x = g->GetLabel(k);
        if (x>=0)
          g->SetLabel(k,x);
     }

     if (doProbe)
     {
        QPBO<int>::ProbeOptions options;
        options.C = 200000;
        options.dilation = 1;
        options.weak_persistencies = 1;
        //options.directed_constraints = 0;

        if (mapping != NULL)
          freeMapping();

        mapping  = (int*) malloc( sizeof(int) * g->GetNodeNum() );
        int *mapping2 = (int*) malloc( sizeof(int) * g->GetNodeNum() );

        for (int i = 0; i < g->GetNodeNum(); i++) 
        {
           mapping[i]  = i * 2;
           mapping2[i] = i * 2;
        }

        g->Probe( mapping, options );
        g->ComputeWeakPersistencies();

        for (int i=1;i < maxProbe; i++)
        {
           g->Probe( mapping2, options );
           g->ComputeWeakPersistencies();
           g->MergeMappings( num_nodes, mapping, mapping2 );
        }
        free(mapping2);
     }

     if (doImprove)
     {
       const int maxInner = 10;//10;
       const int maxReps = 10; int imp = 0;
       while (imp < maxInner && impIts < maxReps)
       {
         impIts++;
         for(imp=0; imp<maxInner; imp++ )
           if ( g->Improve() )  {printf("Improved after probing\n");break;};

         allImpIts += imp;
         if (imp == maxInner) break;// out of while
       }
     }

    return non_sol;
  }


  int GetLabel(int j)
  {
    assert( num_nodes>j );

    if (mapping == NULL)
      return g->GetLabel(j);//==1;
    else
    {
       int x  = g->GetLabel(mapping[j]/2);
       if (x>0)
          x     += mapping[j] % 2;
       return x;
    }
  }

  /// if probing on:
  void freeMapping()
  {
    if (mapping != NULL)
      free( mapping );
    mapping = NULL;
  }

  void AddUnaryTerm( int p, const P2& values )
  {
    unaries[p] += values;
  }

  void AddUnaryTerm( int p, Scalar u0, Scalar u1)
  {
    unaries[p] += P2(u0,u1);
  }

  void AddPairwiseTerm( int p, int q, Scalar f00, Scalar f01, Scalar f10, Scalar f11 )
  {
    if (f00+f11 > f10+f01)
      non_subs++;
    num_edges++;

    g->AddPairwiseTerm( p, q, f00, f01, f10, f11 );
  }

  void AddPairwiseTerm( int p, int q, P4& ff )
  {
    if (ff[0]+ff[3] > ff[1]+ff[2])
      non_subs++;
    num_edges++;
    g->AddPairwiseTerm( p, q, ff[0], ff[1], ff[2], ff[3] );
  }

  void StorePairwiseTerm( int p, int q, P4& ff )
  {
    if (ff[0]+ff[3] > ff[1]+ff[2])
      non_subs++;
    num_edges++;

    g->AddPairwiseTerm( p, q, ff[0], ff[1], ff[2], ff[3] );
    /*
    mapkey mykey( std::min(p,q), std::max(p,q) );

    if ( multiEdges.find(mykey) == multiEdges.end() )
      multiEdges[mykey] = ff;
    else
      multiEdges[mykey] += ff;
      */
  }

  void mergeParallelEdges( )
  {
    g->MergeParallelEdges();
  }

  void mergeParallelEdges_Aux( )
  {
    g->MergeParallelEdges();
  }


  void AddPairwiseTermLSAAux( int p, int q, P4& ff )
  {
    AddPairwiseTerm( p, q, ff );return;
  }

  int get_NonSubs()  {return non_subs;};
  int get_NumEdges() {return num_edges;}; 
  int get_NumNodes() {return num_nodes;}; 
  
  int solve_lsa()
  {
    return solve();
  }

private:

  GraphType *g;

  std::vector<P2> unaries;

  int num_nodes;

  int nNodes;
  int nEdges;
  int non_subs;
  int num_edges;

  int *mapping;

  /// in case of multiple edges - these must be merged first before put into GC
//  edgeMap multiEdges;
};
#endif