// Solver stuff

#include <vector>
#include <iostream>
#include "Math/VectorT.h"

//#define DebugOut 1
// Solver stuff
#include "Math/SparseMatrix/RowSparseMatrixMultiplication.h"

using namespace Math;
// tells the matrix solver types how
// to apply the matrix multiplication
template <int N>
struct Sparse_Config {

  typedef double Scalar;  // matrix element
  typedef VectorT<double, N> Element;

  Sparse_Config(SparseMath::RowSparseMatrixMultiplication<Scalar>& _m)
    : d_m(_m)
  {}

  // dimension of right hand side/solution vector
  int dimension() const { return d_m.rows(); }
  unsigned int cols() const { return (d_m.cols()+1); }

  // apply sparse matrix multiplication
  void mult(Element* _y,const Element* _x,int _n) {
    assert(_n==d_m.rows());    
    //std::cout << d_m.rows() << " rows \n";
    d_m.mult(_y,_x);
  }

  void mult_transposed(Element* _y,const Element* _x,int _n) {
    assert(_n==d_m.cols()+1);    
    //std::cout << d_m.rows() << " rows \n";
    d_m.mult_transposed(_y,_x);
  }

private:
  SparseMath::RowSparseMatrixMultiplication<Scalar>& d_m;
};


namespace SparseMath {
  // need this helper function for "fancy" Element types, e.g. 2-vectors
  static void setComponents(Vec3d& _f, double _v) 
  {_f.vectorize(_v);}
  static void setComponents(Vec2d& _f, double _v) 
  {_f.vectorize(_v);}
  static void setComponents(Vec1d& _f, double _v) 
  {_f.vectorize(_v);}

  //    static void setComponents(Vec1d& _f, double _v) 
  //  {_f[0] = _v;}

}

#include "Math/SparseMatrix/BiCGSTAB2.h"
#include "Math/SparseMatrix/cg.h"

// ------------------   Sparse Math stuff   ------------------------
using namespace SparseMath;

template<int N>
class SolveLaplace
{
  typedef typename Sparse_Config<N>::Scalar Scalar;
  typedef typename Sparse_Config<N>::Element Element;

  typedef SparseMath::BiCGSTAB2< Sparse_Config<N> > My_BiCG;
  typedef SparseMath::CG       < Sparse_Config<N> > My_CG;  


public:

  /// width height, dual solutions, on off map, initial solutions
  SolveLaplace(int _w, int _h, Scalar* _wx, int* _map, Scalar* _wX) 
    : w_(_w), h_(_h), wx_(_wx), wy_(NULL), wz_(NULL), map_(_map), wX_(_wX), wY_(NULL), wZ_(NULL), elems_(0), weightsV_(NULL), weightsH_(NULL), invalid_(NULL)
  {
    int sizeCheck[ N == 1 ? 1 : -1 ]; // Fails if N != 1, since an array with -1 elements is illegal
  };
  
  /// width height, dual solutions, on off map, initial solutions
  SolveLaplace(int _w, int _h, Scalar* _wx, Scalar* _wy, int* _map, Scalar* _wX, Scalar* _wY) 
    : w_(_w), h_(_h), wx_(_wx), wy_(_wy), wz_(NULL), map_(_map), wX_(_wX), wY_(_wY), wZ_(NULL), elems_(0), weightsV_(NULL), weightsH_(NULL), invalid_(NULL)
  {
    int sizeCheck[ N == 2 ? 1 : -1 ]; // Fails if N != 2, since an array with -1 elements is illegal
  };

  /// width height, dual solutions, on off map, initial solutions
  SolveLaplace(int _w, int _h, Scalar* _wx, Scalar* _wy, Scalar* _wz, int* _map, Scalar* _wX, Scalar* _wY, Scalar* _wZ) 
    : w_(_w), h_(_h), wx_(_wx), wy_(_wy), wz_(_wz), map_(_map), wX_(_wX), wY_(_wY), wZ_(_wZ), elems_(0), weightsV_(NULL), weightsH_(NULL), invalid_(NULL)
  {
    int sizeCheck[ N == 3 ? 1 : -1 ]; // Fails if N != 3, since an array with -1 elements is illegal
  };

  ~SolveLaplace() {};

  void setup_Matrix( int iter )
  {
    SparseMath::RowSparseMatrixMultiplication<Scalar> A_RowSparse;

    // right hand side
    ids_.clear();
    ids_.resize(w_*h_);

    elems_ = 0;

    Element zeroElem;
    zeroElem.vectorize(0);

    // determine ids:
    for (int j = 0; j < h_; j++)
      for (int i = 0; i < w_; i++)
      {
        int id = j * w_ + i;

        if ((map_[id] == 0) || ((invalid_ != NULL) && (invalid_[id] == 0)))
          ids_[id] = -1;
        else
        {
          ids_[id] = elems_;
          elems_++;
        }
      }

#ifdef DebugOut
      for (int j = 0; j < h_; j++)
      {
        for (int i = 0; i < w_; i++)
        {
          int id = j * w_ + i;
          printf("Pos: (%d,%d)", id, ids_[id] );
        }
        printf("\n");
      }
      printf("#Elems: %d \n", elems_);
#endif

      b_RowSparse_.clear();
      b_RowSparse_.resize( elems_ );
      x_RowSparse_.clear();
      x_RowSparse_.resize( elems_ );
      //    std::fill(matVals.begin(),matVals.end(), zeroElem );

      A_RowSparse.beginSetup();

      for (int j = 0; j < h_; j++)
        for (int i = 0; i < w_; i++)
        {
          int id  = j * w_ + i;
          if ((map_[id] == 0) || ((invalid_ != NULL) && (invalid_[id] == 0)))
            continue;

          x_RowSparse_[ids_[ id ]]  = getElement(id); //Element( wX_[id], wY_[id], wZ_[id] );
          b_RowSparse_[ids_[ id ]]  = zeroElem;
#ifdef DebugOut
          printf (" XVec: %d , %f \n", ids_[ id ],  x_RowSparse_[ids_[ id ]][0] );
#endif
          Scalar Wij(0.0);
          Scalar Wii(0.0);
          A_RowSparse.beginRow();

          // up
          if ((j > 0) && ( (invalid_ == NULL) || (invalid_[id-w_] != 0) ))
          {
            if ( weightsV_ == NULL )
              Wij = 0.25;
            else
              Wij = weightsV_[id-w_];

            Wii += Wij;
            if ( map_[id-w_] == 0) // invalid but yielding a border constraint:
            {
              // add to right hand side:
              b_RowSparse_[ids_ [ id ]] += Wij * getElement(id-w_);//Element( wX_[id-w_], wY_[id-w_], wZ_[id-w_] );
#ifdef DebugOut
              printf(" up ID: %d (%d %d), x: %f", ids_ [ id-w_ ], id, id-w_, wX_[id-w_]);
#endif
            }
            else
            {
              A_RowSparse.set( ids_ [ id-w_ ] , -Wij );
            }
          }
          // down
          if ((j < h_-1) &&( (invalid_ == NULL) || (invalid_[id+w_] != 0)))
          {
            if ( weightsV_ == NULL )
              Wij = 0.25;
            else
              Wij = weightsV_[id];

            //Wij = (weightsU[idv] + weightsD[idv+1]) * (scales2[idv]*x[idv] + scales[idv]);
//            Wij = 0.25;
            Wii += Wij;
            if ( map_[id+w_] == 0) // invlaid  bt yielding a border constraint:
            {
              // add to right hand side:
              b_RowSparse_[ids_ [ id ]] += Wij * getElement(id+w_);//Element( wX_[id+w_], wY_[id+w_], wZ_[id+w_] );
#ifdef DebugOut
              printf(" down ID: %d (%d %d), x: %f", ids_ [ id+w_ ], id, id+w_, wX_[id+w_]);
#endif
            }
            else
            {
              A_RowSparse.set( ids_ [ id+w_ ] , -Wij );
            }
          }
          // left
          if ((i > 0) && ((invalid_ == NULL) || (invalid_[id-1] != 0)))
          {
            //        Wij = (weightsU[idv-N] + weightsD[idv]) * (scales2[idv]*x[idv] + scales[idv]);
//            Wij = 0.25;
            if ( weightsH_ == NULL )
              Wij = 0.25;
            else
              Wij = weightsH_[id-1];

            Wii += Wij;
            if ( map_[id-1] == 0) // invlaid  bt yielding a border constraint:
            {
              // add to right hand side:
              b_RowSparse_[ids_ [ id ]] += Wij * getElement(id-1);//Element( wX_[id-1], wY_[id-1], wZ_[id-1] );
#ifdef DebugOut
              printf(" left ID: %d (%d %d), x: %f", ids_ [ id-1 ], id, id-1, wX_[id-1]);
#endif
            }
            else
            {
              A_RowSparse.set( ids_ [ id-1 ] , -Wij );
            }
          }
          // right
          if ((i < w_-1) && ((invalid_ == NULL) || (invalid_[id+1] != 0)))
          {
            //        Wij = (weightsU[idv] + weightsD[idv+N]) * (scales2[idv]*x[idv] + scales[idv]);
//            Wij = 0.25;

            if ( weightsH_ == NULL )
              Wij = 0.25;
            else
              Wij = weightsH_[id];


            Wii += Wij;
            if ( map_[id+1] == 0) // invalid  bt yielding a border constraint:
            {
              // add to right hand side:
              b_RowSparse_[ids_ [ id ]] += Wij * getElement(id+1);
#ifdef DebugOut
              printf(" right ID: %d (%d %d), x: %f", ids_ [ id+1 ], id, id+1, wX_[id+1]);
#endif
            }
            else
            {
              A_RowSparse.set( ids_ [ id+1 ], -Wij );
            }
          }
#ifdef DebugOut
          printf (" BVec: %d , %f \n", ids_[ id ],  b_RowSparse_[ids_[ id ]][0] );
#endif
          A_RowSparse.set( ids_ [ id ], Wii );
          A_RowSparse.endRow();
        }

        A_RowSparse.endSetup();

#ifdef DebugOut
        printf("End Setup\n");
#endif
        Sparse_Config<N> config (A_RowSparse);

        My_BiCG bcg(config);
        bcg.setup(&x_RowSparse_.front(), &b_RowSparse_.front());
        bcg.initialize();
#ifdef DebugOut
        printf("End Initialize\n");
#endif

        iterate( iter, bcg );
        //  My_CG cg(config);
        //  cg.setup(&x_RowSparse.front(), &b_RowSparse.front());
        //  cg.initialize();


        //  A_RowSparse.print_triv();
  }

  double iterate( int iter, My_BiCG& bcg )
  {
    for (int l = 0;l< iter; l++)
    {
      if (!bcg.iterate()) break;
      //      cg.iterate();
    }
    
#ifdef DebugOut
    printf("End Iterating\n");
#endif
    for (int j = 0; j < h_; j++)
      for (int i = 0; i < w_; i++)
      {
        int id  = j * w_ + i;
        if (map_[id] == 0)
          continue;
#ifdef DebugOut
        printf("Current index is %d", ids_[id]);
        printf("Result : %f", bcg.X()[ids_[id]][0]);
#endif
        setElement(id, bcg);
        //        optVal[j] = cg.X()[j][0];


      }  
#ifdef DebugOut
      printf("\n");
#endif
      return 0.0; // Whatever
  }


void setWeights( Scalar* _wH, Scalar* _wV )
{
  weightsV_ = _wV;
  weightsH_ = _wH;
  //if (_wH != NULL)
  //  printf("weights are set %f %f\n", _wV[0], _wH[1]);
}

void setInvalid( int* _invalid )
{
  invalid_ = _invalid;
  //if (_invalid != NULL)
  //  printf("_invalids are set %d %d\n", _invalid[0], _invalid[1]);
}


  ////////////////////////////////////////////////////////////////
private:

  inline void setElement(int id, My_BiCG& bcg)
  {
    Element e = bcg.X()[ids_[id]];
    wX_[id] = e[0];

    if (N > 1)
      wY_[id] = e[1];

    if (N > 2)
      wZ_[id] = e[2];

#ifdef DebugOut
    if (N == 1)
      printf("Result : %f ", e[0]);
    if (N == 2)
      printf("Result : %f, %f ", e[0], e[1] );
    if (N == 3)
      printf("Result : %f, %f %f ", e[0], e[1], e[2] );
#endif

  }

  inline Element getElement(int id)
  {

#ifdef DebugOut
    if (N == 1)
      printf("Get Element : %f", wX_[0]);
    if (N == 2)
      printf("Get Element : %f", wY_[0]);
    if (N == 3)
      printf("Get Element : %f", wZ_[0]);
#endif

    if (N == 1)
      return Element(wx_[id]);

    if (N == 2)
      return Element(wx_[id], wy_[id]);

    if (N == 3)
      return Element(wx_[id], wy_[id], wz_[id]);
  }

  int w_;
  int h_;

  int elems_;

  /// map from id in the image to id in the solver
  std::vector<int> ids_;

  std::vector<Element> x_RowSparse_;
  std::vector<Element> b_RowSparse_;

  // the right hand sides:
  Scalar* wx_; 
  Scalar* wy_; 
  Scalar* wz_; 

  // the initial values:
  Scalar* wX_; 
  Scalar* wY_; 
  Scalar* wZ_; 

  /// horizontal weights
  Scalar* weightsH_;
  /// vertical weights
  Scalar* weightsV_;

  // points tha should be computed are marked with 1 -> dirichlet points : 0
  int*    map_;

  /// points marked with 0 here are not even considered in the equation system -> vNeuman
  int*    invalid_;
};
