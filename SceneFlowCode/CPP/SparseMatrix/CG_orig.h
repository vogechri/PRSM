//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//  CLASS CG_orig
//
//-----------------------------------------------------------------------------
//=============================================================================

#ifndef CG_orig_HH
#define CG_orig_HH


//== INCLUDES =================================================================

#include <vector>

//== CLASS DEFINITION =========================================================
namespace SparseMath {
	      
/** 

/** \class CG_orig.h CG_orig.cpp
    \brief Conjugate gradient solver CG, the given matrix must be symmetrix positive definite.
    \ingroup sparse_math
    
    Implementation of the CG algorithm taken from
    <em>Van der Vorst, Chan. Linear System Solvers: Sparse Iterative Methods.
    1997</em>.

    Solves the (sparse) linear system <tt> Ax=b </tt>, where \c A is given
    in form of a \c Configuration object. It defines the following:
    \arg The \c Configuration object most importantly defines/performs the
    the linear operator <tt>y=Ax</tt>, i.e. the matrix multiplication.
    This is done in the form
    \arg \t Configuration::Scalar, \c Configuration::Element<br>
    These two types define scalar numbers and the elements of the \c y and
    \c x vectors.
    \arg <tt>int Configuration::dimension() const</tt><br>
    returns the dimension (#rows) of the \c y and \c x vectors, this 
    is evaluated by setup().
    \arg <tt>void Configuration::mult(Element* y,const Element* x,int n)</tt><br>
    Where \c y and \c x are the respective vectors and \c n ist their
    dimension (#rows).
    
    \a Hint: SparseMath::RowSparseMatrixMultiplication can be used as matrix
    representation and for performing the multiplication.

    \a Note: Currently no preconditioning is implemented!
    

    \sa SparseMath, math::cg, SparseMath::RowSparseMatrixMultiplication
   

   $Date: 2007//08//27 
   $Revision: 1.0
   $Author: c.vogel $
*/

#ifdef _WIN32
#define FINITE _finite
#else
#define FINITE finite
#endif
  
template <typename _Configuration> class CG_orig {
public:
  /// <i></i>
  typedef _Configuration Configuration;
  /// <i></i>
  typedef typename Configuration::Scalar Scalar;
  /// <i></i>
  typedef typename Configuration::Element Element;
   
  //
  // CONSTRUCTORS
  //
  
  /// Constructor.
  CG_orig(Configuration& _config) : d_n(-1), d_config(_config), d_nIter(0) {
    d_X=d_B=d_R=0; 
  }
 
  /// Destructor.
  ~CG_orig() {}

    
  //
  // METHODS
  //
  
  /// returns Configuration (including types and matrix multiplication)
  Configuration& config() { return d_config; }
  /// same but \c const
  const Configuration& config() const { return d_config; }

  /** Setup memory for solver.
      Initialization of _x and _B should be done afterwards (X(),B()).
      If the vectors are provided by the user then \a he has to make sure that
      their size fits!
      \param _x initial value/guess
      \param _b right hand side
      \param _r residuum
   */
  void setup(Element* _x=0, Element* _b=0, Element* _r=0);
  
  void set_x_new(Element* _x)
  {
    d_X=_x;
    safe_X = _x;
  };

  /// get solution vector (use for initial guess)
  Element* X() { return safe_X; }
  /// same but \c const
  const Element* X() const { return safe_X; }

  /// get right hand side B (use for initialization)
  Element* B() { return d_B; }
  /// same but \c const
  const Element* B() const { return d_B; }

  /// get residuum
  const Element* R() const { return d_R; }
  
  /// calculate squared Euklidean norm of R() (componentwise in Element)
  Element sqrnormR() const { return dot(d_R,d_R); }
  
  /** Initialize solver.
      X() and B() must be set, and config() must be initialized.
      Call just before iteration().
   */
  void initialize();

  /// do one single iteration
  bool iterate();

  /// return number of iterations since last call to initialize()
  int nIterations() const { return d_nIter; }


private:
  /// Copy constructor. Never used!
  CG_orig(const CG_orig&);
  /// Assignment operator. Never used!
  CG_orig& operator=(const CG_orig&);

  /// <i></i>
  typedef std::vector<Element> VECTOR;

  /// dot product, pointer to first element in VECTOR
  Element dot(const Element* _a,const Element* _b) const {
    const Element* a=_a, *ae=_a+d_n, *b=_b;
    Element sum; GMU::Math::setComponents(sum,Scalar(0));
    while (a!=ae) {
      sum+=(*a)*(*b);
      ++a; ++b;
    }
    return sum;
  }
  
  /// dot product, pointer to first element in VECTOR
  Element dot_m(const Element* _a,const Element* _b) const {
    const Element* a=_a, *ae=_a+d_m, *b=_b;
    Element sum; GMU::Math::setComponents(sum,Scalar(0));
    while (a!=ae) {
      sum+=(*a)*(*b);
      ++a; ++b;
    }
    return sum;
  }

  /// VECTOR addition, _c=_a + _b*_v, returns first ELEMNT in _C
  template <typename WEIGHT>
  Element* add(Element* _c,const Element* _a,const Element* _b,WEIGHT _v) const {
    const Element* a=_a, *ae=_a+d_n, *b=_b;
    Element* c=_c;
    while (a!=ae) {
      *c=(*a)+(*b)*_v;
      ++a; ++b; ++c;
    }
    return _c;
  }

  //template <typename Configuration> 
  Element get_best_of(Element &e) const{
    Element e_;
    for (int i = 0; i < e.dim();i++)
      FINITE(e[i]) ? e_[i] = e[i] : e_[i] = Scalar (0.0);
    return e_;
    //return true;
  }


  //Element get_best_of_e(Element &e);

  int      d_n;   //!< dimension
  int      d_m;   //!< rows
  VECTOR   d_tmp; //!< temporary vectors
  VECTOR   d_pX;  //!< private X
  VECTOR   d_pR;  //!< private R 
  VECTOR   d_pB;  //!< private B
  Element* d_X,* d_B, *d_R, *d_TMP, *d_TMP_m; //!< pointer to VECTOR

  //Element* d_RR, *d_RH, *d_S, *d_T, *d_U, *d_V, *d_W; //!< pointer into d_tmp
  Element null_;
  Element* d_D, *d_Q;//, *d_RR;
  Element  d_new, d_old, d_alpha, d_beta, d_zero, d_gamma;
  Configuration& d_config; //!< types and matrix multiplication
  int d_nIter; //!< number of iterations
  Element* safe_X;
  Element res_sq;
  Scalar  full_res; 
};

//=============================================================================
} //namespace SparseMath
//=============================================================================
// include templates in header file (g++ only)
/*
# if defined(__GNUG__) && !defined(GMU_BICGSTAB2_CC)
#  include "BiCGSTAB2.cc"
#  define GMU_BICGSTAB2_TEMPLATE
# endif
*/
# if !defined(CG_orig_CC)
#  include "CG_orig.cpp"
# endif

//=============================================================================
#endif // CG_orig_HH defined

