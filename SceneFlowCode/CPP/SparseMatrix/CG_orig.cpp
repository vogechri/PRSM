//=============================================================================
// $TEMPLATE_HEADLINE$
//-----------------------------------------------------------------------------
//
//  CLASS CG_orig - IMPLEMENTATION
//
//-----------------------------------------------------------------------------
//
//=============================================================================

//-***************************************************************************
// for g++ (template defs only compiled through header file)

#if !defined(CG_orig_CC)
#define CG_orig_CC
//-***************************************************************************

//== INCLUDES =================================================================

#include "CG_orig.h"

#include <float.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <limits>

//-----------------------------------------------------------------------------
namespace SparseMath {

//== IMPLEMENTATION ===========================================================


//-----------------------------------------------------------------------------

template <typename Configuration> void
  CG_orig<Configuration>::setup(Element* _x, Element* _b, Element* _r) {
  d_n=d_config.dimension();  // mit A^t * A

  if (_x==0) { d_pX.resize(d_n); d_X=&d_pX.front(); } else d_X=_x; 
  if (_b==0) { d_pB.resize(d_n); d_B=&d_pB.front(); } else d_B=_b; 
  if (_r==0) { d_pR.resize(d_n); d_R=&d_pR.front(); } else d_R=_r; 

  GMU::Math::setComponents(null_,Scalar(0));

  d_tmp.resize(d_n*4, null_);// mit A^t * A 

  d_D     = &d_tmp.front();
  d_Q     = d_D    + d_n;
  d_TMP   = d_Q    + d_n;
  safe_X  = d_TMP  + d_n;
  
}

//-----------------------------------------------------------------------------

template <typename Configuration> void
CG_orig<Configuration>::initialize() {
  ASSERT_EXPR(d_n>0);
  d_nIter=0;
  
  //  mit A^t * A

  d_config.mult(d_TMP,d_X,d_n); // Mat-Mult A*x = d_U
  add(d_R,d_B,d_TMP,Scalar(-1));   // r_0=b-Ax_0

#pragma warning ( push )
#pragma warning(disable:4996)
  std::copy(d_X,d_X+d_n,safe_X); // safe old

  std::copy(d_R,d_R+d_n, d_D);   // d_RR = d_D = d_R
  // residuum
#pragma warning ( pop )

  std::fill(d_TMP,d_TMP+d_n, null_);  // d_tMP = 0

  d_new = d_zero = res_sq = dot(d_R,d_R);
//  mit A^t * A
}

//-----------------------------------------------------------------------------

template <typename Configuration> bool
CG_orig<Configuration>::iterate() {

  // mit A^t*A
  Element dummy;

  std::fill(d_Q,d_Q+d_n,null_);  // noetig ??

  d_config.mult(d_Q,d_D,d_n); // d_Q = A * d_D
  d_gamma  = dot(d_D,d_Q);
  dummy    = (d_new/d_gamma);

  d_alpha  = get_best_of(dummy); 
  // lustiger-weise nicht getestet ob division durch null !!

  add(d_X,d_X,d_D,d_alpha);//-- dx = dx + alpha * dd
  
  if ( (d_nIter % 20)  == 0)
    {
      d_config.mult(d_TMP,d_X,d_n); // Mat-Mult A*x = d_U
      add(d_R,d_B,d_TMP,Scalar(-1));   // r_0=b-Ax_0
      std::fill(d_TMP,d_TMP+d_n,null_);  // d_TMP = 0
    }
  else
    {
      add(d_R,d_R,d_Q, (- d_alpha));   // r_0=b-Ax_0
    }
  
  d_old = d_new;
  d_new = dot(d_R,d_R);
  dummy = (d_new/d_old);
  d_beta = get_best_of(dummy);
  
  add(d_D,d_R,d_D,d_beta);
  
  d_nIter++;

   // mit A^t*A

//   std::cout << "Iteration: "  << d_nIter << "\n";
//   std::cout << "d_old, d_new, d_beta, d_alpha : "  << d_old << " " 
// 	    << d_new << " " << d_beta << " " << d_alpha  << "\n";
//   std::cout << "1st x: " << *d_X << "\n";
  
  bool impr[res_sq.dim()];  // nur gcc
  bool all_improved = true;
  Element res_akt =  d_new;
  bool too_small = true;

  for (int i = 0; (i < res_sq.dim())&&(too_small);i++)
    if (res_sq[i] > std::numeric_limits<Scalar>::epsilon()) too_small = false;

  if (too_small) 
    {
      std::cout << "Solver : tiny residuum:\n";
      for (int i = 0; (i < res_sq.dim())&&(too_small);i++)
	std::cout << "res old " << i << ": " << res_sq[i]  << "\n";

      return false; // und tschuess
    }
  // wenn keine weitere iteration notwendig

  int i;
  //#pragma omp parallel for private(i)
  for (i = 0; i < res_akt.dim();i++)
  {
    impr[i] = false;

    if ((res_sq[i] <= res_akt[i]) || (res_sq[i] < DBL_MIN * 1000 ) )
      all_improved = false;
    else
      impr[i] = true;
      //{ res_sq[i] = res_akt[i];impr[i] = true; };
  };

  // nur wenn alle besser akzeptieren
  if (all_improved) 
    {
#pragma warning ( push )
#pragma warning(disable:4996)
      std::copy (d_X,d_X+d_n,safe_X); // copy all
#pragma warning ( pop )
      for (i = 0; i < res_akt.dim();i++)
	{
	  std::cout << "res old " << i << ": " << res_sq[i]  << "\n";
	  res_sq[i] = res_akt[i];
	}
    }
  else // einzeln verbessern
  {
    bool nonimpr(true);
    int i;
    Element* ix, *ex, *is;
    //#pragma omp parallel for private (i,ix,ex,is)
    for (i = 0; i < res_akt.dim();i++)
      if (impr[i])
      // akt system i
      {
	std::cerr << i << "  improved | ";
	
	ix =d_X, ex =d_X+d_n, is = safe_X;
	//copy whole vektor for system
	while (ix!=ex)
	  {
	    (*is)[i]=(*ix)[i];
	    ++ix;++is;//next element
	  }
	nonimpr = false;
      }
//     if (nonimpr) return false;
//     else std::cerr << std::endl;
  }
  
  if (all_improved){ 
    std::cerr << "all  improved\n";
    for (i = 0; i < res_akt.dim();i++)
      {
	std::cout << "res_akt " << i << ": " << res_akt[i] << "\n";  
      }
    return true;
  }
  // now safe_X stores best solution in every dimension of solution
  //  return false; 
  return true;
  // nur noch false bei residuum < epsilon
}

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace SparseMath
//=============================================================================
//=============================================================================

//-***************************************************************************
#endif // CG_orig_CC
//-***************************************************************************
