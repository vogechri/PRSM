//=============================================================================
// $TEMPLATE_HEADLINE$
//-----------------------------------------------------------------------------
//
//  CLASS BiCGSTAB2 - IMPLEMENTATION
//
//-----------------------------------------------------------------------------
// 
// $Id: BiCGSTAB2.hh,v 1.2 2002/11/08 09:52:53 cvogel Exp $
// $Date: 2007//08//27 
// $Revision: 1.0
// Created by $Author: c.vogel $
//
//
//=============================================================================

//-***************************************************************************
// for g++ (template defs only compiled through header file)

// Set to print out residuum information (first find out how many iteration
// are necessary, then dsiable printing (costs time) )
//#define BiCG_PrintResiduum

#if !defined(BICGSTAB2_CC)
#define BICGSTAB2_CC
//-***************************************************************************

//== INCLUDES =================================================================

#include "BiCGSTAB2.h"

#include <float.h>
#include <math.h>
#include <algorithm>
#include <limits>

//-----------------------------------------------------------------------------
namespace SparseMath {

//== IMPLEMENTATION ===========================================================


//-----------------------------------------------------------------------------

template <typename Configuration> void
  BiCGSTAB2<Configuration>::setup(Element* _x, Element* _b, Element* _r) {
   d_n=d_config.dimension();
   if (_x==0) { d_pX.resize(d_n); d_X=&d_pX.front(); } else d_X=_x; 
   if (_b==0) { d_pB.resize(d_n); d_B=&d_pB.front(); } else d_B=_b; 
   if (_r==0) { d_pR.resize(d_n); d_R=&d_pR.front(); } else d_R=_r; 
   d_tmp.resize(d_n*8);
   d_RR=&d_tmp.front();
   d_RH=d_RR+d_n;
   d_S= d_RH+d_n;
   d_T=d_S+d_n;
   d_U=d_T+d_n;
   d_V=d_U+d_n;
   d_W=d_V+d_n;
   safe_X = d_W+d_n;
}

//-----------------------------------------------------------------------------

template <typename Configuration> void
BiCGSTAB2<Configuration>::initialize() {
  assert(d_n>0);
  d_nIter=0;
  

  d_config.mult(d_U,d_X,d_n); // Mat-Mult A*x = d_U
  add(d_R,d_B,d_U,Scalar(-1));   // r_0=b-Ax_0

#pragma warning ( push )
#pragma warning(disable:4996)
#define _SCL_SECURE_NO_WARNINGS
  std::copy(d_X,d_X+d_n,safe_X); // safe old

  std::copy(d_R,d_R+d_n, d_RR);   // \hat{r}=r=r_0  
  std::copy(d_R,d_R+d_n, d_RH);  // d_RR = d_RH = d_R
#undef _SCL_SECURE_NO_WARNINGS
#pragma warning ( pop )

  //ASSERT_EXPR(dot(d_RH,d_RR)>SparseMath::epsilonValue(Scalar(4)));
  
  SparseMath::setComponents(d_rho0,Scalar(1));

  Element null; SparseMath::setComponents(null,Scalar(0));
  std::fill(d_U,d_U+d_n,null);  // d_U = 0
  SparseMath::setComponents(d_alpha,Scalar(0));
  SparseMath::setComponents(d_omega2,Scalar(1));

  res_sq = dot(d_R,d_R);

}

//-----------------------------------------------------------------------------

template <typename Configuration> bool
BiCGSTAB2<Configuration>::iterate() {

  Element dummy;

  d_rho0=d_omega2*d_rho0*Scalar(-1);

  // even Bi-CG step
  d_rho1=dot(d_RH,d_R);
  dummy = (d_alpha*d_rho1/d_rho0);
  d_beta=get_best_of(dummy); 
  d_rho0=d_rho1; 
  
  add(d_U,d_R,d_U,d_beta*Scalar(-1));

  d_config.mult(d_V,d_U,d_n); // d_V = A * d_U  //start : d_U = d_R + 0 * d_U
  d_gamma=dot(d_V,d_RH);
  dummy = (d_rho0/d_gamma);
  d_alpha = get_best_of(dummy); 
  // lustiger-weise nicht getestet ob division durch null !!
  //  daher meine Aufgabe !!
  


  add(d_RR,d_R,d_V,d_alpha*Scalar(-1));
  d_config.mult(d_S,d_RR,d_n);
  add(d_X,d_X,d_U,d_alpha);//--
  
  // odd Bi-CG step
  d_rho1=dot(d_RH,d_S); 
  
  dummy = d_alpha*d_rho1/d_rho0 ;
  d_beta = get_best_of(dummy);
  d_rho0=d_rho1;

  add(d_V,d_S,d_V,d_beta*Scalar(-1));
  d_config.mult(d_W,d_V,d_n);

  d_gamma=dot(d_W,d_RH);
  dummy = d_rho0/d_gamma;
  d_alpha=get_best_of(dummy);

  add(d_U,d_RR,d_U,d_beta*Scalar(-1)); //--
  add(d_RR,d_RR,d_V,d_alpha*Scalar(-1)); //--
  add(d_S,d_S,d_W,d_alpha*Scalar(-1));
  d_config.mult(d_T,d_S,d_n);
  
  // GCG(2)-part
  d_omega1=dot(d_RR,d_S);
  Element mu=dot(d_S,d_S), nu=dot(d_S,d_T), tau=dot(d_T,d_T);
  
  d_omega2=dot(d_RR,d_T);
  dummy = nu*nu/mu;
  tau-=get_best_of(dummy); 
  dummy = (d_omega2-nu*d_omega1/mu)/tau;
  d_omega2= get_best_of (dummy);

  dummy = (d_omega1-nu*d_omega2)/mu;
  d_omega1= get_best_of( dummy );


  Element* ix=d_X, *ex=d_X+d_n, *ir=d_R, *irr=d_RR, *is=d_S, *it=d_T, *iu=d_U;
  while (ix!=ex) {
    (*ix) += (*irr)*d_omega1 + (*is)*d_omega2 + (*iu)*d_alpha; // new Vec X
    (*ir) =  (*irr) - (*is)*d_omega1 - (*it)*d_omega2;  // new residuum
    ++ix; ++ir; ++irr; ++is; ++it; ++iu;
  }

  Element* iv=d_V, *iw=d_W;
  for (iu=d_U,ex=d_U+d_n;iu!=ex;++iu,++iv,++iw) 
    (*iu) -= (*iv)*d_omega1 + (*iw)*d_omega2;

  ++d_nIter;

  std::vector<bool> impr(res_sq.dim(), false);
  bool all_improved = true;
  Element res_akt =  dot(d_R,d_R);
  bool too_small = true;

  for (int i = 0; (i < res_sq.dim())&&(too_small);i++)
    if (res_sq[i] > std::numeric_limits<Scalar>::epsilon()) too_small = false;

  if (too_small) 
    {
#ifdef BiCG_PrintResiduum
      std::cout << "\nSolver : tiny residuum:\n";
      for (int i = 0; (i < res_sq.dim())&&(too_small);i++)
	std::cout << "res old " << i << ": " << res_sq[i]  << "\n";
#endif
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
#ifdef BiCG_PrintResiduum
	  std::cout << "\rres old " << i << ": " << res_sq[i]  << " | ";
#endif
	  res_sq[i] = res_akt[i];
	}
    }
  else // einzeln verbessern
  {
    bool nonimpr(true);
    int i;
    //#pragma omp parallel for private (i,ix,ex,is)
    for (i = 0; i < res_akt.dim();i++)
      if (impr[i])
      // akt system i
      {
//	std::cerr << i << "  improved | ";
	
	ix =d_X, ex =d_X+d_n, is = safe_X;
	//copy whole vektor for system
	while (ix!=ex)
	  {
	    (*is)[i]=(*ix)[i];
	    ++ix;++is;//next element
	  }
	nonimpr = false;
      }
//    if (nonimpr) return false;
    if (nonimpr) return true;

    else std::cerr << std::endl;
  }
  
  if (all_improved){ 
//    std::cerr << "all  improved\n";
#ifdef BiCG_PrintResiduum
    for (i = 0; i < res_akt.dim();i++)
      {
	std::cout << "res_akt " << i << ": " << res_akt[i] << "\t\t\t";  
      }
#endif
    return true;
  }
  // now safe_X stores best solution in every dimension of solution
  return true; 
}

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace SparseMath
//=============================================================================
//=============================================================================

//-***************************************************************************
#endif // BICGSTAB2_CC
//-***************************************************************************
