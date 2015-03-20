//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//  CLASS RowSparseMatrixMultiplication
//
//-----------------------------------------------------------------------------
// 
// $Id: RowSparseMatrixMultiplication.h,v 1.0 2007/08/27 c.vogel Exp $
// $Date: 2007//08//27 
// $Revision: 1.0
// Created by $Author: c.vogel $
// 
//=============================================================================

#ifndef ROWSPARSEMATRIXMULTIPLICATION_HH
#define ROWSPARSEMATRIXMULTIPLICATION_HH


//== INCLUDES =================================================================

#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <algorithm>

#include "Math/VectorT.h"


//== CLASS DEFINITION =========================================================
namespace SparseMath {
	      
/** \class RowSparseMatrixMultiplication RowSparseMatrixMultiplication.hh
    \brief Sparse matrix (row-wise storage) with mult() operation.
    \ingroup gmu_math
    
    This class represents a sparse matrix \c M. The only operation defined
    is the multiplication with a vector \c x, i.e. <tt>y=Mx</tt>.
   
    The matrix is setup row by row using nested calls like
    \code
    M.beginSetup();
     M.beginRow(); M.set(0,1); M.endRow();
     M.beginRow(); M.set(1,1); M.endRow();
     M.beginRow(); M.set(2,1); M.endRow();
    M.endSetup();
    \endcode
    This would setup a 3x3 identity matrix. The matrix cannot be modified 
    after that, unless it is completely setup again by beginSetup(), ...

    The mult() method multiplies \c M a Element vector \c x and writes
    the result to the Element vector \c _y.
    
    \a Notes:
    \arg There is no means of accessing specific elements of the matrix!
    \arg Dimension: rows() is the number of calls to beginRow(),
    cols() is the maximum number of columns provided to set()! The vectors
    \c x and \c y (mult()) must have  at least dimension cols()!
   
    \sa SparseMath::BiCGSTAB2 SparseMath::CG

    $Revision: 1.0
    Created by $Author: c.vogel $

*/
  
#define ASSERT_EXPR assert

template <typename _Element> class RowSparseMatrixMultiplication {
public:
  ///<i></i>
  typedef _Element Element;

  //
  // CONSTRUCTORS
  //
  
  /// Default constructor. Calls cleear().
  RowSparseMatrixMultiplication() { d_curRowIdx=d_maxCol=0; }

  /// Destructor.
  ~RowSparseMatrixMultiplication() {}

    
  //
  // METHODS
  //

  /// reset matrix storage
  void beginSetup() { 
    d_m.clear();
    d_rows.clear();
    d_curRowIdx=d_maxCol=0;
  }
  /// start setup of a row
  void beginRow() {
    d_rows.push_back(d_curRowIdx);
  }
  /// set element in current row
  void set(int _column,const Element& _element) {
    MElt e; e.e=_element; e.c=_column;
    d_m.push_back(e);
    if (_column>d_maxCol) d_maxCol=_column;
    ++d_curRowIdx;
  }
  /// end setup of current row
  void endRow() {
    CmpMElt cmp;
    std::sort(d_m.begin()+d_rows.back(),d_m.end(),cmp);
  }
  /// end setup of matrix
  void endSetup() { d_rows.push_back(d_curRowIdx); } 
  // => d_rows.size()==rows()+1 !


  /// returns number of rows
  int rows() const { return d_rows.size()-1; }
  /// returns number of columns
  int cols() const { return d_maxCol; }

  /** Matrix multiplication <tt>_Y=A x _X</tt>.
      Make sure that the range of the iterators is at least rows()!
      \param _Y iterator, result is stored here
      \param _X \c const iterator, multiply matrix with this vector
   */

  void print()
  {
    typename std::vector<MElt>::const_iterator mi=d_m.begin();
    std::vector<int>::const_iterator  ri=d_rows.begin();
    
    int row = 0;
    while (ri!=d_rows.end()-1) { 
      typename std::vector<MElt>::const_iterator me=mi+(*(ri+1)- *ri);
      //::SparseMath::setComponents(*_y,Scalar(0));   
      // --- g++ does not accept this (nor *_y *= Scalar(0))
      
      while (mi!=me) 
	{
	  std::cout << "| e:" << row << "-" << mi->c << " " <<  mi->e; 
	  ++mi;
	}
      std::cout << "\n";
      row++;
      ++ri;
    }
  }

  void print_triv()
  {
    std::cout << "Anzahl Elemente:" <<  d_m.size() << "\n";
  };

#ifdef _OPENMP
  template <typename Iterator,typename ConstIterator>
    void mult(Iterator _y, ConstIterator _x) const { 

//    int i,j,k,p,q;
    int i,j,s = d_rows.size()-1;
//    p=0;q=0;
#pragma omp parallel for private (j,i)
    for (i=0;i<s;i++)
    {
      //q = omp_get_thread_num();
      //p++;
      //if (p % 100 == 0) std::cerr << p << " : p  " << q << "q\n";
      // works with export OMP_NUM_THREADS = x (krag:2)
      // check env | grep "OMP"
      SparseMath::setComponents(*(_y+i),0);
      for (j = d_rows[i];j <d_rows[i+1];j++)
	*(_y+i) += d_m[j].e * *(_x + d_m[j].c);//x[d_m[j].c];
    }
  }
#else
  template <typename Iterator,typename ConstIterator>
  void mult(Iterator _y, ConstIterator _x) const { 

    typename std::vector<MElt>::const_iterator mi=d_m.begin();
    std::vector<int>::const_iterator  ri=d_rows.begin();
    while (ri!=d_rows.end()-1) { 
      typename std::vector<MElt>::const_iterator me=mi+(*(ri+1)- *ri);
      //::SparseMath::setComponents(*_y,Scalar(0));   
      // --- g++ does not accept this (nor *_y *= Scalar(0))
      if (mi!=me) {
	*_y = mi->e * *(_x+mi->c); 
	// first iteration --- split due to g++ problems
	++mi;	
	
	while (mi!=me) {
	  *_y += mi->e * *(_x+mi->c);
	  ++mi;
	} 
      } 
      ++ri; ++_y;
    }
  }
#endif

  template <typename Iterator,typename ConstIterator>
  void mult_transposed(Iterator _y,ConstIterator _x) const { 
    
    typename std::vector<MElt>::const_iterator mi=d_m.begin();
    std::vector<int>::const_iterator  ri=d_rows.begin();

    // ri row iterator ---   mi matrix iterator
    int cols = 0;
    while (ri!=d_rows.end()-1) {
      typename std::vector<MElt>::const_iterator me=mi+(*(ri+1)- *ri);
      //::SparseMath::setComponents(*_y,Scalar(0));
      //--- g++ does not accept this (nor *_y *= Scalar(0))

      while (mi!=me) // ende erreicht bei iteration
	{
	  *(_y+mi->c) += mi->e * *(_x+cols);  
	  // mi->e : matrix element 
	  //first iteration --- split due to g++ problems
	  ++mi;	
	} 
      ++ri; cols++;
    }
  }

  /*
template <typename Iterator,typename ConstIterator>
  void mult_transposed_normal(Iterator _y,ConstIterator _x) const { 
    typename std::vector<MElt>::const_iterator mi=d_m.begin();
    std::vector<int>::const_iterator  ri=d_rows.begin();

    // ri row iterator ---   mi matrix iterator
    unsigned int cols = 0;
    while (ri!=d_rows.end()-1) {
      typename std::vector<MElt>::const_iterator me=mi+(*(ri+1)- *ri);
      //::SparseMath::setComponents(*_y,Scalar(0));
      //--- g++ does not accept this (nor *_y *= Scalar(0))

      while (mi!=me) // ende erreicht bei iteration
	{
	  *(_y+mi->c) += mi->e * *(_x+cols);

	  // mi->e : matrix element 
	  //first iteration --- split due to g++ problems
	  ++mi;	
	} 
      ++ri; cols++;
    }
  }
  */

  /// Matrix multiplication <tt>_Y=A x _X</tt>. Provided for convenience.
  template <typename ElementY,typename ElementX>
  void mult(std::vector<ElementY>& _y,const std::vector<ElementX>& _x) const {
    if (int (_y.size())!=( rows() )) std::cout << _y.size() <<" :y-size, " << rows() <<":rows\n";
    ASSERT_EXPR(int (_y.size())==( rows() ));
    if (int(_x.size())!=( (cols()+1) )) std::cout << _x.size() << " x-size\n";
    ASSERT_EXPR(int (_x.size())==( cols() +1));
    mult(_y.begin(),_x.begin());
  }

  template <typename ElementY,typename ElementX>
    void mult_transposed(std::vector<ElementY>& _y,
			 const std::vector<ElementX>& _x) const 
    {
      if (int(_x.size())!=rows() )
	  std::cout << _x.size() << "  " << rows() <<" y-size,rows";
      ASSERT_EXPR(int (_x.size())==(rows()));
      if (int(_y.size())!= (cols()+1))
	  std::cout << _y.size() << " " << cols()+1 << " x-size,cols\n";
      ASSERT_EXPR(int(_y.size())==(cols()+1));


      //ElementY null; SparseMath::setComponents(null,Scalar(0));
      //std::fill(d_U,d_U+d_n,null);  // d_U = 0

      mult_transposed( _y.begin(), _x.begin() );  
    }

protected:
  /// matrix element
  struct MElt {
    Element e; //!< value
    int     c; //!< column index
  };
  /// compare MElt objects
  struct CmpMElt { 
    bool operator()(const MElt& _a,const MElt& _b) const {
      return _a.c<_b.c;
    }
  };

  int               d_curRowIdx;
  int               d_maxCol;
  std::vector<MElt> d_m;     //!< matrix storage
  std::vector<int>  d_rows;  //!< row indices
};

//=============================================================================
} //namespace SparseMath
//=============================================================================
// include templates in header file (g++ only)
# if defined(__GNUG__) && !defined(ROWSPARSEMATRIXMULTIPLICATION_CC)
#  define ROWSPARSEMATRIXMULTIPLICATION_TEMPLATE
#  include "RowSparseMatrixMultiplication.cpp"
# endif
//=============================================================================
#endif // ROWSPARSEMATRIXMULTIPLICATION_HH defined

