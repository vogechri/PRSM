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
#include <vector>
#include <algorithm>
#include "lapack.h"

/// apears to be worse - is it ? or just roughly the same
//#define infogain
#define _thresh_pix_ 0.1
// sucks:
//#define _newVariant_
#define _maxDeviation_ 1.

void print_matrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.4f", a[i+j*lda] );
                printf( "\n" );
        }
}

void print_vector_norm( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        double norm;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) {
                norm = 0.0;
                for( i = 0; i < m; i++ ) norm += a[i+j*lda] * a[i+j*lda];
                printf( " %6.2f", norm );
        }
        printf( "\n" );
}

void print_err_norm( int n, double* b1, double* b2, double* bb ) {
        int j;
        double norm;
        printf( "\n %s\n", "error" );
        norm = 0.0;
        for( j = 0; j < n; j++ ) {
           norm += b1[j] * b2[j];
        }
        printf( " %6.2f", bb[0] - norm );
        printf( "\n" );
}


template <class T, int N>
class Matmatvecl2
{
	public:

	// well the matrix and other stuff is local/own property
	// i cannot solve later with k rhs == neighs if i want to sum matrices etc, so
	static const int rows = N;
	static const int cols = N;	
	
	Matmatvecl2( const double* M_, const double* V_, const double* b_ )
    : M(N*N,0), M_h(N*N,0),  V(N,0), V_h(N*2,0), b(*b_), solution(N,0), own_solution(N,0), pps(1), ppsOn(false), ownScore(0), avScore (0)
	{
		//copy: b already done
		std::copy ( M_, M_+N*N , M.begin());
		std::copy ( V_, V_+N   , V.begin());
    redo_ownscore( );
	}
	
	Matmatvecl2( const double* M_, const double* V_, const double* b_, const double* pps_)
    : M(N*N,0), M_h(N*N,0),  V(N,0), V_h(N*2,0), b(*b_), solution(N,0), own_solution(N,0),pps(*pps_), ppsOn(true), ownScore(0), avScore (0)
	{
		//copy: b already done
		std::copy ( M_, M_+N*N , M.begin());
		std::copy ( V_, V_+N   , V.begin());
    redo_ownscore( );
	}

	Matmatvecl2( ): M(N*N,0), M_h(N*N,0),  V(N,0), V(N*2,0), b(0), solution(N,0), own_solution(N,0), ownScore(0), avScore (0)
	{};	

    /// copy constructor
	Matmatvecl2 ( const Matmatvecl2<T,N>& mmv_): M(N*N,0), M_h(N*N,0),  V(N,0), V_h(N*2,0), b(0), solution(N,0), own_solution(N,0), ownScore(0), avScore (0)
	{
		operator=(mmv_);
	}

  ~Matmatvecl2( ){};

  void usePPS( bool onOff) {ppsOn = onOff;}

  T getAVScore () {return avScore;};

	/// assignment from other matrix type
	Matmatvecl2<T,N>& operator=( const Matmatvecl2<T, N>& mmv_)
	{
		// copy
		std::copy ( mmv_.M.begin(), mmv_.M.end() , M.begin());
		std::copy ( mmv_.V.begin(), mmv_.V.end() , V.begin());
		b = mmv_.b;
		std::copy ( mmv_.solution.begin(), mmv_.solution.end() , solution.begin());
    pps = mmv_.pps;
		ownScore = mmv_.ownScore;
		std::copy ( mmv_.own_solution.begin(), mmv_.own_solution.end() , own_solution.begin());
		return *this;
	}

	Matmatvecl2<T,N> operator+( const Matmatvecl2<T,N>& mmv_) const 
	{ 
		Matmatvecl2<T,N> mmv (mmv_);
		
		std::transform(  mmv.M.begin(), mmv.M.end(), M.begin(), mmv.M.begin(), std::plus<T>()  );
		std::transform(  mmv.V.begin(), mmv.V.end(), V.begin(), mmv.V.begin(), std::plus<T>()  );
		mmv.b   += b;
		mmv.pps += pps;
    mmv.redo_ownscore( );
		return mmv;
	}
		
	Matmatvecl2<T,N>& operator+=( const Matmatvecl2<T,N>& mmv_)
	{ 
//		Matmatvecl2<T,N> mmv (mmv_);
		
		std::transform(  M.begin(), M.end(), mmv_.M.begin(), M.begin(), std::plus<T>()  );
		std::transform(  V.begin(), V.end(), mmv_.V.begin(), V.begin(), std::plus<T>()  );
		b   += mmv_.b;
    pps += mmv_.pps;
    redo_ownscore( );
		return *this;
	}

  /// finds own solution
	void solve(  ) // should be specialized ! float, dounle , etc.
	{
		std::copy ( M.begin(), M.end(), M_h.begin());
		std::copy ( V.begin(), V.end(), V_h.begin());
		solve_p();
	}
	
	///	solves the equations by combining both problems and puts solution, but does not store combination
#ifdef _newVariant_
	T solve_old( const Matmatvecl2<T,N>& mmv_ )
#else
	T solve( const Matmatvecl2<T,N>& mmv_ )
#endif
	{
		std::transform(  M.begin(), M.end(), mmv_.M.begin(), M_h.begin(), std::plus<T>()  );
		std::transform(  V.begin(), V.end(), mmv_.V.begin(), V_h.begin(), std::plus<T>()  );
		solve_p(); // stored in solution

    // reconstruct . dgels overwrites it
		std::transform(  V.begin(), V.end(), mmv_.V.begin(), V_h.begin(), std::plus<T>()  );

	  T norm(0);
	  for( int j = 0; j < N; j++ )
      norm += solution[j] * V_h[j];

    if ( ppsOn ) // average per pixel
      avScore = (mmv_.b + b - norm) / (pps + mmv_.pps);

//    if (!ppsOn )
    //avscore not updated == 0 -> last part not used
  		return (mmv_.b + b - norm) - ownScore - mmv_.ownScore + (avScore>_thresh_pix_) * 1000000.;
 //   else
 // 		return (mmv_.b + b - norm)/ sqrt(pps + mmv_.pps) - ownScore - mmv_.ownScore;

//  		return (mmv_.b + b - norm)/ (pps + mmv_.pps);//fails
	}	

  ///	solves the equations by trying the solution of each other at own problem
#ifdef _newVariant_
	T solve( const Matmatvecl2<T,N>& mmv_ )
#else
	T solve_new( const Matmatvecl2<T,N>& mmv_ )
#endif
	{
    T score1 = mmv_.score( own_solution );
    T score2 = score( mmv_.own_solution );

    T md = std::max( score1/mmv_.pps, score2/pps);
    return md + (md > _maxDeviation_) * 1000000.;

/*
    if ( ppsOn ) 
      avScore = (mmv_.b + b - norm) / (pps + mmv_.pps);

//    if (!ppsOn )
  		return (mmv_.b + b - norm) - ownScore - mmv_.ownScore + (avScore>_thresh_pix_) * 1000000.;
      */
	}

	/// computes residual of current solution - never called so far
	T score( )
	{ 
	  T norm(0);
	  for( int j = 0; j < N; j++ ) {
           norm += solution[j] * V[j];
        }
        printf( " %6.2f", b - norm );
//    if (!ppsOn )
  		return b - norm ;
//    else
//  		return (b - norm)/ sqrt(pps) ;
//  		return (b - norm)/ (pps) ; // fails
	}

	T score( const std::vector<T>& sol ) const
	{ 
    // <M*x-b> = x' M'M x + b'b - 2 b'Mx
    // solution = M'M^-1 M'b == x -> 
    // x'M'b + b'b - 2b'Mx = b'b - b'Mx -> bb-norm
    // in general though: 
    // compute M*x
	  T norm(0);

  	std::vector<T> Mx(N,0);
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
        Mx[i] += M[i*N+j] * sol[j];

	  for( int j = 0; j < N; j++ )
       norm += sol[j]*Mx[j] - 2.0*sol[j] * V[j];
   //printf( " %6.2f", b - norm );
//    if (!ppsOn )
    return b + norm ;
//    else
//  		return (b - norm)/ sqrt(pps) ;
//  		return (b - norm)/ (pps) ; // fails
	}

	private:

  // normal case:
#ifndef infogain
  void redo_ownscore( ) 
  {
#ifdef _newVariant_
    solve();
	  T norm(0);
	  for( int j = 0; j < N; j++ )
      norm += solution[j] * V[j];
//    printf( " %6.2f", b - norm );
    ownScore =  b - norm;
    std::copy(solution.begin(), solution.end(), own_solution.begin() );
#endif
  };
#else
	void redo_ownscore( )
	{ 
    solve();
	  T norm(0);
	  for( int j = 0; j < N; j++ )
      norm += solution[j] * V[j];
//    printf( " %6.2f", b - norm );
    ownScore =  b - norm;
    std::copy(solution.begin(), solution.end(), own_solution.begin() );
	}
#endif

	/// solves the current situation given by M^{-1} V and puts solution
	void solve_p(  ) // should be specialized ! float, double , etc.
	{
		ptrdiff_t n(N), m(N), nrhs(1), lda(N), ldb(N), lwork (-1), info;
		T wkopt;
		T* work;
		T* a = (T*) (&(M_h[0]));
		T* b = (T*) (&(V_h[0]));
		dgels( "No transpose", &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info );
		lwork = (ptrdiff_t) wkopt;
		work = (T*) malloc( lwork*sizeof(T) );

		// overrides b and xtra elements -- ???? and the matrix ! so i need to copy it
		dgels( "No transpose", &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork,&info );
	
		free ( work );
		std::copy(b, b+N, solution.begin() );
	}

	
	std::vector<T> M;
	std::vector<T> V;
	T b;
	std::vector<T> solution;
	std::vector<T> own_solution;	

	// local copies with more space for solver
	std::vector<T> M_h;
	std::vector<T> V_h;

  bool ppsOn; 
  T pps;
  /// own surrent score to be update : + , constructor
  T ownScore;
  T avScore;
};

typedef Matmatvecl2 <double,6> Matmatvecl2_d6;
typedef Matmatvecl2 <double,9> Matmatvecl2_d9;
typedef Matmatvecl2 <double,3> Matmatvecl2_d3;
