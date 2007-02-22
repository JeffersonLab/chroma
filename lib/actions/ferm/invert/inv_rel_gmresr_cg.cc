#include "chromabase.h"
#include "actions/ferm/invert/inv_rel_cg1.h"
#include "actions/ferm/invert/inv_rel_gmresr_cg.h"

namespace Chroma {

    
template<typename T>
void InvRelGMRESR_CG_a(const LinearOperator<T>& PrecMM,
			 const LinearOperator<T>& UnprecMM,
			 const T& b,
			 T& x,
			 const Real& epsilon, 
			 const Real& epsilon_prec,
			 int MaxGMRESR, 
			 int MaxGMRESRPrec,
			 int& n_count)
{
  const Subset&  s= UnprecMM.subset();

  x[s] = zero;
  T r; 
  r[s] = b;

  // Arrays of pointers for C and U, up to the maximum number
  multi1d<T*> C(MaxGMRESR);
  multi1d<T*> U(MaxGMRESR);
  int C_size=0;

  // Get || r ||
  Double norm_r = sqrt(norm2(r,s));
  Double terminate = sqrt(norm2(b,s))*epsilon;

  int iter=0;

  int prec_count;
  int prec_count_total = 0;

  // Do the loop
  while( toBool( norm_r > terminate) && iter < MaxGMRESR ) {
    // First we have to get a new u vector (in U)
    T* u = new T;
    (*u)[s] = zero;

    
    // Now solve with the preconditioned system to preconditioned accuracy
    // THis is done with CG
    // solve      MM u = r
    InvRelCG1(PrecMM, r, *u, epsilon_prec, MaxGMRESRPrec, prec_count);

    // Get C
    T* c = new T;
    (*c)[s] = zero;

    // Now Apply Matrix A, to u, to precision epsilon || b ||/ || r ||
    Real inner_tol = terminate / norm_r;

    // Apply zeta + rho U
    // First apply U
    UnprecMM( *c, *u, PLUS, inner_tol);


    // Now there is some orthogonalisation to do
    for(int i =0; i < C_size; i++) { 
      Complex beta = innerProduct( *(C[i]), *(c), s );
      (*c)[s] -= beta*(*(C[i]));
      (*u)[s] -= beta*(*(U[i]));
    }

    // Normalise c
    Double norm_c = sqrt(norm2( *c, s));

    (*c)[s] /= norm_c;
    (*u)[s] /= norm_c;

    // Hook vectors into the U, C arrays
    C[C_size] = c;
    U[C_size] = u;
    C_size++;

    Complex alpha = innerProduct( *c, r, s);
    
    x[s] += alpha*(*u);
    r[s] -= alpha*(*c);

    iter++;
    norm_r = sqrt(norm2(r,s));
    QDPIO::cout << "Inv Rel GMRESR: iter "<< iter <<" || r || = " << norm_r << endl;
  }

  // Cleanup
  for(int i=0; i < C_size; i++) { 
    delete C[i];
    delete U[i];
  }

  if( iter == MaxGMRESR ) { 
    QDPIO::cout << "Nonconvergence warning " << endl;
  }
  
  n_count = iter;
}

template<>
void InvRelGMRESR_CG(const LinearOperator<LatticeFermion>& PrecMM,
		     const LinearOperator<LatticeFermion>& UnprecMM,
		     const LatticeFermion& b,
		     LatticeFermion& x,
		     const Real& epsilon, 
		     const Real& epsilon_prec,
		     int MaxGMRESR, 
		     int MaxGMRESRPrec,
		     int& n_count)
{
  InvRelGMRESR_CG_a(PrecMM, UnprecMM, b, x, epsilon, epsilon_prec, MaxGMRESR, MaxGMRESRPrec, n_count);
}

}  // end namespace Chroma
