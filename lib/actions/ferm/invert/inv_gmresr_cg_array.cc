#include "chromabase.h"
#include "actions/ferm/invert/invcg1_array.h"
#include "actions/ferm/invert/inv_gmresr_cg_array.h"


namespace Chroma 
{
    
template<typename T>
void InvGMRESR_CG_a(const LinearOperatorArray<T>& PrecMM,
		    const LinearOperatorArray<T>& PrecM,
		    const LinearOperatorArray<T>& UnprecM,
		    const multi1d<T>& b,
		    multi1d<T>& x,
		    const Real& epsilon, 
		    const Real& epsilon_prec,
		    int MaxGMRESR, 
		    int MaxGMRESRPrec,
		    int& n_count)
{
  START_CODE();

  const Subset&  s= UnprecM.subset();
  int N5 = UnprecM.size();

  // Sanity checks
  if( b.size() != N5 ) { 
    QDPIO::cerr << "b has wrong size in 5th dimension. N5 = " << N5 << " b.size() = " << b.size() << endl;
  }

  if( x.size() != N5 ) { 
    QDPIO::cerr << "x has wrong size in 5th dimension. N5 = " << N5 << " x.size() = " << x.size() << endl;
  }


  multi1d<T> r(N5);
  // x_0 = 0
  // r_0 = b
  for(int n5=0; n5 < N5; n5++ ) {
    x[n5][s] = zero;
    r[n5][s] = b[n5];
  }

  // Arrays of pointers for C and U, up to the maximum number
  multi1d<multi1d<T>*> C(MaxGMRESR);
  multi1d<multi1d<T>*> U(MaxGMRESR);
  int C_size=0;

  // Get || r ||
  Double norm_r = sqrt(norm2(r,s));
  Double terminate = sqrt(norm2(b,s))*epsilon;
  int iter=0;

  int prec_count;
  // Do the loop
  while( toBool( norm_r > terminate) && iter < MaxGMRESR ) {
    // First we have to get a new u vector (in U)
    multi1d<T>* u = new multi1d<T>;
    (*u).resize(N5);
     
    /*
    for(int n5=0; n5 < N5; n5++) {
      (*u)[n5][s] = zero;
    }
    */

    multi1d<T> tmp(N5);
    for(int n5=0;n5 < N5; n5++) { 
      (*u)[n5][s]= zero;
    }
    // Now solve with the preconditioned system to preconditioned accuracy
    // THis is done with CG
    // solve      MM u = r

    PrecM(tmp,r , MINUS);
    InvCG1(PrecMM, tmp, *u, epsilon_prec, MaxGMRESRPrec, prec_count);

  
    // Get C
    multi1d<T>* c = new multi1d<T>;
    (*c).resize(N5);

    // c = A u
    UnprecM( *c, *u, PLUS);


    // Now there is some orthogonalisation to do
    for(int i =0; i < C_size; i++) { 

      // Compute < C[i], c > in 5D
      Complex beta = Real(0);
      for(int n5=0; n5 < N5; n5++) {
	beta += innerProduct( (*(C[i]))[n5], (*c)[n5], s );
      }
      
      // Compute  c -= beta*C[i]
      //          u -= beta*U[i]
      for(int n5=0; n5 < N5; n5++) {
	(*c)[n5][s] -= beta*(*(C[i]))[n5];
	(*u)[n5][s] -= beta*(*(U[i]))[n5];
      }
    }

    
    // norm2 works in 5D
    Double norm_c = sqrt(norm2( *c, s));

    // Now normalise c and u
    for(int n5=0; n5 < N5; n5++) {
      (*c)[n5][s] /= norm_c;
      (*u)[n5][s] /= norm_c;
    }
    // Hook vectors into the U, C arrays
    C[C_size] = c;
    U[C_size] = u;
    C_size++;

    Complex alpha = Real(0);
    for(int n5=0; n5 < N5; n5++ ) {
      alpha += innerProduct( (*c)[n5], r[n5], s);
    }

    for(int n5=0; n5 < N5; n5++) { 
      x[n5][s] += alpha*(*u)[n5];
      r[n5][s] -= alpha*(*c)[n5];
    }
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
  END_CODE();
}

template<>
void InvGMRESR_CG(const LinearOperatorArray<LatticeFermion>& PrecMM,
		  const LinearOperatorArray<LatticeFermion>& PrecM,
		  const LinearOperatorArray<LatticeFermion>& UnprecM,
		  const multi1d<LatticeFermion>& b,
		  multi1d<LatticeFermion>& x,
		  const Real& epsilon, 
		  const Real& epsilon_prec,
		  int MaxGMRESR, 
		  int MaxGMRESRPrec,
		  int& n_count)
{
  InvGMRESR_CG_a(PrecMM, PrecM, UnprecM, b, x, epsilon, epsilon_prec, MaxGMRESR, MaxGMRESRPrec, n_count);
}

}  // end namespace Chroma
