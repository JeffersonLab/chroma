#include "inv_multiprec_richardson.h"

namespace Chroma { 


  /*! This is the multi precision inverter
    Solves systems via iterative Refinement.
    The interface specifies explicit single/double precision references
    These should be set up by the caller. 
    Currently I have made the SyssolverLinOpRichardson and the SyssolverMdagMRichardson 
    do this 
  */
  void InvMultiPrecRichardson( const SystemSolver< LatticeFermionF >& Dinv,
			       const LinearOperator< LatticeFermionD >& D,
			       const LatticeFermionD& b, 
			       LatticeFermionD& x,
			       int MaxIter,
			       Real RsdTarget,
			       SystemSolverResults_t& res)
  {
    START_CODE();

    
    LatticeFermionD r;

    const Subset& s = D.subset();

    int iter = 0;

    // Target Residue
    Double rsd_t = Double(RsdTarget)*Double(RsdTarget)*norm2(b,s);

    // Compute Initial residue:
    {
      LatticeFermionD tmp;
      D(tmp, x, PLUS);
      r[s] = b;
      r[s] -= tmp;      // r = b-Ax
    }

    Double rnorm = norm2(r, s);   // || r ||^2
    QDPIO::cout << "Initial Norm: " << rnorm << endl;

    if( toBool( rnorm <= rsd_t ) ) { 

      // We are done
      QDPIO::cout << "Iter 0: || r || = "<< rnorm << "  || r_target || = " << rsd_t << endl;

      res.n_count = 0;
      res.resid = Real(sqrt(rnorm));

      return;
    }
    
    for(int i=1; i <= MaxIter ; i++) { 

      LatticeFermionD delta_x;         // The Richardson Correction
      LatticeFermionD Ddelta_x;        // D delta_x

      // Compute Delta_x = D^{-1} r
      LatticeFermionF r_single(r);     // Downcast r to single
      LatticeFermionF dx_single;   // Initial guess for dx solve
      

      
      dx_single[s] = zero;

      Dinv(dx_single, r);              // Solve in single

      delta_x[s] = dx_single;          // Upcast back to double
      
      D(Ddelta_x, delta_x, PLUS);      // Get D delta_x


      x[s] += delta_x;
      r[s] -= Ddelta_x;

      rnorm = norm2(r, s); 

      // Convergence check
      if( toBool( rnorm <= rsd_t ) ) { 
	// Compute Initial residue:
	{
	  LatticeFermionD tmp;
	  D(tmp, x, PLUS);
	  r[s] = b;
	  r[s] -= tmp;      // r = b-Ax
	}
	// We are done
	
	res.n_count = i;
	res.resid = Real(sqrt(norm2(r)/norm2(b)));
	return;
      }

    }

    QDPIO::cout << "Richardson Multi Prec Solver: NONCONVERGENCE" << endl;
    res.n_count = MaxIter;
    res.resid = Real(sqrt(rnorm));

    END_CODE();
  }
			       



}
