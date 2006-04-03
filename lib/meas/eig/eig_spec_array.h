/* ! $Id: eig_spec_array.h,v 3.0 2006-04-03 04:58:57 edwards Exp $ */

#ifndef __eig_spec_bj_array_w_h__
#define __eig_spec_bj_array_w_h__

#include "chromabase.h"
#include "linearop.h"

namespace Chroma {

void EigSpecRitzCG(const LinearOperatorArray<LatticeFermion>& H, // Herm pos def operator
		   multi1d<Real>& lambda_H,          // E-values
		   multi2d<LatticeFermion>& psi,     // E-vectors
		   int n_eig,                       // no of eig wanted
		   int n_renorm,                    // renorm frequency
		   int n_min,                       // minimum iters / e_value
		   int MaxCG,                       // Max no of CG iters
		   const Real& Rsd_r,               // relative residuum of each 
		                                  // e-value
		   const Real& Rsd_a,               // absolute target residuum
		                                    // (for small ev-s)
		   const Real& zero_cutoff,         // if evalue slips below this 
                                                    // we consider it as zero
		   const bool ProjApsiP,            // Project in Ritz?
		   
		   int& n_cg_tot,                   // Total no of CG iters
		   XMLWriter& xml_out         // Diagnostics
	      );

void EigSpecRitzKS(const LinearOperatorArray<LatticeFermion>& H, // Herm pos def operator
		   multi1d<Real>& lambda_H,          // E-values
		   multi2d<LatticeFermion>& psi,     // E-vectors
		   int n_eig,                       // no of eig wanted
		   int n_dummy,                     // No of Dummy e-vals to use

		   int n_renorm,                    // renorm frequency
		   int n_min,                       // minimum iters / e_value      
		   int n_max,                       // max iters / e_value
		   int n_max_KS,                    // max KS cycles
		   const Real& gamma_factor,        // the KS gamma factor

		   int MaxCG,                       // Max no of CG iters
		   const Real& Rsd_r,               // relative residuum of each 
		                                  // e-value
		   const Real& Rsd_a,               // absolute residuum (for small ev's)
		   const Real& zero_cutoff,         // if an ev slips below this
		                                    // we consider it zero
		   const bool ProjApsiP,            // Project in Ritz?
		   
		   int& n_cg_tot,                   // Total no of CG iters
		   int& n_KS,                       // No of KS cycles
		   int& n_jacob_tot,                // No of Jacobi Diag
		   XMLWriter& xml_out         // Diagnostics
	      );


void fixMMev2Mev(const LinearOperatorArray<LatticeFermion>& M,  // The Op to fix to
		 multi1d<Real>& lambda,                    // The Evals of M^{dag}M on input
		                                           // The Evals of M on output        
		 multi2d<LatticeFermion>& ev_psi,          // The Evecs corresponding to lambda
		 const int n_eig,                          // The no of evals/evecs to deal with
		 const Real& Rsd_r,                       // Relative error
		 const Real& Rsd_a,                       // Absolute error
		 const Real& zero_cutoff,                  // if EV slips below this we consider
		                                           // it to be zero
		 multi1d<bool>& valid_eig,                 // Validity mask (Write)
		 int& n_valid,                             // No of valids  (Write)
		 int& n_jacob                              // How many Jacobis were done
		 );

}  // end namespace Chroma

#endif
