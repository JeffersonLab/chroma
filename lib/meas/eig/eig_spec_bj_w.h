/* ! $Id: eig_spec_bj_w.h,v 1.4 2004-01-19 17:58:26 bjoo Exp $ */

#ifndef __eig_spec_bj_w_h__
#define __eig_spec_bj_w_h__

#include "chromabase.h"
#include "linearop.h"

void EigSpecRitzCG(const LinearOperator<LatticeFermion>& H, // Herm pos def operator
		   multi1d<Real>& lambda_H,          // E-values
		   multi1d<LatticeFermion>& psi,     // E-vectors
		   int n_eig,                       // no of eig wanted
		   int n_renorm,                    // renorm frequency
		   int n_min,                       // minimum iters / e_value
		   int MaxCG,                       // Max no of CG iters
		   const Real& Rsd_r,               // relative residuum of each 
		                                  // e-value
		   const bool ProjApsiP,            // Project in Ritz?
		   
		   int& n_cg_tot,                   // Total no of CG iters
		   XMLBufferWriter& xml_out         // Diagnostics
	      );

void EigSpecRitzKS(const LinearOperator<LatticeFermion>& H, // Herm pos def operator
		   multi1d<Real>& lambda_H,          // E-values
		   multi1d<LatticeFermion>& psi,     // E-vectors
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
		   const Real& Rsd_jacobi,          // Separate tolerance for Jacobi
		   const bool ProjApsiP,            // Project in Ritz?
		   
		   int& n_cg_tot,                   // Total no of CG iters
		   int& n_KS,                       // No of KS cycles
		   int& n_jacob_tot,                // No of Jacobi Diag
		   XMLBufferWriter& xml_out         // Diagnostics
	      );

#endif
