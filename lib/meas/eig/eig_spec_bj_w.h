/* ! $Id: eig_spec_bj_w.h,v 1.3 2004-01-16 16:51:19 bjoo Exp $ */

#ifndef __eig_spec_bj_w_h__
#define __eig_spec_bj_w_h__

#include "chromabase.h"
#include "linearop.h"

void EigSpecRitzCG(const LinearOperator<LatticeFermion>& H, // Herm pos def operator
		   multi1d<Real> lambda_H,          // E-values
		   multi1d<LatticeFermion> psi,     // E-vectors
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


#endif
