/* ! $Id: eig_spec_bj_w.h,v 1.1 2004-01-16 15:38:37 bjoo Exp $ */

#ifndef __eig_spec_bj_w_h__
#define __eig_spec_bj_w_h__

#include "chromabase.h"


void Eig_Spec(LinearOperator<LatticeFermion>& H, // Herm pos def operator
	      multi1d<Real> lambda_H,            // E-values
	      multi1d<LatticeFermion> psi,       // E-vectors
	      int n_eig,                         // no of eig wanted
	      int n_dummy,                       // no of dummy values
	      int n_renorm,                      // renorm frequency
	      int n_min,                         // minimum iters / e_value
	      int MaxCG,                         // Max no of CG iters
	      const Real& Rsd_r,                 // relative residuum of each 
                                                 // e-value
	      const bool ProjAPsiP,              // Project in Ritz?

	      const bool Kalk_Sim,               // KS mode
	      int n_max,                         // max iters / KS cycle
	      int n_KS_max,                      // max KS cycles
	      const Real& gamma_factor,          // Convergence factor gamma
	     
	      int& n_cg_tot,                     // Total no of CG iters
	      int& n_KS_tot,                     // Total no of KS cycles
	      multi1d<bool> n_valid,             // which ev-s are valid  
	      XMLBufferWriter& xml_out           // Diagnostics
	      );


#endif
