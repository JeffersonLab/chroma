// $Id: eig_spec_bj_w.cc,v 1.1 2004-01-16 15:38:37 bjoo Exp $
/*! \file
 *  \brief Compute low lying eigenvalues of the hermitian H
 */

#error "CONVERSION NOT COMPLETE: NEED TO MAKE APPROPRIATE SUBCLASSING"

#include "chromabase.h"

#include "meas/eig/eig_spec_bj_w.h"
#include "meas/eig/ritz.h"
#include "meas/eig/sn_jacobi.h"


using namespace QDP;

//! Compute low lying eigenvalues of the hermitian H 
/*!
 * \ingroup eig
 *
 *  Compute low lying eigenvalues of the hermitian H 
 *  using the Ritz functional minimization routine,
 *  if desired in the Kalkreuter-Simma algorithm
 *  H		The operator        		(Read)
 *  lambda_H    Eigenvalues                     (Modify)
 *  psi		Eigenvectors			(Modify)
 *  N_eig	No of eigenvalues		(Read)
 *  N_min	Minimal CG iterations		(Read)
 *  N_max	Maximal CG iterations		(Read)
 *  Kalk_Sim	Use Kalkreuter-Simma criterion	(Read)
 *  n_renorm	Renormalize every n_renorm iter.	(Read)
 *  N_KS_max	Max number of Kalkreuter-Simma iterations	(Read)
 *  RsdR_r	(relative) residue		(Read)
 *  Cv_fact	"Convergence factor" required	(Read)
 *  NCG_tot	Total number of CG iter		(Write)
 *  n_KS	total number of Kalkreuter-Simma iterations	(Write)
 *  n_valid	number of valid eigenvalues	(Write)
 *  ProjApsiP	flag for projecting A.psi	(Read)
 */


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
	      )
{
  END_CODE("Eig_spec");
}
