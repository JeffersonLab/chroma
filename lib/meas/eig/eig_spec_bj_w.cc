// $Id: eig_spec_bj_w.cc,v 1.3 2004-01-16 16:51:19 bjoo Exp $
/*! \file
 *  \brief Compute low lying eigenvalues of the hermitian H
 */


#include "chromabase.h"

#include <sstream>

using namespace std;

#include "meas/eig/eig_spec_bj_w.h"
#include "meas/eig/ritz.h"
#include "meas/eig/sn_jacob.h"


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


void EigSpecRitzCG(const LinearOperator<LatticeFermion>& M, // Herm pos def operator
		 multi1d<Real> lambda_H,            // E-values
		 multi1d<LatticeFermion> psi,       // E-vectors
		 int n_eig,                        // No of e-values to find
		 int n_renorm,                      // renorm frequency
		 int n_min,                         // minimum iters / e_value
		 int MaxCG,                         // Max no of CG iters
		 const Real& Rsd_r,                 // relative residuum of each 
	     // e-value
		 const bool ProjApsiP,              // Project in Ritz?
	     
		 int& n_cg_tot,                     // Total no of CG iters
		 XMLBufferWriter& xml_out           // Diagnostics
	     )
{
  START_CODE("EigSpecRitzCG");
  
  
  push(xml_out, "EigSpecRitzCG");


  n_cg_tot = 0;
  int n_count = 0;
  int n_max = MaxCG+1;

  int n, i;

  multi1d<Real> resid_rel(n_eig);
  for(n = 1, i = 0; n <= n_eig;  n++, i++) {

    // Initialise lambda[i] = 1
    lambda_H[i] = Real(1);
    
    // Do the Ritz
    Ritz(M, lambda_H[i], psi, n, Rsd_r, n_renorm, n_min, n_max, 	   
	 MaxCG, ProjApsiP, n_count, false, Real(1), Real(1));

    // Add n_count
    n_cg_tot += n_count;
    
    // Check e-value
    LatticeFermion   D_e;
    LatticeFermion  lambda_e;
    M(D_e, psi[i], PLUS);
    lambda_e = lambda_H[i]*psi[i];
    D_e -= lambda_e;
    Double r_norm = sqrt(norm2(D_e));
    resid_rel[i] = Real(r_norm)/lambda_H[i];
    QDPIO::cout << "Evalue["<<n<<"]: eigen_norm = " << r_norm << " resid_rel = " << resid_rel[i] << endl << endl;
    // Dump info 
    ostringstream s;
    s << "eval" << n;
    push(xml_out, s.str());
    write(xml_out, "n_count", n_count);
    write(xml_out, "e_norm", r_norm);
    write(xml_out, "rel_e_norm", resid_rel[i]);
    pop(xml_out);
  }

  // All EV-s done. Dump-em
  push(xml_out, "eValues");
  write(xml_out, "lambda", lambda_H);
  pop(xml_out);
  write(xml_out, "n_cg_tot", n_cg_tot);

  pop(xml_out); // EigSpecRitz

      
  END_CODE("EigSpecRitzCG");
}
