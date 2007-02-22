// $Id: eig_spec_array.cc,v 3.1 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Compute low lying eigenvalues of the hermitian H
 */


#include "chromabase.h"

#include <sstream>

#include "meas/eig/eig_spec_array.h"
#include "meas/eig/ritz_array.h"
#include "meas/eig/sn_jacob_array.h"

namespace Chroma {

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


void EigSpecRitzCG(const LinearOperatorArray<LatticeFermion>& M, // Herm pos def operator
		   multi1d<Real>& lambda_H,          // E-values
		   multi2d<LatticeFermion>& psi,     // E-vectors
		   int n_eig,                        // No of e-values to find
		   int n_renorm,                     // renorm frequency
		   int n_min,                        // minimum iters / e_value
		   int MaxCG,                        // Max no of CG iters
		   const Real& Rsd_r,                // relative residuum of each 
		                                     // e-value
		   const Real& Rsd_a,                // absolute residuum
		   const Real& zero_cutoff,          // if ev slips below this
		                                     // we consider it zero
		   const bool ProjApsiP,             // Project in Ritz?
		   
		   int& n_cg_tot,                    // Total no of CG iters
		   XMLWriter& xml_out          // Diagnostics
	     )
{
  QDPIO::cout << "EigSpecArray" << endl;

  START_CODE();
  
  const Subset& sub = M.subset(); // Subset over which M acts
  
  push(xml_out, "EigSpecRitzCG");

  n_cg_tot = 0;
  int n_count = 0;
  int n_max = MaxCG+1;

  int n, i;
  Real final_grad;

  multi1d<Real> resid_rel(n_eig);
  for(n = 1, i = 0; n <= n_eig;  n++, i++) 
  {
    // Initialise lambda[i] = 1
    lambda_H[i] = Real(1);
    
    QDPIO::cout << "Call ritz n=" << n << endl;

    // Do the Ritz
    Ritz(M, lambda_H[i], psi, n, Rsd_r, Rsd_a, zero_cutoff, n_renorm, n_min, 
	 n_max, MaxCG, ProjApsiP, n_count, final_grad, false, Real(1), 
	 Real(1));

    // Add n_count
    n_cg_tot += n_count;
    
    // Check e-value
    ostringstream s;
    s << "eval" << n;
    push(xml_out, s.str());
    write(xml_out, "n_count", n_count);
    pop(xml_out);

    if( toBool( fabs(lambda_H[i]) < zero_cutoff ) ) { 
      QDPIO::cout << "Evalue["<< n << "] = " << lambda_H[i] << " is considered zero" << endl;
    }
    else 
    {
      multi1d<LatticeFermion> D_e(M.size());
      multi1d<LatticeFermion> lambda_e(M.size());

      M(D_e, psi[i], PLUS);

      for(int j=0; j < M.size(); ++j)
      {
	lambda_e[j][sub] = lambda_H[i]*psi[i][j];
	D_e[j][sub] -= lambda_e[j];
      }
      Double r_norm = sqrt(norm2(D_e,sub));
      resid_rel[i] = Real(r_norm)/lambda_H[i];
      QDPIO::cout << "Evalue["<<n<<"]: eigen_norm = " << r_norm << " resid_rel = " << resid_rel[i] << endl << endl;
    }
  }

  // All EV-s done. Dump-em
  push(xml_out, "eValues");
  write(xml_out, "lambda", lambda_H);
  pop(xml_out);
  write(xml_out, "n_cg_tot", n_cg_tot);

  pop(xml_out); // EigSpecRitz

      
  END_CODE();
}


void EigSpecRitzKS(const LinearOperatorArray<LatticeFermion>& M, // Herm pos def operator
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
		   const Real& Rsd_a,               // Absolute residuum
		   const Real& zero_cutoff,         // if e-value slips below this, we 
		                                    // consider it zero
		   const bool ProjApsiP,            // Project in Ritz?
		   
		   int& n_cg_tot,                   // Total no of CG iters
		   int& n_KS,                       // Total no of KS cycles
		   int& n_jacob_tot,
		   XMLWriter& xml_out         // Diagnostics
	      )
{
  START_CODE();

  const Subset& s = M.subset(); // Subset over which M acts
  
  // Sanity Checks: 
  // Make sure lambda_H is large enough
  if( lambda_H.size() < (n_eig+n_dummy) ) { 
    QDP_error_exit("lambda_H is too small to hold n_eig + n_dummy values\n");
  }

  // Make sure psi is large enough
  if( psi.size2() < (n_eig + n_dummy) ) {
    QDP_error_exit("psi is too small to hold n_eig + n_dummy values\n");
  }

  if( n_eig <=0 ) { 
    QDP_error_exit("n_eig must be > 0. it is %d\n", n_eig);
  }

  int N5 = psi.size1();

  if( n_eig ==1 ) { 
    // if n_eig is one, KS algorithm is not applicable. We revert to the 
    // Normal CG method
    EigSpecRitzCG(M, lambda_H, psi, n_eig, n_renorm, n_min, MaxCG, Rsd_r, 
		  Rsd_a, zero_cutoff, ProjApsiP, n_cg_tot, xml_out);
    return;
  }

  // Internal lambda's
  // lambda_H holds initial guesses. Copy these over.
  int n_working_eig = n_eig+n_dummy;
  multi1d<Real> lambda_intern(n_working_eig);
  multi1d<Real> lambda_previous(n_working_eig);
  multi1d<Real> delta_cycle(n_working_eig);
  int i;
  for(i = 0; i < n_working_eig; i++) { 
    lambda_intern[i] = lambda_H[i];
    lambda_previous[i] = Real(1);
    delta_cycle[i] = Real(1);
  }

  // Off diag elements of hermitian matrix to diagonalise
  multi1d<Complex> off_diag((n_working_eig)*(n_working_eig-1)/2);

  
  multi1d<Real> final_grad(n_working_eig);

  int n_count = 0;
  n_cg_tot = 0;
  n_KS = 0;
  n_jacob_tot = 0;

  int ev, j, ij;
  int n_jacob;

  push(xml_out, "EigSpecRitzKS");
  
  bool convP = false;

  for( int KS_iter = 0; KS_iter < n_max_KS; KS_iter++) 
  {
    n_KS++;

    push(xml_out, "KS_iter");

    // KS Step 1: for each k=0, n-1 in succession compute
    // approximations to the eigenvectors of M by performing
    // only a certain number of CG iterations (governed by gamma_factor)
    for(ev = 1, i=0; ev <= n_working_eig; ev++, i++) {

      // Do the Ritz for all working eigenvalues
      Ritz(M, lambda_intern[i], psi, ev, Rsd_r, Rsd_a, zero_cutoff,  n_renorm, 
	   n_min, n_max, MaxCG, ProjApsiP, n_count, final_grad[i], true, 
	   Real(1), gamma_factor);

      // Count the CG cycles
      n_cg_tot += n_count;
      push(xml_out, "ev");
      write(xml_out, "n_count", n_count);
      write(xml_out, "final_grad", final_grad[i]);
      pop(xml_out);

    }
    
    // Construct the hermitian matrix M to diagonalise.
    // Note only the off diagonal elements are needed
    // the lambda_intern[i] are the diagonal elements

    // M_ij = (psi_j, M psi_i )
    //
    multi1d<LatticeFermion> tmp(N5);
    for(i=0, ij=0; i < n_working_eig; i++) { 
      M(tmp, psi[i], PLUS);
      for(j=0; j < i; j++) { 
	off_diag[ij] = innerProduct(psi[j][0], tmp[0], s);
	for(int n = 1; n < N5; n++) { 
	  off_diag[ij] += innerProduct(psi[j][n], tmp[n], s);
	}
	ij++;
      }
    }



    // Now diagonalise it, rotate evecs, and sort
    // 
    // Jacobi at the moment works with the absolute error
    n_jacob = SN_Jacob_Array(psi, n_working_eig, lambda_intern, off_diag, Rsd_a, 50, s);
    n_jacob_tot += n_jacob;

    write(xml_out, "n_jacob", n_jacob);

    
    // Now we should check convergence
    // Pessimistic but safe criterion || g ||/lambda  < omega
    // but can also be converged if we consider something to have a zero ev
    // we only check convergence of the wanted n_eig
    bool convP = true;
    for(i=0; i < n_eig; i++) {
      bool zeroReachedP = toBool( fabs(lambda_intern[i]) < zero_cutoff );

      // This ev is converged if either 
      //   i)  Pessimistic absolute error is reached 
      //  ii)  Pessimistic relative error is reached 
      // iii)  if the ev is below the zero cutoff
      // 
      convP &= ( toBool( final_grad[i] < Rsd_a ) 
		 || toBool( final_grad[i] < Rsd_r*fabs(lambda_intern[i]) ) 
		 || zeroReachedP );
    }

    pop(xml_out); // KS_iter

    // All the wanted e-values converged
    if( convP ) {

      // Do final Jacobi without the dummies.
      // Make the matrix
      for(i=0, ij=0; i < n_eig; i++) { 
	M(tmp, psi[i], PLUS);
	for(j=0; j < i; j++) { 
	  off_diag[ij] = innerProduct(psi[j][0], tmp[0], s);
	  for(int n=1; n < N5; n++) { 
	    off_diag[ij] += innerProduct(psi[j][n], tmp[n], s);
	  }
	  ij++;
	}
      }

      // Diagonalise, rotate, sort
      n_jacob = SN_Jacob_Array(psi, n_eig, lambda_intern, off_diag, Rsd_a, 50, s);
      write(xml_out, "final_n_jacob", n_jacob);
      write(xml_out, "n_cg_tot", n_cg_tot);
      write(xml_out, "n_KS", n_KS);

      
      // Copy lambda_intern back into lambda_H and return
      for(i=0; i < n_eig; i++) { 
	lambda_H[i] = lambda_intern[i];
      }

      write(xml_out, "lambda_H", lambda_H);
      pop(xml_out); // EigSpecRitzKS
      END_CODE();
      return;
    }

  }

  write(xml_out, "n_cg_tot", n_cg_tot);
  write(xml_out, "n_KS", n_KS);
  pop(xml_out); //EigSpecRitzKS
  
  // If we reached here then we have done more than n_max KS
  QDP_error_exit("n_max_KS reached with no convergence");
  END_CODE();
}


void fixMMev2Mev(const LinearOperatorArray<LatticeFermion>& M,  // The Op to fix to
		 multi1d<Real>& lambda,       // The Evals of M^{dag}M on input
		                             // The Evals of M on output 
		 multi2d<LatticeFermion>& ev_psi,  // The Evecs 
		 const int n_eig,             // The no of evals/evecs 
		 const Real& Rsd_r,           // Relative error for validity
		 const Real& Rsd_a,           // Absolute error for validity
		 const Real& zero_cutoff,     // Zero cutoff
		 multi1d<bool>& valid_eig,    // Validity mask (Write)
		 int& n_valid,                // No of valids  (Write)
		 int& n_jacob                 // How many Jacobis were done
		 )
{
  START_CODE();

  const Subset& s = M.subset(); // Subset over which M acts
  

  // Sanity checking
  if( n_eig > lambda.size() ) { 
    QDP_error_exit("n_eig greater than size of lambda array\n");
  }

  if( lambda.size() != ev_psi.size2() ) { 
    QDP_error_exit("lambda and ev_psi arrays must have same size\n");
  }

  // Make the valid eig the right size
  valid_eig.resize(n_eig);

  int N5 = ev_psi.size1();

  // Temporaries 
  multi1d<LatticeFermion> tmp(N5);
  Double lambda_H_sq;
  Double delta_lambda;
  bool zeroMatchedP;
  bool convP;
  // We are all set -- lets do it

  if( n_eig == 1 ) {                  // only one eigenvalue
    M(tmp, ev_psi[0], PLUS);
    Double lambda_fix_single = innerProductReal(ev_psi[0][0], tmp[0], s);
    for(int n=1; n < N5; n++) { 
      lambda_fix_single += innerProductReal(ev_psi[0][n], tmp[n], s);
    }
    // No diagonalisation needed -- only 1 eval

    // We square lambda_fix_single.
    lambda_H_sq = lambda_fix_single*lambda_fix_single;


    // Different validity criteria, depending on whether lambda_fix_single^2 
    // is less than our zero cutoff
    if( toBool( lambda_H_sq < zero_cutoff ) ) 
    {
      // Yes, the square of our estimate for lambda is less than the cutoff
      // Check whether the original lambda was also less than the cutoff
      // if so, our ev is considered a valid zero e-value
      if( toBool( lambda[0] < zero_cutoff ) ) {
	valid_eig[0] = true;
	n_valid = 1;
      }
      else { 
	// No, the square of our estimate was not zero. Mismatch
	// BTW what else can we do here?
	valid_eig[0] = false;
	n_valid = 1;
      }
    }
    else { 
      // Check relative error or absolute error is satisfied
      delta_lambda = fabs( lambda_H_sq - Double( lambda[0] ) );
      
      convP = toBool( delta_lambda < Double(Rsd_a) )
	|| toBool( delta_lambda < Double(Rsd_r)*fabs(Double(lambda[0])) );
      
      if ( convP ) { 
	// Condition is satisfied
	valid_eig[0] = true;
	n_valid = 1;
      } else {
	// Condition is not satisfied
	valid_eig[0] = false;
	n_valid = 0;
      }
    }
 
    // Whatever happened lambda_fix_single is our estimate for
    // lambda[0]
    lambda[0] = Real(lambda_fix_single);
    return;
  }  // end if (n_eig == 1)

    
  // Interesting case multiple ev-s
  // Construct new ev-s and off diagonal matrix
  multi1d<Complex> off_diag(n_eig*(n_eig-1)/2 );
  multi1d<Real>    lambda_fix(n_eig);
  for(int i = 0; i < n_eig; i++) { 
    valid_eig[i] = false;
  }

  for(int i=0, ij=0; i < n_eig; i++) { 
    M(tmp, ev_psi[i], PLUS); 
    lambda_fix[i] = innerProductReal(ev_psi[i][0], tmp[0], s);
    for(int n=1; n < N5; n++) { 
      lambda_fix[i] += innerProductReal(ev_psi[i][n], tmp[n], s);
    }

    for(int j=0; j < i; j++) { 
      off_diag[ij] = innerProduct(ev_psi[j][0], tmp[0], s);
      for(int n=1; n < N5; n++) { 
	off_diag[ij] += innerProduct(ev_psi[j][n], tmp[n], s);
      }
      ij++;
    }
  }
    
  // Diagonalise and sort according to size
  n_jacob = SN_Jacob_Array(ev_psi, n_eig, lambda_fix, off_diag, Rsd_a, 50, s);

  // Now the tricky business of matching up with ev-s of H^2
  n_valid = 0;
    
  for(int i=0; i < n_eig; i++) 
  { 
    lambda_H_sq = Double(lambda_fix[i])*Double(lambda_fix[i]);

    if( toBool(lambda_H_sq < zero_cutoff ) ) 
    {
      // Now compare with all the invalid ev-s
      for(int j = n_valid; ( j < n_eig ) && ( valid_eig[i] == false ); j++ ) 
      {
	// If we find an original zero lambda, we match it with this 
	// and consider this matched
	if( toBool( lambda[j] < zero_cutoff ) ) { 
	  
	  // lambda_fix[i] matches lambda[j] 
	  valid_eig[i] = true;
	  
	  // swap lambda[j] with lambda[valid_eig];
	  // so we don't search it again
	  Real tmp = lambda[j];
	  lambda[j] = lambda[n_valid];
	  lambda[n_valid] = tmp;
	  
	  // Increase the validity count
	  n_valid++;
	}

	// Else compare with next j
      }
      // At this point either valid_eig[i] is true, or we were unable to match
      // in which case valid_eig[i] is false from initialisation
      // in any case time for next_eig

    }
    else 
    {
      // Now compare with all the invalid ev-s
      for(int j = n_valid; ( j < n_eig ) && ( valid_eig[i] == false ); j++ ) 
      {
	delta_lambda = fabs( lambda_H_sq - Double(lambda[j]) );
	
	convP = toBool( delta_lambda < Double(Rsd_a) ) 
	  || toBool( delta_lambda < Double(Rsd_r) * fabs( Double(lambda[j]) ) );
	if( convP ) 
	{
	  // lambda_fix[i] matches lambda[j] 
	  valid_eig[i] = true;
	  
	  // swap lambda[j] with lambda[valid_eig];
	  // so we don't search it again
	  Real ftmp = lambda[j];
	  lambda[j] = lambda[n_valid];
	  lambda[n_valid] = ftmp;

	  // Increase the validity count
	  n_valid++;
	}

	// Else compare with next j
      }
      // At this point either valid_eig[i] is true, or we were unable to match
      // in which case valid_eig[i] is false from initialisation
      // in any case time for next_eig
    }
  } // We are done

  // Now we can overwrite lambda with lambda_fix
  // A desirable feature would be to move the invalid ev-s to the end
  // but I cannot be bothered
  for(int i = 0; i < n_eig; i++ ) { 
    lambda[i] = lambda_fix[i];
  }
  
  END_CODE();

  // we are done 
  return;
}

}  // end namespace Chroma
