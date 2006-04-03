/* $Id: unprec_ovlap_contfrac5d_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
*  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
*/

#include "chromabase.h"
#include "linearop.h"

#include "actions/ferm/linop/unprec_ovlap_contfrac5d_linop_array_w.h"


namespace Chroma 
{ 
  //! Apply the operator onto a source vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecOvlapContFrac5DLinOpArray::operator() (multi1d<LatticeFermion>& chi,
					       const multi1d<LatticeFermion>& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);
    int G5 = Ns*Ns - 1;

    // This is the upper limit for the index of the 5th dimension, i.e.
    //  it runs from 0 to TwoN. (altogether TwoN + 1 numbers. This
    // makes sense since N5 is ALWAYS odd. )
    int TwoN = N5 - 1;

    // This is our scaling: we evaluate 
    //  ( 1 + m_q ) / ( 1 - m_q ) gamma5 + eps
    //
    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);

    // The sign of the diagonal terms.
    // We flip Hsign BEFORE using it 
    // This makes the first element always have a PLUS sign
    int Hsign=-1;

    /* N5 is ALWAYS ODD hence the explicit setting above */
    /*
      if( N5%2 == 0 ) {
      Hsign = 1;
      }
      else {
      Hsign = -1;
      }
    */

    // Run through all the pseudofermion fields:
    //   chi(0) = beta(0)*H*psi(0) + alpha(0)*psi(1)
    //   chi(n) = alpha(n-1)*psi(n-1)
    //              + (-)^n*beta(n)*H*psi(n) + alpha(n)*psi(n+1)
    //   chi(TwoN) = alpha(TwoN-1)*psi(TwoN-1)
    //                      + (-)^TwoN*beta(TwoN)*H*psi(TwoN)         
    
    LatticeFermion tmp1, tmp2;
    Real pmscale;

    for(int n = 0; n < TwoN; ++n) {
      (*M)(tmp1, psi[n], PLUS);      // tmp1 = M psi[n]
      tmp2 = Gamma(G5)*tmp1;        // tmp2 = gamma_5 M psi

      // Scale factor and sign for the diagonal term proportional to H
      // The scale factor should be chosen in the fermact call such
      //  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
      Hsign = -Hsign;
      pmscale = beta[n]*Hsign*scale_fac;
      chi[n] = pmscale*tmp2;

      if( N5 > 1 ) { 
	chi[n] += alpha[n] * psi[n+1];
      }

      if( n > 0 ) {
	chi[n] += alpha[n-1]*psi[n-1];
      }
    }

    // Last Component
    // chi[N] = mass*gamma5*psi[N] + alpha[N-1]*psi[N-1]
    tmp1 = Gamma(G5)*psi[TwoN];
    chi[TwoN] = mass*tmp1;;
    chi[TwoN] += alpha[TwoN-1]*psi[ TwoN -1 ];

    // Project out eigenvectors from Source if desired 
    //
    // (What does this do then? Urs????)
    LatticeFermion psi_proj;
    psi_proj = zero;


    if(  NEig > 0 ) {
      Complex cconsts;
      for(int i = 0 ; i < NEig; i++) {
      
	// Project out eigenvectors in psi_{2N} */
	cconsts = innerProduct(EigVec[i],psi[TwoN]);
	psi_proj -= EigVec[i]*cconsts;

	// The vectors are added into chi_2N
	chi[TwoN] += cconsts*EigValFunc[i]*EigVec[i];

	// Project out the eigenvectors from psi_{2N-1} and subtract from chi[N]
	cconsts = innerProduct(EigVec[i],psi[TwoN-1]);
	chi[TwoN] -= alpha[TwoN-1]*cconsts*EigVec[i];
      }

      // Subtract out projected psi_{N} from chi_{N-1} */
      chi[TwoN-1] += psi_proj * alpha[TwoN-1];
    }
                            
    // Complete the last component
    // chi(N) = chi(N) + beta_{N}*gamma_5*M*psi_{N}
    // The other two contributions
    //     chi[N] = mass * gamma[5]*psi[N] + alpha[N-1] * psi[N-1]
    // have been already calculated above 
   
    //  The contribution beta_{N}*gamma_5*M*psi_{N} is
    //   caluclated only if betaa_{N} != 0 */


    if( !isLastZeroP ) {
    
      // Complete psi_proj
      psi_proj += psi[TwoN];

      (*M)(tmp1, psi_proj, PLUS);
      pmscale = beta[TwoN]*scale_fac;
      tmp2 = Gamma(G5)*tmp1;
      chi[TwoN] += pmscale*tmp2;

    }

    getFermBC().modifyF(chi);

    END_CODE();
  }


  //! Derivative
  void 
  UnprecOvlapContFrac5DLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
					 const multi1d<LatticeFermion>& chi, 
					 const multi1d<LatticeFermion>& psi, 
					 enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    int G5 = Ns*Ns - 1;

    // This is the upper limit for the index of the 5th dimension, i.e.
    //  it runs from 0 to TwoN. (altogether TwoN + 1 numbers. This
    // makes sense since N5 is ALWAYS odd. )
    int TwoN = N5 - 1;

    // The sign of the diagonal terms.
    // We flip Hsign BEFORE using it 
    // This makes the first element always have a PLUS sign
    int Hsign=-1;

    /* N5 is ALWAYS ODD hence the explicit setting above */
    /*
      if( N5%2 == 0 ) {
      Hsign = 1;
      }
      else {
      Hsign = -1;
      }
    */

    // Run through all the pseudofermion fields:
    //   chi(0) = beta(0)*H*psi(0) + alpha(0)*psi(1)
    //   chi(n) = alpha(n-1)*psi(n-1)
    //              + (-)^n*beta(n)*H*psi(n) + alpha(n)*psi(n+1)
    //   chi(TwoN) = alpha(TwoN-1)*psi(TwoN-1)
    //                      + (-)^TwoN*beta(TwoN)*H*psi(TwoN)         
    
    LatticeFermion tmp1, tmp2;
    Real pmscale;

    switch (isign)
    {
    case PLUS:
      for(int n = 0; n < N5; ++n) 
      {
	if (n == N5-1 && isLastZeroP) continue;

	tmp2 = Gamma(G5)*chi[n];        // tmp2 = gamma_5 M chi

	// Scale factor and sign for the diagonal term proportional to H
	// The scale factor should be chosen in fermact call such
	//  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
	Hsign = -Hsign;
	pmscale = beta[n]*Hsign*scale_fac;
	tmp1 = pmscale*tmp2;

	M->deriv(ds_tmp, tmp1, psi[n], PLUS);
	ds_u += ds_tmp;
      }
      break;

    case MINUS:
      for(int n = 0; n < N5; ++n) 
      {
	if (n == N5-1 && isLastZeroP) continue;

	tmp2 = Gamma(G5)*psi[n];        // tmp2 = M^dag gamma_5 psi

	// Scale factor and sign for the diagonal term proportional to H
	// The scale factor should be chosen in conszolotarev5d_w.m such
	//  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
	Hsign = -Hsign;
	pmscale = beta[n]*Hsign*scale_fac;
	tmp1 = pmscale*tmp2;

	M->deriv(ds_tmp, chi[n], tmp1, MINUS);
	ds_u += ds_tmp;
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }

    // Last Component
    if(  NEig > 0 ) 
    {
      QDPIO::cerr << "contfrac5d deriv - projection not supported" << endl;
      QDP_abort(1);
    }

    getFermBC().zero(ds_u);

    END_CODE();
  }

}; // End Namespace Chroma


