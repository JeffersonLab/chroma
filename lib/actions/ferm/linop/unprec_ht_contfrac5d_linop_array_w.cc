// $Id: unprec_ht_contfrac5d_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned H_T kernel continued fraction (5D) operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_ht_contfrac5d_linop_array_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"


namespace Chroma 
{ 
  // Initialize
  void
  UnprecHTContFrac5DLinOpArray::init(Handle< FermState<T,P,Q> > state,
				     const Real& b5, const Real& c5)
  {
    scale_fac = b5 + c5;
    a5 = b5 - c5;

    D_w = new UnprecWilsonLinOp(state, Real(-OverMass));
    D_denum = new UnprecDWFTransfDenLinOp(a5, D_w);
    fbc = state->getFermBC();
  }


  //! Apply the operator onto a source vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi     Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecHTContFrac5DLinOpArray::operator() (multi1d<LatticeFermion>& chi,
					    const multi1d<LatticeFermion>& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);

    int G5 = Ns*Ns - 1;
    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

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
    
    LatticeFermion tmp1, tmp2;     moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);

    Real pmscale;

    for(int n = 0; n < TwoN; ++n) {
      (*D_w)(tmp1, psi[n], PLUS);   // tmp1 = D_w psi[n]
      tmp2 = Gamma(G5)*tmp1;        // tmp2 = gamma_5 D_w psi

      // Scale factor and sign for the diagonal term proportional to H
      // The scale factor should be chosen in conszolotarev5d_w.m such
      //  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
      Hsign = -Hsign;
      pmscale = beta[n]*Hsign*scale_fac;
      chi[n] = pmscale*tmp2;

      if( N5 > 1 ) { 
	(*D_denum)(tmp1, psi[n+1], msign);   // tmp1 = D_denum^dag psi[n+1]
	chi[n] += alpha[n] * tmp1;
      }

      if( n > 0 ) {
	(*D_denum)(tmp1, psi[n-1], msign);   // tmp1 = D_denum^dag psi[n-1]
	chi[n] += alpha[n-1]*tmp1;
      }
    }

    // Last Component
    // chi[N] = D_denum^msign*(mass*gamma5*psi[N] + alpha[N-1]*psi[N-1])
    tmp1 = Gamma(G5)*psi[TwoN];
    (*D_denum)(tmp2, tmp1, MINUS);   // NOTE FIXED SIGN HERE
    chi[TwoN] = mass*tmp2;

    (*D_denum)(tmp1, psi[TwoN-1], msign);   // tmp2 = alpha[TwoN-1]*D_denum^msign * psi[TwoN-1]
    chi[TwoN] += alpha[TwoN-1]*tmp1;

    // Complete the last component
    // chi(N) += beta_{N}*gamma_5*M*psi_{N}
    // The other two contributions
    //     chi[N] = mass * gamma[5]*psi[N] + alpha[N-1] * psi[N-1]
    // have been already calculated above 
   
    //  The contribution beta_{N}*gamma_5*M*psi_{N} is
    //   calculated only if beta_{N} != 0 */

    if( !isLastZeroP ) 
    {
      (*D_w)(tmp1, psi[TwoN], PLUS);
      pmscale = beta[TwoN]*scale_fac;
      tmp2 = Gamma(G5)*tmp1;
      chi[TwoN] += pmscale*tmp2;
    }

    getFermBC().modifyF(chi);

    END_CODE();
  }

}; // End Namespace Chroma


