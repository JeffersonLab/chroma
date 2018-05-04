/*! \file
 *  \brief Unpreconditioned Pauli-Villars Continued Fraction 5D
 */

#include "chromabase.h"
#include "linearop.h"

#include "actions/ferm/linop/unprec_ovlap_contfrac5d_pv_linop_array_w.h"

namespace Chroma 
{ 
  //! Apply the operator onto a source std::vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecOvlapContFrac5DPVLinOpArray::operator() (multi1d<LatticeFermion>& chi,
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

    // The sign of the diagonal terms.
    // We flip Hsign BEFORE using it 
    // This makes the first element always have a PLUS sign
    int Hsign=-1;

    // Run through all the pseudofermion fields:
    //   chi(n) = beta(n)*H*psi(0)
    //   chi(TwoN) = beta(TwoN)*H*psi(TwoN)         
    
    LatticeFermion tmp1, tmp2;
    Real pmscale;

    for(int n = 0; n < TwoN; ++n) 
    {
      (*M)(tmp1, psi[n], PLUS);      // tmp1 = M psi[n]
      tmp2 = Gamma(G5)*tmp1;        // tmp2 = gamma_5 M psi

      // Scale factor and sign for the diagonal term proportional to H
      // The scale factor should be chosen in the fermact call such that
      //  that scale_fac * gamma5 * M has eigenvalues between -1 and 1  
      Hsign = -Hsign;
      pmscale = beta[n]*Hsign*scale_fac;
      chi[n] = pmscale*tmp2;

      if( n < TwoN-1 ) { 
	chi[n] += alpha[n] * psi[n+1];
      }

      if( n > 0 ) {
	chi[n] += alpha[n-1]*psi[n-1];
      }
    }

    // Last Component
    // chi[N] = psi[N]
    chi[TwoN] = psi[TwoN];

    // Project out eigenvectors from Source if desired 
    if(  NEig > 0 ) 
    {
      QDPIO::cerr << "project in contfrac5d PV not supported" << std::endl;
      QDP_abort(1);
    }
                            
    // Complete the last component
    // chi(N) = chi(N) + beta_{N}*psi_{N}
    //  The contribution beta_{N}*gamma_5*M*psi_{N} is
    //   caluclated only if betaa_{N} != 0 */

//    if( !isLastZeroP ) 
//    {
//      // This term does not contribute to PV
//    }

    getFermBC().modifyF(chi);

    END_CODE();
  }


  //! Derivative
  void 
  UnprecOvlapContFrac5DPVLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
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

    // Run through all the pseudofermion fields:
    //   chi(n) = beta(n)*H*psi(0)
    //   chi(TwoN) = beta(TwoN)*H*psi(TwoN)         
    
    LatticeFermion tmp1, tmp2;
    Real pmscale;

    switch (isign)
    {
    case PLUS:
      for(int n = 0; n < TwoN; ++n) 
      {
	tmp2 = Gamma(G5)*chi[n];        // tmp2 = gamma_5 M chi

	// Scale factor and sign for the diagonal term proportional to H
	// The scale factor should be chosen in fermact call such
	//  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
	Hsign = -Hsign;
	pmscale = beta[n]*Hsign*scale_fac;
	tmp1 = pmscale*tmp2;

	M->deriv(ds_tmp, tmp1, psi[n], PLUS);      // tmp1 = M psi[n]
	ds_u += ds_tmp;
      }
      break;

    case MINUS:
      for(int n = 0; n < TwoN; ++n) 
      {
	tmp2 = Gamma(G5)*psi[n];        // tmp2 = M^dag gamma_5 psi

	// Scale factor and sign for the diagonal term proportional to H
	// The scale factor should be chosen in fermact call such
	//  that scale_fac * gamma5 * M has eigenvalues between -1 and 1 
	Hsign = -Hsign;
	pmscale = beta[n]*Hsign*scale_fac;
	tmp1 = pmscale*tmp2;

	M->deriv(ds_tmp, chi[n], tmp1, MINUS);      // tmp1 = M psi[n]
	ds_u += ds_tmp;
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }

    // Project out eigenvectors from Source if desired 
    if(  NEig > 0 ) 
    {
      QDPIO::cerr << "project in contfrac5d PV not supported" << std::endl;
      QDP_abort(1);
    }
                            
    getFermBC().zero(ds_u);

    END_CODE();
  }


} // End Namespace Chroma


