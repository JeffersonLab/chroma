// $Id: prec_ovlap_contfrac5d_pv_linop_array_w.cc,v 1.3 2005-01-21 05:26:22 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Pauli-Villars Continued Fraction 5D
 */

#include "chromabase.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/prec_ovlap_contfrac5d_pv_linop_array_w.h"


namespace Chroma 
{ 
  EvenOddPrecOvlapContFrac5DPVLinOpArray::EvenOddPrecOvlapContFrac5DPVLinOpArray(
    Handle<const ConnectState> state,
    const Real& _m_q,
    const Real& _OverMass,
    int _N5,
    const Real& _scale_fac,
    const multi1d<Real>& _alpha,
    const multi1d<Real>& _beta,
    const bool _isLastZeroP ) :
    m_q(_m_q), OverMass(_OverMass), N5(_N5), scale_fac(_scale_fac), 
    alpha(_alpha), beta(_beta), isLastZeroP(_isLastZeroP)
  {
    START_CODE();

    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Ds(new WilsonDslash(state->getLinks()));
    Dslash  = Ds;  // Copy Handle -- M now owns dslash

    // Now compute some coefficients.
    // First the beta_tilde_i
    // Basically this is beta[n]*scale_fac
    // Now N5 is always odd. So the first Hsign is +
    // and the last one should also be
    // Hence at the end of this loop Hsign should be flipped from +->-
    beta_tilde.resize(N5);
    for(int i=0; i < N5; i++) 
    {
      beta_tilde[i] = beta[i]*scale_fac; 

      /*
	QDPIO::cout << "beta["<<i<<"]=" << beta[i]
	<< "  scale_fac=" << scale_fac 
	<< "  beta_tilde["<<i<<"]=" << beta_tilde[i] << endl;

      */
    }

    // Now the a_i's and b_i's
    a.resize(N5);
    for(int i=0; i < N5-1; i++) { 
      a[i] = beta_tilde[i]*(Nd + OverMass);
    }
    a[N5-1] = 1;         // CHECK THIS: WHAT IS THE NORM OF THE 1 IN THE PV TERM

    /*
      QDPIO::cout << "Nd + OverMass = " << Nd+ OverMass << endl;
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "a["<<i<<"]= " << a[i] << endl;
      }
    */

    END_CODE();
  }




  //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecOvlapContFrac5DPVLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
						    const multi1d<LatticeFermion>& psi, 
						    enum PlusMinus isign,
						    const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    // We don't care about the isign because our operator is Hermitian
    // Apply matrix
    //   [ A_0  B_0   0     ...                       ]  [ psi_0    ]
    //   [ B_0  A_1  B_1                   ...   ...  ]  [ psi_1    ]
    //   [  0   ...  ...     ...                 ...  ]  [ psi_2    ]
    //   [  ...    ...    0    B_N5-3  A_N5-2  B_N5-2 ]  [ psi_N5-2 ]
    //   [  ...    ...    ...   0      B_N5-2  A_N5-1 ]  [ psi_N5-1 ]

    // With A_i = gamma_5 a_i = a_i gamma_5
    // and  B_i = b_i I = alpha_i I

    int G5=Ns*Ns-1;

    // All in one shot
    for(int i=0; i < N5-1; i++) 
    {
      // A_{i} psi[i] = a_{i} g_5 psi[i]
      chi[i][rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] *= a[i];
    }

    chi[N5-1][rb[cb]] = psi[N5-1];

    END_CODE();
  }


  //! Apply the inverse even-even (odd-odd)
  /*!
   * \ingroup linop
   *
   * Here we apply the LDU decomposition
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecOvlapContFrac5DPVLinOpArray::applyDiagInv(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    Real coeff;

    const int G5 = Ns*Ns-1;

    // Invert diagonal piece  chi <- D^{-1} psi
    for(int i = 0; i < N5-1; i++) 
    {
      chi[i][rb[cb]] = Gamma(G5)*psi[i];
      coeff = Real(1)/a[i];
      chi[i][rb[cb]] *= coeff;
    }

    chi[N5-1][rb[cb]] = psi[N5-1];

    END_CODE();
  }


  //! Apply the off diagonal block
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )   	   (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void EvenOddPrecOvlapContFrac5DPVLinOpArray::applyOffDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi,
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    LatticeFermion tmp;
    Real coeff;
    int G5 = Ns*Ns-1;

    // Optimisation... do up to the penultimate block...
    for(int i=0; i < N5-1; i++) 
    {
      // CB is CB of TARGET
      // gamma_5 Dslash is hermitian so I can ignore isign

      // Apply g5 Dslash
      Dslash->apply(tmp, psi[i], PLUS, cb);
      chi[i][rb[cb]] = Gamma(G5)*tmp;

      // Multiply coefficient
      coeff = -Real(0.5)*beta_tilde[i];

      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[i][rb[cb]] *= coeff;
    }

 
    // Only do last block if beta_tilde[i] is not zero
//    if( !isLastZeroP ) {
//    }
  
    chi[N5-1][rb[cb]] = zero;

    END_CODE();
  }


  // Derivative of even-odd linop component
  /* 
   * This is a copy of the above applyOffDiag with the D.apply(...) replaced
   * by  D.deriv(ds_tmp,...) like calls.
   */
  void 
  EvenOddPrecOvlapContFrac5DPVLinOpArray::applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u,
							    const multi1d<LatticeFermion>& chi, 
							    const multi1d<LatticeFermion>& psi, 
							    enum PlusMinus isign,
							    int cb) const 
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
						   
    LatticeFermion tmp;
    Real coeff;
    int G5 = Ns*Ns-1;

    switch (isign)
    {
    case PLUS:
      // Optimisation... do up to the penultimate block...
      for(int i=0; i < N5-1; i++) 
      {
	// CB is CB of TARGET
	// consider gamma_5 Dslash
	tmp[rb[cb]] = Gamma(G5)*chi[i];

	// Multiply coefficient
	coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[cb]] *= coeff;

	// Apply g5 Dslash
	Dslash->deriv(ds_tmp, tmp, psi[i], PLUS, cb);
	ds_u += ds_tmp;
      }
      break;

    case MINUS:
      // Optimisation... do up to the penultimate block...
      for(int i=0; i < N5-1; i++) 
      {
	// CB is CB of TARGET
	// consider Dslash^dag gamma_5
	tmp[rb[1-cb]] = Gamma(G5)*psi[i];

	// Multiply coefficient
	coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[1-cb]] *= coeff;

	// Apply g5 Dslash
	Dslash->deriv(ds_tmp, chi[i], tmp, MINUS, cb);
	ds_u += ds_tmp;
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }
    

//    if( !isLastZeroP ) 
//    {
//      // This term does not contribute to PV
//    }

 
    END_CODE();
  }


}; // End Namespace Chroma


