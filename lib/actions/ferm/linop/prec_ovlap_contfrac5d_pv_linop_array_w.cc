// $Id: prec_ovlap_contfrac5d_pv_linop_array_w.cc,v 1.7 2005-03-02 18:32:05 bjoo Exp $
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
    // Basically this is beta[n]*Hsign*scale_fac
    // Now N5 is always odd. So the first Hsign is +
    // and the last one should also be
    // Hence at the end of this loop Hsign should be flipped from +->-
    beta_tilde.resize(N5);
    int Hsign = 1;
    for(int i=0; i < N5; i++) 
    {
      // Flip Hsign
      beta_tilde[i] = beta[i]*Hsign*scale_fac; 

      /*
	QDPIO::cout << "beta["<<i<<"]=" << beta[i]
	<< "  scale_fac=" << scale_fac 
	<< "  beta_tilde["<<i<<"]=" << beta_tilde[i] << endl;

      */
      Hsign = -Hsign;
    }

    // Sanity Check
    if ( Hsign > 0 ) {
      QDPIO::cerr << "Something is wrong. At the end of this loop"
		  << " Hsign should be -ve" << endl;
    }

    // Now the a_i's and b_i's
    a.resize(N5);
    for(int i=0; i < N5-1; i++) { 
      a[i] = beta_tilde[i]*(Nd - OverMass);
    }
    a[N5-1] = 1;         // CHECK THIS: WHAT IS THE NORM OF THE 1 IN THE PV TERM

    /*
      QDPIO::cout << "Nd - OverMass = " << Nd - OverMass << endl;
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "a["<<i<<"]= " << a[i] << endl;
      }
    */

    // Now the d-s
    d.resize(N5);
    d[0] = a[0];
    d[N5-1] = a[N5-1];
    for(int i=1; i < N5-1; i++) { 
      d[i] = a[i] - (alpha[i-1]*alpha[i-1])/d[i-1];
    }
    
    /*
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "d["<<i<<"]=" << d[i] << endl;
      }
    */

    // Now the u-s
    u.resize(N5-1);
    u[N5-2] = 0.0;
    for(int i=0; i < N5-2; i++) { 
      u[i] = alpha[i]/d[i];
    }

    /*
      for(int i=0; i < N5-1; i++) { 
      QDPIO::cout << "u["<<i<<"] = " << u[i] << endl;
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

    if( chi.size() != N5 ) chi.resize(N5);

    // We don't care about the isign because our operator is Hermitian
    // Apply matrix
    //   [ A_0  B_0   0     ...                       ]  [ psi_0    ]
    //   [ B_0  A_1  B_1                   ...   ...  ]  [ psi_1    ]
    //   [  0   ...  ...     ...                 ...  ]  [ psi_2    ]
    //   [  ...    ...    0    B_N5-3  A_N5-2  B_N5-2 ]  [ psi_N5-2 ]
    //   [  ...    ...    ...   0      B_N5-2  A_N5-1 ]  [ psi_N5-1 ]

    // With A_i = gamma_5 a_i = a_i gamma_5
    // and  B_i = b_i I = alpha_i I

    LatticeFermion tmp;
    int G5=Ns*Ns-1;

    // First 0ne 
    tmp[rb[cb]] = Gamma(G5)*psi[0];
    chi[0][rb[cb]] = a[0]*tmp;
    if( N5 > 1 ) { 
      chi[0][rb[cb]] += alpha[0]*psi[1];
    }

    // All the rest
    for(int i=1; i < N5-1; i++) 
    { 
      // B_{i-1} psi_[i-1]
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      tmp[rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] += a[i]*tmp;

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-2) 
	chi[i][rb[cb]] += alpha[i]*psi[i+1];
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

    if( chi.size() != N5 )  chi.resize(N5);

    multi1d<LatticeFermion> y(N5);

    LatticeFermion tmp;
    Real coeff;

    const int G5 = Ns*Ns-1;

    // Solve L y = psi
    y[0][rb[cb]] = psi[0];

    for(int i = 1; i < N5-1; i++) 
    {
      tmp[rb[cb]] = Gamma(G5)*y[i-1];
      y[i][rb[cb]] = psi[i] - u[i-1]*tmp;
    } 

    // Invert diagonal piece  y <- D^{-1} y
    for(int i = 0; i < N5-1; i++) 
    {
      tmp[rb[cb]] = Gamma(G5)*y[i];
      coeff = Real(1)/d[i];
      y[i][rb[cb]] = coeff*tmp;
    }

    // Backsubstitute U chi = y
    chi[N5-1][rb[cb]] = psi[N5-1];

    for(int i = N5-2; i >= 0; i--) {
      tmp = Gamma(G5)*chi[i+1];
      chi[i][rb[cb]] = y[i] - u[i]*tmp;
    }

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

    if( chi.size() != N5 )  chi.resize(N5);

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


