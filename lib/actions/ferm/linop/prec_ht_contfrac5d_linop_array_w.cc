// $Id: prec_ht_contfrac5d_linop_array_w.cc,v 1.2 2005-02-13 18:16:28 edwards Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/prec_ht_contfrac5d_linop_array_w.h"


namespace Chroma 
{ 
  EvenOddPrecHtContFrac5DLinOpArray::EvenOddPrecHtContFrac5DLinOpArray(
    Handle<const ConnectState> state,
    const Real& _m_q,
    const Real& _OverMass,
    int _N5,
    const Real& _scale_fac,
    const multi1d<Real>& _alpha,
    const multi1d<Real>& _beta,
    const Real& b5_,
    const Real& c5_,
    const bool _isLastZeroP ) :
    m_q(_m_q), OverMass(_OverMass), N5(_N5), scale_fac(_scale_fac), 
    alpha(_alpha), beta(_beta), isLastZeroP(_isLastZeroP), b5(b5_), c5(c5_)
  {
    START_CODE();

    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Ds(new WilsonDslash(state->getLinks()));
    Dslash  = Ds;  // Copy Handle -- M now owns dslash

    // The mass ratio
    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);

    // f_+ = (b5 + c5)
    f_plus = b5_ + c5_;

    // f_- = (b5 - c5)
    f_minus = b5_ - c5_;

    // Now compute some coefficients.
    // First the beta_tilde_i
    // Basically this is beta[n]*Hsign*scale_fac
    // Now N5 is always odd. So the first Hsign is +
    // and the last one should also be
    // Hence at the end of this loop Hsign should be flipped from +->-
    beta_tilde.resize(N5);
    int Hsign = 1;
    for(int i=0; i < N5; i++) { 
      
      beta_tilde[i] = beta[i]*Hsign*scale_fac*f_plus; 

      // Flip Hsign
      Hsign = -Hsign;

   
    }

    // Sanity Check
    if ( Hsign > 0 ) {
      QDPIO::cerr << "Something is wrong. At the end of this loop"
		  << " Hsign should be -ve" << endl;
    }

    alpha_tilde.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 

      alpha_tilde[i] = alpha[i]*(Real(2) + f_minus*(Nd - OverMass));
      
    }


    // Now the a_i's and b_i's
    a.resize(N5);
    for(int i=0; i < N5-1; i++) { 
      a[i] = beta_tilde[i]*(Nd - OverMass);
    }
    a[N5-1] = mass*(Real(2) + f_minus*(Nd - OverMass))
      + (beta_tilde[N5-1]*(Nd - OverMass));

    /*
      QDPIO::cout << "Nd - OverMass = " << Nd - OverMass << endl;
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "a["<<i<<"]= " << a[i] << endl;
      }
    */
    // Now the d-s
    d.resize(N5);
    d[0] = a[0];
    for(int i=1; i < N5; i++) { 
      d[i] = a[i] - (alpha_tilde[i-1]*alpha_tilde[i-1])/d[i-1];
    }
    
    /*
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "d["<<i<<"]=" << d[i] << endl;
      }
    */

    // Now the u-s
    u.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 
      u[i] = alpha_tilde[i]/d[i];
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
  EvenOddPrecHtContFrac5DLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
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
    // and  B_i = b_i I = alpha_tilde_i I

    LatticeFermion tmp;
    int G5=Ns*Ns-1;

    // First 0ne 
    tmp[rb[cb]] = Gamma(G5)*psi[0];
    chi[0][rb[cb]] = a[0]*tmp;
    if( N5 > 1 ) { 
      chi[0][rb[cb]] += alpha_tilde[0]*psi[1];
    }

    // All the rest
    for(int i=1; i < N5; i++) { 

      // B_{i-1} psi_[i-1]
      chi[i][rb[cb]] = alpha_tilde[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      tmp[rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] += a[i]*tmp;

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-1) {
	chi[i][rb[cb]] += alpha_tilde[i]*psi[i+1];
      }
    }

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
  EvenOddPrecHtContFrac5DLinOpArray::applyDiagInv(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    multi1d<LatticeFermion> y(N5);

    LatticeFermion tmp;
    Real coeff;

    const int G5 = Ns*Ns-1;

    // Solve L y = psi
    y[0][rb[cb]] = psi[0];

    for(int i = 1; i < N5; i++) { 
      tmp[rb[cb]] = Gamma(G5)*y[i-1];
      y[i][rb[cb]] = psi[i] - u[i-1]*tmp;
    } 

    // Invert diagonal piece  y <- D^{-1} y
    for(int i = 0; i < N5; i++) { 
      tmp[rb[cb]] = Gamma(G5)*y[i];
      coeff = Real(1)/d[i];
      y[i][rb[cb]] = coeff*tmp;
    }

    // Backsubstitute U chi = y
    chi[N5-1][rb[cb]] = y[N5-1];

    for(int i = N5-2; i >= 0; i--) {
      tmp = Gamma(G5)*chi[i+1];
      chi[i][rb[cb]] = y[i] - u[i]*tmp;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }

  //! Apply the off diagonal block
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )   	   (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void EvenOddPrecHtContFrac5DLinOpArray::applyOffDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi,
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);
    LatticeFermion tmp;


    int G5 = Ns*Ns-1;
    Real ftmp;

    switch(isign) { 
    case PLUS : 
      {
	// [ A_0    B_0   0    .......        ]
        // [ B_0    A_1   B_1  .......        ]
	// [  0     B_1   A_2  B_2 .......    ]
	// [  ......................... B_N5-2]
        // [              0    B_N5-2   A_N5-1]
	//
	//
	//  WIth 
	// A[i]    = beta_tilde[i] gamma_5 (-1/2) Dslash        i < N5-1
	// A[N5-1] = mass*f_minus (-1/2) Dslash^dagger gamma_5 
        //        +    beta_tilde[N5-1] gamma_5 (-1/2) Dslash
	//
	// (beta_tilde has f_plus folded into it already)
	// 
	// B_i = alpha_i * f_minus (-1/2) D^{\dagger} 

	//
	// Evaluate (-1/2) Dslash and (-1/2) Dslash^dagger up front
	multi1d<LatticeFermion> D_psi(N5);
	multi1d<LatticeFermion> Ddag_psi(N5);
	Real ftmp_mhalf = Real(-0.5);

	for(int i=0; i < N5; i++) {
	  D_psi[i] = zero;
	  Dslash->apply(D_psi[i], psi[i], PLUS, cb);
	  D_psi[i][rb[cb]] *= ftmp_mhalf;
	}

	for(int i=0; i < N5; i++) { 
	  Ddag_psi[i] = zero;
	  Dslash->apply(Ddag_psi[i], psi[i], MINUS, cb);
	  Ddag_psi[i][rb[cb]] *= ftmp_mhalf;
	}
	
	
	// First component only 2 terms 

	// term1 : A_0 psi 0
	chi[0][rb[cb]]  = Gamma(G5)*D_psi[0];
	chi[0][rb[cb]] *= beta_tilde[0];     

	// term2: B_0 psi_1 
	ftmp = alpha[0]*f_minus;
	chi[0][rb[cb]] += ftmp*Ddag_psi[1];  


	// Other components
	for(int i=1; i < N5-1; i++) {
		  
	  // term1: A[i-1] psi[i-1]
	  ftmp = alpha[i-1]*f_minus;
	  chi[i][rb[cb]] = ftmp*Ddag_psi[i-1];
	  
	  // term2: B[i] psi[i]
	  ftmp = beta_tilde[i];
	  tmp[rb[cb]] = Gamma(G5)*D_psi[i];
	  chi[i][rb[cb]] += ftmp*tmp;
	  
	  // term 3 A[i] psi[i+1]
	  ftmp = alpha[i]*f_minus;
	  chi[i][rb[cb]] += ftmp*Ddag_psi[i+1];
	  
	}

	// Last component
	// gamma_5 is on the unpleasant side of Dslash^dagger
	// but can convert to a gamma_5 D since
	// gamma_5 Dslash = Dslash^dagger Gamma_5 

	// Start with the: mass fminus (-1/2)  Dslash^dagger gamma_5 psi[N5-1]
        //          = mass fminus gamma_5 (-1/2) Dslash psi[N5-1];   
	tmp[rb[cb]] = Gamma(G5)*D_psi[N5-1];
	ftmp = mass * f_minus;
	chi[N5-1][rb[cb]] = ftmp*tmp;

	// Term number 2: alpha[N5-2]*f_minus*(-1/2) Dslash^dagger psi[N5-2]
	ftmp = alpha[N5-2]*f_minus;
	chi[N5-1][rb[cb]] += ftmp*Ddag_psi[N5-2];

	// Term involving beta is last is not zero
	if( !isLastZeroP ) { 
	  tmp[rb[cb]] = Gamma(G5)*D_psi[N5-1];
	  chi[N5-1][rb[cb]] += beta_tilde[N5-1]*tmp;
	}

      }
      break;

    case MINUS:
      {
	multi1d<LatticeFermion> D_psi(N5);
	Real ftmp_mhalf = Real(-0.5);
	
	for(int i=0; i < N5; i++) {
	  D_psi[i] = zero;
	  Dslash->apply(D_psi[i], psi[i], PLUS, cb);
	  D_psi[i][rb[cb]] *= ftmp_mhalf;
	}

	// First bit
	// Bits involving beta_tilde are hermitian and do not change
	chi[0][rb[cb]] = Gamma(G5)*D_psi[0];
	chi[0][rb[cb]] *= beta_tilde[0];

	// Bits involving alpha: Dslash^dagger -> Dslash
	ftmp = alpha[0]*f_minus;
	chi[0][rb[cb]] += ftmp*D_psi[1];


	for(int i=1; i < N5-1; i++) {
		  
	  ftmp = alpha[i-1]*f_minus;
	  chi[i][rb[cb]] = ftmp*D_psi[i-1];
	  
	  ftmp = beta_tilde[i];
	  tmp[rb[cb]] = Gamma(G5)*D_psi[i];
	  chi[i][rb[cb]] += ftmp*tmp;

	  ftmp = alpha[i]*f_minus;
	  chi[i][rb[cb]] += ftmp*D_psi[i+1];
	  
	}

	tmp[rb[cb]] = Gamma(G5) * D_psi[N5-1];
	ftmp = mass*f_minus;
	chi[N5-1][rb[cb]] = ftmp*tmp;

	ftmp = alpha[N5-2]*f_minus;
	chi[N5-1][rb[cb]] += ftmp*D_psi[N5-2];

	if( !isLastZeroP ) { 
	  tmp[rb[cb]] = Gamma(G5)*D_psi[N5-1];
	  chi[N5-1][rb[cb]] += beta_tilde[N5-1]*tmp;
	}



      }
      break;
    default:
      {
	QDPIO::cerr << "Should never reach here. Isign is only 2 valued" << endl;
	QDP_abort(1);
      }
      break;
    }
	  
    END_CODE();
  }


  // Derivative of even-odd linop component
  /* 
   * This is a copy of the above applyOffDiag with the D.apply(...) replaced
   * by  D.deriv(ds_tmp,...) like calls.
   */
  void 
  EvenOddPrecHtContFrac5DLinOpArray::applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u,
							  const multi1d<LatticeFermion>& chi, 
							  const multi1d<LatticeFermion>& psi, 
							  enum PlusMinus isign,
							  int cb) const 
  {
    START_CODE();

#if 0
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
      for(int i=0; i < N5; i++) 
      {
	if (i == N5-1 && isLastZeroP) continue;

	// CB is CB of TARGET
	// consider case of gamma_5 Dslash
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
      for(int i=0; i < N5; i++) 
      {
	if (i == N5-1 && isLastZeroP) continue;

	// CB is CB of TARGET
	// consider case of Dslash^dag gamma_5
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

#else
    QDPIO::cerr << "NOt yet implemented " << endl;
    QDP_abort(1);
#endif

    END_CODE();
  }




  // THIS IS AN OPTIMIZED VERSION OF THE DERIVATIVE
  void 
  EvenOddPrecHtContFrac5DLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u,
					      const multi1d<LatticeFermion>& chi, 
					      const multi1d<LatticeFermion>& psi, 
					      enum PlusMinus isign) const
  {
    START_CODE();
#if 0
    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

    ds_u.resize(Nd);

    multi1d<LatticeFermion>  tmp1, tmp2, tmp3;
    multi1d<LatticeColorMatrix> ds_tmp;

    //  ds_u   =  chi^dag * D'_oe * Ainv_ee * D_eo * psi_o
    evenOddLinOp(tmp1, psi, isign);
    evenEvenInvLinOp(tmp2, tmp1, isign);
    derivOddEvenLinOp(ds_u, chi, tmp2, isign);

    //  ds_u  +=  chi^dag * D_oe * Ainv_ee * D'_eo * psi_o
    evenOddLinOp(tmp1, chi, msign);
    evenEvenInvLinOp(tmp3, tmp1, msign);
    derivEvenOddLinOp(ds_tmp, tmp3, psi, isign);
    ds_u += ds_tmp;
    
    for(int mu=0; mu < Nd; mu++)
      ds_u[mu] *= Real(-1);
#else
    QDPIO::cerr << "Not yet implemented " << endl;
    QDP_abort(1);
#endif

    END_CODE();
  }

}; // End Namespace Chroma


