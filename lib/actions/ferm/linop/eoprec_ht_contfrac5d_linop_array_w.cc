// $Id: eoprec_ht_contfrac5d_linop_array_w.cc,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/eoprec_ht_contfrac5d_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  // Full constructor
  EvenOddPrecHtContFrac5DLinOpArray::EvenOddPrecHtContFrac5DLinOpArray(
    Handle< FermState<T,P,Q> > fs,
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

    Dslash.create(fs, _N5);

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
    invd.resize(N5);

    d[0] = a[0];
    invd[0] = Real(1)/d[0];
    for(int i=1; i < N5; i++) { 
      d[i] = a[i] - (alpha_tilde[i-1]*alpha_tilde[i-1])/d[i-1];
      invd[i] = Real(1)/d[i];
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
   *
   * 6*Nc*Ns + (N5-2)*10*Nc*Ns + 6Nc*Ns = 12*NcNs + 10*N5NcNs - 20*NcNs
   * = (10N5 - 8)NcNs
   */
  void 
  EvenOddPrecHtContFrac5DLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
					       const multi1d<LatticeFermion>& psi, 
					       enum PlusMinus isign,
					       const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);

    // We don't care about the isign because our operator is Hermitian
    // Apply matrix
    //   [ A_0  B_0   0     ...                       ]  [ psi_0    ]
    //   [ B_0  A_1  B_1                   ...   ...  ]  [ psi_1    ]
    //   [  0   ...  ...     ...                 ...  ]  [ psi_2    ]
    //   [  ...    ...    0    B_N5-3  A_N5-2  B_N5-2 ]  [ psi_N5-2 ]
    //   [  ...    ...    ...   0      B_N5-2  A_N5-1 ]  [ psi_N5-1 ]

    // With A_i = gamma_5 a_i = a_i gamma_5
    // and  B_i = b_i I = alpha_tilde_i I

    // LatticeFermion tmp;
    int G5=Ns*Ns-1;

    // First 0ne 
    //tmp[rb[cb]] = Gamma(G5)*psi[0];

    // 2*Nc*Ns flops/cbsite
    chi[0][rb[cb]] = a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    if( N5 > 1 ) { 

      // 4*Nc*Ns flops/cbsite
      chi[0][rb[cb]] += alpha_tilde[0]*psi[1];
    }

    // All the rest
    for(int i=1; i < N5; i++) { 

      // B_{i-1} psi_[i-1]
      // chi[i][rb[cb]] = alpha_tilde[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      // tmp[rb[cb]] = Gamma(G5)*psi[i];
      // 6 NcNs flops/cbsite
      chi[i][rb[cb]] = alpha_tilde[i-1]*psi[i-1] + a[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-1) {
	// 4NcNs flops/cbsite
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
   *
   * (N5-2)*6NcNs + 2NcNs + (N5-1)*4NcNs
   * = 6N5 NcNs - 12 NcNs + 2NcNs + 4 N5 NcNs - 4NcNs
   * = 10*N5*Nc*Ns -14NcNs
   */
  void 
  EvenOddPrecHtContFrac5DLinOpArray::applyDiagInv(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    multi1d<LatticeFermion> y(N5);
    moveToFastMemoryHint(y);

    //    LatticeFermion tmp;
    //    Real coeff;

    const int G5 = Ns*Ns-1;

    // Solve L y = psi
    //       y=D^{-1}y 
    // together in one loop
    y[0][rb[cb]] = psi[0];

    /* (N5 - 2)* 6NcNs */
    for(int i = 1; i < N5; i++) { 
      y[i][rb[cb]] = psi[i] - u[i-1]*(GammaConst<Ns,Ns*Ns-1>()*y[i-1]);
      y[i-1][rb[cb]] = invd[i-1]*(GammaConst<Ns,Ns*Ns-1>()*y[i-1]);
    } 
    /* 2NcNs */
    y[N5-1][rb[cb]] = invd[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*y[N5-1]);

    // Backsubstitute U chi = y
    chi[N5-1][rb[cb]] = y[N5-1];

    // (N5-1)*4NcNs
    for(int i = N5-2; i >= 0; i--) {
      // tmp = Gamma(G5)*chi[i+1];
      chi[i][rb[cb]] = y[i] - u[i]*(GammaConst<Ns,Ns*Ns-1>()*chi[i+1]);
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
   *
   * 10*(N5-2)*NcNs+12*NcNs+ N5*1320
   * (10N5 - 8)*NcNs + N5*1320
   */
  void EvenOddPrecHtContFrac5DLinOpArray::applyOffDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi,
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);


    int G5 = Ns*Ns-1;
    Real ftmp;
    Real ftmp2;
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
      //         = Dslash^dagger [ (-1/2) gamma_5 beta_tilde[i] ]

      // A[N5-1] = mass*f_minus (-1/2) Dslash^dagger gamma_5 
      //        +    beta_tilde[N5-1] gamma_5 (-1/2) Dslash
      //        = Dslash^dagger [ (-1/2)gamma_g5{
      //                         mass*f_minus + beta_tilde[N5-1] }
 
      // (beta_tilde has f_plus folded into it already)
      // 
      // B_i = alpha_i * f_minus (-1/2) D^{\dagger} 
      //     = D^{dagger} [ (-1/2)alpha_i*f_minus ]

      // Work everyhting out first. Apply D^dagger at the end.
	

      multi1d<LatticeFermion> tmp5(N5);
      moveToFastMemoryHint(tmp5);
      Real coeff_1, coeff_2, coeff_3;
      int otherCB = (cb + 1)%2;
	
      // First comomnent A_0 psi_0 + B_0 psi_1
      // 
      // without D^dagger we get:
      //
      // [ (-1/2)gamma_g5 beta_tilde[0] psi_0 + (-1/2)alpha_0 f_minus psi_1 ]
      // = (-1/2) [ gamma_g5 beta_tilde[0] psi_0 + alpha_0 f_minus_psi_1 ]
	

#if 0 
      tmp[rb[otherCB]] = Gamma(G5)*psi[0];
      coeff_1 = Real(-0.5)*beta_tilde[0];
      coeff_2 = Real(-0.5)*alpha[0]*f_minus;

      tmp5[0][rb[otherCB]] = coeff_1*tmp + coeff_2*psi[1];
#else
      coeff_1 = Real(-0.5)*beta_tilde[0];
      coeff_2 = Real(-0.5)*alpha[0]*f_minus;

      // 6NcNs flops
      tmp5[0][rb[otherCB]] = coeff_2*psi[1] + coeff_1*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
#endif	
      // N5-2*(10NcNs)
      for(int i=1; i < N5-1; i++) { 

	// i=1 .. N5-2
	//
	// B_{i-1} psi_{i-1} + A_i psi_i + B_{i} psi_{i+1}
	//   (-1/2)alpha_{i-1} f_minus psi_{i-1}
	// + (-1/2)alpha_{i} f_minus psi_{i+1}
	// + (-1/2)gamma_g5 beta_tilde[i] psi_i
	coeff_1 = Real(-0.5)*alpha[i-1]*f_minus;
	coeff_2 = Real(-0.5)*alpha[i]*f_minus;
	coeff_3 = Real(-0.5)*beta_tilde[i];
	//	  tmp[rb[otherCB]] = Gamma(G5)*psi[i];
	  
	tmp5[i][rb[otherCB]] = coeff_1*psi[i-1] + coeff_3*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);
	//tmp5[i][rb[otherCB]] += coeff_3*tmp;
	tmp5[i][rb[otherCB]] += coeff_2*psi[i+1];
      }

      // i=N5-1
      //
      // B_{i-1} psi_{i-1} + A_[N5-1] psi_[N5-1]
      // B_{i-1} psi_{i-1} = (-1/2) alpha[N5-1]*f_minus
	
      // A_[N5-1] psi[N5-1] = (-1/2)gamma_g5{
      //                         mass*f_minus + beta_tilde[N5-1] }
      coeff_1 = Real(-0.5)*alpha[N5-2]*f_minus;
      coeff_2 = Real(-0.5)*(mass*f_minus + beta_tilde[N5-1]);
      //tmp[rb[otherCB]] = Gamma(G5)*psi[N5-1];

      // 6NcNs
      tmp5[N5-1][rb[otherCB]] = coeff_1 * psi[N5-2] + coeff_2*(GammaConst<Ns,Ns*Ns-1>()*psi[N5-1]);

      // Now apply Dslash^{DAGGER} !!!!! 
      //
      // (the dagger comes from the denominator. Along the diagonal
      //  we have gamma_5 D = D^{dagger} gamma_g5
      // This could be done with a vec op. 
      // Alternatively, I could infuse it with other loops
      // don't know what's best.
      // for(int i=0; i < N5; i++) {
      //  Dslash->apply(chi[i], tmp5[i], MINUS, cb);
      //}
      Dslash.apply(chi,tmp5,MINUS, cb);

    }
    break;

    case MINUS:
    {
      multi1d<LatticeFermion> D_psi(N5); moveToFastMemoryHint(D_psi);
      Real ftmp_mhalf = Real(-0.5);

      Dslash.apply(D_psi, psi, PLUS, cb);

      for(int i=0; i < N5; i++) {
	// Do I need this? D_psi[i] = zero; Seemingly not!

	//	  Dslash.apply(D_psi[i], psi[i], PLUS, cb);
	D_psi[i][rb[cb]] *= ftmp_mhalf;
      }

      // First bit
      // Bits involving beta_tilde are hermitian and do not change
      //	chi[0][rb[cb]] = Gamma(G5)*D_psi[0];
      //      chi[0][rb[cb]] *= beta_tilde[0];
      ftmp = alpha[0]*f_minus;
      chi[0][rb[cb]] = ftmp*D_psi[1] + beta_tilde[0]*(GammaConst<Ns,Ns*Ns-1>()*D_psi[0]);

      for(int i=1; i < N5-1; i++) {
		  
	ftmp = alpha[i-1]*f_minus;
	ftmp2 = alpha[i]*f_minus;
	chi[i][rb[cb]] = ftmp*D_psi[i-1] + beta_tilde[i]*(GammaConst<Ns,Ns*Ns-1>()*D_psi[i]);
	  
	//	  ftmp = beta_tilde[i];
	//      tmp[rb[cb]] = Gamma(G5)*D_psi[i];
	// chi[i][rb[cb]] += ftmp*tmp;


	chi[i][rb[cb]] += ftmp2*D_psi[i+1];
	  
      }

      // tmp[rb[cb]] = Gamma(G5) * D_psi[N5-1];
      ftmp = mass*f_minus;
      ftmp2 = alpha[N5-2]*f_minus;
      chi[N5-1][rb[cb]] = ftmp2*D_psi[N5-2] + ftmp*(GammaConst<Ns,Ns*Ns-1>()*D_psi[N5-1]);




      if( !isLastZeroP ) { 
	LatticeFermion tmp; moveToFastMemoryHint(tmp);
	tmp[rb[cb]] = beta_tilde[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*D_psi[N5-1]);
	chi[N5-1][rb[cb]] += tmp;
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
						   
    LatticeFermion tmp; moveToFastMemoryHint(tmp);
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
	Dslash.deriv(ds_tmp, tmp, psi[i], PLUS, cb);
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
	Dslash.deriv(ds_tmp, chi[i], tmp, MINUS, cb);
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

    getFermBC().zero(ds_u);

    END_CODE();
  }

} // End Namespace Chroma


