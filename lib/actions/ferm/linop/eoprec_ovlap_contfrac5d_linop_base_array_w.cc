// $Id: eoprec_ovlap_contfrac5d_linop_base_array_w.cc,v 3.1 2006-10-19 16:01:30 edwards Exp $
/*! \file
 *  \brief Base class for even-odd prec. 5D continued fraction linop
 */

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_linop_base_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  EvenOddPrecOvlapContFrac5DLinOpBaseArray::EvenOddPrecOvlapContFrac5DLinOpBaseArray(
    Handle< FermState<T,P,Q> > state,
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

    int dslash_length;
    if( isLastZeroP ) {
      dslash_length = N5-1;
    }
    else {
      dslash_length = N5;
    }

    Dslash.create(state,dslash_length);

    // The mass ratio
    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);

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
	<< "  Hsign=" << Hsign
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
    a[N5-1] = mass + (beta_tilde[N5-1]*(Nd - OverMass));

    /*
      QDPIO::cout << "Nd - OverMass = " << Nd- OverMass << endl;
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "a["<<i<<"]= " << a[i] << endl;
      }
    */
    // Now the d-s
    multi1d<Real> d(N5);
    invd.resize(N5);

    d[0] = a[0];
    invd[0] = Real(1)/d[0];

    for(int i=1; i < N5; i++) { 
      d[i] = a[i] - (alpha[i-1]*alpha[i-1])/d[i-1];
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
      u[i] = alpha[i]/d[i];
    }

    off_diag_coeff.resize(N5);
    for(int i=0; i < N5; i++) { 
      off_diag_coeff[i] = -Real(0.5)*beta_tilde[i];
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
   *
   *  Flopcount: N5*6NcNs + (N5-2)*4NcNs = NcNs( 6N5 +4(N5-2) ) = (10N5-8) Nc Ns / cb_site
   */
  void 
  EvenOddPrecOvlapContFrac5DLinOpBaseArray::applyDiag(
    multi1d<LatticeFermion>& chi, 
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

    LatticeFermion tmp;  moveToFastMemoryHint(tmp);
    int G5=Ns*Ns-1;

    // First 0ne 
    // Operation: chi[0][rb[cb]] = a[0] G5 psi[0] + alpha[0]*psi[1]
    //
    //  Useful flops: 6Nc Ns / site 
    // tmp[rb[cb]] = Gamma(G5)*psi[0];
    // chi[0][rb[cb]] = a[0]*tmp;
    
    if( N5 > 1 ) { 
      chi[0][rb[cb]] = alpha[0]*psi[1] + a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }
    else {
      chi[0][rb[cb]] = a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }

    // All the rest
    for(int i=1; i < N5; i++) { 

      // Operation: 
      //   N5 - 1 times:
      //    chi[i]  = alpha[i-1]*psi[i-1] + a[i] Gamma_5 *psi[i]
      //   N5 - 2 times:
      //    chi[i] += alpha[i]*psi[i+1];
      //  Useful flops = (N5-1) * 6NcNs + (N5-2)*4Nc*Ns

      /*
      // B_{i-1} psi_[i-1]
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      tmp[rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] += a[i]*tmp;
      */
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1] + a[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-1) {
	chi[i][rb[cb]] += alpha[i]*psi[i+1];
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
   * Total flopcount: (N5-1)*4NcNs + 2NcNs + (N5-1)*6NcNs
   *                 = (N5-1)*10NcNs + 2NcNs
   *                 = (10N5-8) Nc Ns per_cb_site
   */
  void 
  EvenOddPrecOvlapContFrac5DLinOpBaseArray::applyDiagInv(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);

    multi1d<LatticeFermion> y(N5);   moveToFastMemoryHint(y);

    Real coeff;

    const int G5 = Ns*Ns-1;

    // Solve L y = psi
    y[0][rb[cb]] = psi[0];

    // (N5-1)*4NcNs
    for(int i = 1; i < N5; i++) { 
      // tmp[rb[cb]] = Gamma(G5)*y[i-1];
      // y[i][rb[cb]] = psi[i] - u[i-1]*tmp;
      y[i][rb[cb]] = psi[i] - u[i-1]*(GammaConst<Ns,Ns*Ns-1>()*y[i-1]);
    } 

    // Invert diagonal piece  y <- D^{-1} y
    // N5 times: y = (1/d_i) gamma_5 y[i] 
    // Useful flops: N5 * 2NcNs flops / site
    // Rolled this into the bottom loop

    // Backsubstitute U chi = y

    // 2NcNs

    chi[N5-1][rb[cb]] = invd[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*y[N5-1]);

    // N5-1 * 6NcNs
    for(int i = N5-2; i >= 0; i--) {
      // tmp[rb[cb]] = Gamma(G5)*chi[i+1]
      // chi[i][rb[cb]] = y[i] - u[i]*tmp;
      // y[i][rb[cb]] = invd[i]*(GammaConst<Ns,Ns*Ns-1>()*y[i]);
      // chi[i][rb[cb]] = y[i] - u[i]*(GammaConst<Ns,Ns*Ns-1>()*chi[i+1]);
      chi[i][rb[cb]] = GammaConst<Ns,Ns*Ns-1>()*(invd[i]*y[i]-u[i]*chi[i+1]); 
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
   * Total flopcount: (N5-1)*(Dslash_cb_flops + 2NcNs) per cb_site (isLastZeroP==true)
   *                   N5*(Dslash_cb_flops+2NcNs) per cb_site      (isLastZeroP==false)
   *
   * 
   */
  void EvenOddPrecOvlapContFrac5DLinOpBaseArray::applyOffDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi,
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    multi1d<LatticeFermion> tmp(N5);   moveToFastMemoryHint(tmp);

    int G5 = Ns*Ns-1;

    // Optimisation... do up to the penultimate block...

    // (N5-1)( Dslash + 2NcNs) flops/site
    Dslash.apply(tmp, psi, PLUS, cb);

    for(int i=0; i < N5-1; i++) { 
      /*
      // CB is CB of TARGET
      // gamma_5 Dslash is hermitian so I can ignore isign

      // Apply g5 Dslash
      Dslash.apply(tmp, psi[i], PLUS, cb);
  
      // chi[i][rb[cb]] = Gamma(G5)*tmp;
      */
      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[i][rb[cb]] = off_diag_coeff[i]*(GammaConst<Ns,Ns*Ns-1>()*tmp[i]);
    }

 
    // Only do last block if beta_tilde[i] is not zero
    // ( Dslash + 2NcNs) flops per site if done.
    if( !isLastZeroP ) {

      // Vector Dslash gets this appropriately
      //Dslash.apply(tmp, psi[N5-1], PLUS, cb);
      // chi[N5-1][rb[cb]] = Gamma(G5)*tmp;

      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[N5-1][rb[cb]] = off_diag_coeff[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*tmp[N5-1]);
    }
    else { 
      chi[N5-1][rb[cb]] = zero;
    }
  
    END_CODE();
  }


  // Derivative of even-odd linop component
  /* 
   * This is a copy of the above applyOffDiag with the D.apply(...) replaced
   * by  D.deriv(ds_tmp,...) like calls.
   */
  void 
  EvenOddPrecOvlapContFrac5DLinOpBaseArray::applyDerivOffDiag(
    multi1d<LatticeColorMatrix>& ds_u,
    const multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    int cb) const 
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
						   
    LatticeFermion tmp;             moveToFastMemoryHint(tmp);

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

    END_CODE();
  }

} // End Namespace Chroma
