// $Id: eoprec_ovlap_contfrac5d_pv_linop_array_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief Even-odd preconditioned Pauli-Villars Continued Fraction 5D
 */

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_pv_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  EvenOddPrecOvlapContFrac5DPVLinOpArray::EvenOddPrecOvlapContFrac5DPVLinOpArray(
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

    // It is always N5-1...
    Dslash.create(state,N5-1);

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
    multi1d<Real> d(N5);
    dinv.resize(N5);
    d[0] = a[0];
    dinv[0] = Real(1)/d[0];

    d[N5-1] = a[N5-1];
    dinv[N5-1] = Real(1)/d[N5-1];

    for(int i=1; i < N5-1; i++) { 
      d[i] = a[i] - (alpha[i-1]*alpha[i-1])/d[i-1];
      dinv[i] = Real(1)/d[i];
    }

    off_diag_coeff.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 
      off_diag_coeff[i] = -Real(0.5)*beta_tilde[i];
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
   *   6NcNs +  N5-2*[ 6NcNs ] + (N5-3)*4Nc Ns
   * = N5-1*(6NcNs) + (N5-3)*4NcNs
   * = 6*N5*NcNs - 6NcNs + 4N5NcNs - 12 NcNs
   *  =10*N5 NcNs - 18NcNs = (10*N5-18) NcNs
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

    LatticeFermion tmp;                moveToFastMemoryHint(tmp);
    int G5=Ns*Ns-1;

    // First 0ne 
    /*
    tmp[rb[cb]] = Gamma(G5)*psi[0];
    chi[0][rb[cb]] = a[0]*tmp;
    */
    if( N5 > 1 ) { 
      // 6NcNs flops
      chi[0][rb[cb]] = alpha[0]*psi[1] + a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }
    else {
      chi[0][rb[cb]] = a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }
    // All the rest

    // N5-2*[ 6NcNs ] + (N5-3)*4Nc Ns
    for(int i=1; i < N5-1; i++) 
    { 
      // B_{i-1} psi_[i-1]
      /*
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      tmp[rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] += a[i]*tmp;
      */
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1] + a[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

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
   *
   *  10NcNs(N5-2)+2NcNs
   * = (10N5 - 18)NcNs flops per cb_site
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

    multi1d<LatticeFermion> y(N5);      moveToFastMemoryHint(y);

    LatticeFermion tmp;                 moveToFastMemoryHint(tmp);
    Real coeff;

    const int G5 = Ns*Ns-1;

    // Solve L y = psi
    y[0][rb[cb]] = psi[0];

    // (N5-2) * 4NcNs flops
    for(int i = 1; i < N5-1; i++) 
    {
      /* tmp[rb[cb]] = Gamma(G5)*y[i-1]; */
      y[i][rb[cb]] = psi[i] - u[i-1]*(GammaConst<Ns,Ns*Ns-1>()*y[i-1]);
    } 

    // Backsubstitute U chi = y
 
    // Special note. chi[N5-1] = y[N5-1] = psi[N5-1] -- lowest
    // corner is a unit matrix
    chi[N5-1][rb[cb]] = psi[N5-1];

    
    // Real backsubstitutions starts from chi[N5-2]
    // 2NcNs flops
    chi[N5-2][rb[cb]] =  dinv[N5-2]*(GammaConst<Ns,Ns*Ns-1>()*y[N5-2]);

    // N5-2 * 6NcNs flops
    for(int i = N5-3; i >= 0; i--) {
      //tmp = Gamma(G5)*chi[i+1];
      // y[i][rb[cb]] = dinv[i]*(GammaConst<Ns,Ns*Ns-1>()*y[i]);
      chi[i][rb[cb]] = GammaConst<Ns,Ns*Ns-1>()*( dinv[i]*y[i] - u[i]*chi[i+1]);
    }

    END_CODE();
  }


  //! Apply the off diagonal block
  /*!
   *  (N5-1)*(Dslash + 2NcNs) 
   *
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

    multi1d<LatticeFermion> tmp(N5);      moveToFastMemoryHint(tmp);
    Real coeff;
    int G5 = Ns*Ns-1;

    // (N5-1)*(Dslash + 2NcNs)
    // Now donw with Dslash vector... 
    Dslash.apply(tmp, psi, PLUS, cb);

    for(int i=0; i < N5-1; i++) 
    {
      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[i][rb[cb]] = off_diag_coeff[i]*(GammaConst<Ns,Ns*Ns-1>()*tmp[i]);
    }
  
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
						   
    LatticeFermion tmp;                   moveToFastMemoryHint(tmp);
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
	// tmp[rb[cb]] = Gamma(G5)*chi[i];

	// Multiply coefficient
	// coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[cb]] = off_diag_coeff[i]*(GammaConst<Ns,Ns*Ns-1>()*chi[i]);

	// Apply g5 Dslash
	Dslash.deriv(ds_tmp, tmp, psi[i], PLUS, cb);
	ds_u += ds_tmp;
      }
      break;

    case MINUS:
      // Optimisation... do up to the penultimate block...
      for(int i=0; i < N5-1; i++) 
      {
	// CB is CB of TARGET
	// consider Dslash^dag gamma_5
	// tmp[rb[1-cb]] = Gamma(G5)*psi[i];

	// Multiply coefficient
	// coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[1-cb]] = off_diag_coeff[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

	// Apply g5 Dslash
	Dslash.deriv(ds_tmp, chi[i], tmp, MINUS, cb);
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


} // End Namespace Chroma


