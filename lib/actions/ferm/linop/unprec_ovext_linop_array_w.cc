/* $Id: unprec_ovext_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
*  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
*/

#include "chromabase.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   */
  void 
  UnprecOvExtLinOpArray::create(Handle< FermState<T,P,Q> > state,
				const int Npoles_,
				const Real& coeffP_,
				const multi1d<Real>& resP_,
				const multi1d<Real>& rootQ_,
				const multi1d<Real>& beta_,
				const Real& OverMass_,
				const Real& Mass_,
				const Real& b5_,
				const Real& c5_)
  {
    Npoles = Npoles_;
    N5 = 2*Npoles_ + 1;
    R = (Real(1) + Mass_)/(Real(1)-Mass_);
    alpha = b5_ + c5_;
    a5 = b5_ - c5_;

    Dw.create(state, -OverMass_);
    fbc = state->getFermBC();
    coeffP = coeffP_;

    p_by_beta_sqrt.resize(Npoles);
    q_sqrt.resize(Npoles);
    beta.resize(Npoles);

    for(int i=0; i < Npoles; i++) { 
      beta[i] = beta_[i];
      q_sqrt[i] = sqrt(rootQ_[i]*beta[i]);
      p_by_beta_sqrt[i] = sqrt(resP_[i]*beta[i]);
    }
      

  }


  //! Apply the operator onto a source vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecOvExtLinOpArray::operator() (multi1d<LatticeFermion>& chi,
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);
    int G5 = Ns*Ns - 1;

    LatticeFermion tmp4; moveToFastMemoryHint(tmp4);

    Real q_sqrt_afive, two_q_sqrt, ftmp, Two;
    Two = Real(2);
	
    switch (isign) { 
    case PLUS: 
      {
	multi1d<LatticeFermion> tmp5(N5);
	moveToFastMemoryHint(tmp5);

	
	// Get a vector of wilson dirac op applications
	for(int i=0; i < N5; i++) { 
	  Dw(tmp5[i], psi[i], PLUS);
	}

	// Start off Chi[N5-1] 
	//         = [ 2 R gamma_5 + (R a5 + p0 alpha )gamma_5 Dw ] psi[N5-1]

	// 2 R gamma_5 psi_[N5-1]
	chi[N5-1] = Gamma(G5)*psi[N5-1];
	ftmp = Two*R;
	chi[N5-1] *= ftmp;

	// gamma_5 Dw psi[N5-1]
	tmp4 = Gamma(G5)*tmp5[N5-1];

	// (R a5 + p0) gamma_5 Dw psi[N5-1]
	ftmp=(R*a5 + coeffP*alpha);
	chi[N5-1] += ftmp*tmp4;
	


	// Loop through the length of the 5th D in steps of 2
	// Keep count of the pole we are on with p
	int p=0;
	for(int i=0; i < N5-2; p++, i+=2) {

	  q_sqrt_afive = q_sqrt[p]*a5;
	  two_q_sqrt = Two*q_sqrt[p];

	  // first row of block 	  
	  // chi[i] = -beta_p alpha H_w psi_i + sqrt(q_p)(2+a5 D_w) psi_i+1
	  //        = 2 sqrt(q_p) psi_{i+1} 
	  //        + a5 sqrt(q_p) Dw psi{i+1}
	  //        - beta_{i} alpha gamma_5 Dw psi{i}

	  // First this part: 
	  //        = 2 sqrt(q_p) psi_{i+1} 
	  //        + a5 sqrt(q_p) Dw psi{i+1}
	  chi[i] = two_q_sqrt*psi[i+1] + q_sqrt_afive * tmp5[i+1];

	  // Now  - beta_{p} alpha gamma_5 Dw psi{i}
	  ftmp = alpha;
	  tmp4 = Gamma(G5)*tmp5[i];
	  chi[i] -= ftmp*tmp4;


	  // Second row of block:	  

	  // chi[i+1] = sqrt(q_p)(2 + a5 D) psi_i 
	  //          + (alpha/beta) g5 D psi_{i+1}
	  //          + sqrt(p_p/beta_p)(2 + a5 Dw) psi_{N5-1}
	  //
	  //   = 2 sqrt(q_p)  psi_i + 2 sqrt(p_p/beta_p)psi_{N5-1}
	  //   + a5 sqrt(q_p) Dpsi_i + a5 sqrt(p_p/beta_p) Dpsi_{N5-1}
	  //   + alpha/beta_p gamma_5 Dpsi_{i+1}
	  
	  // sqrt(p_p/beta_p) precomputed in p_by_beta_sqrt

	  // First: 2 sqrt(q_p)  psi_i + 2 sqrt(p_p/beta_p)psi_{N5-1}
	  ftmp = Two*p_by_beta_sqrt[p];
	  chi[i+1] = two_q_sqrt * psi[i] + ftmp*psi[N5-1];

	  // Second: a5 sqrt(q_p) Dpsi_i + a5 sqrt(p_p/beta_p) Dpsi_{N5-1}
	  chi[i+1] += q_sqrt_afive*tmp5[i];
	  ftmp = a5*p_by_beta_sqrt[p];
	  chi[i+1] += ftmp*tmp5[N5-1];

	  // Third:  alpha/beta_p gamma_5 Dpsi_{i+1}
	  ftmp  = alpha*beta[p];
	  tmp4 = Gamma(G5)*tmp5[i+1];
	  chi[i+1] += ftmp*tmp4;


	  // Update last row: -sqrt(p_p/beta_p)*(2+a5Dw)psi[i+1]
	  //        =      -2 sqrt(p_p/beta_p) psi_[i+1]
	  //              -a5 sqrt(p_p/beta_p) Dw psi[i+1];
	  ftmp = Two*p_by_beta_sqrt[p];
	  chi[N5-1] -= ftmp*psi[i+1];
	  ftmp = a5*p_by_beta_sqrt[p];
	  chi[N5-1] -= ftmp*tmp5[i+1];

	}
      }
      break;

    case MINUS:
      {
	multi1d<LatticeFermion> tmp5_1(N5);
	multi1d<LatticeFermion> tmp5_2(N5);

	moveToFastMemoryHint(tmp5_1);
	moveToFastMemoryHint(tmp5_2);

	// Will need to multiply by D^{dagger]. Pull it outside.
	// We collect 2 5D vectors. One to be multiplied by Dw^{dag}
	// and one not.
	
	// We start with the N5-1th term.
	// this is:
	//
	// [ R(2 + a5 Dw^dag)g5 + p0 alpha Hw = R(2 + a5D^dag)g5 + p0 alpha Dw^dag g5 ] psi[N5-1]
	// = 2Rg5 psi[N5-1] + Dw^{dag}( R a5 + p0 alpha ) g5 psi[N5-1];
	// tmp5_1[N5-1] = 2Rg5 psi[N5-1]
	// tmp5_2[N5-1] = (R a5 + p0 alpha) g5 psi[N5-1]
	tmp5_1[N5-1] = Gamma(G5)*psi[N5-1];
	tmp5_2[N5-1] = tmp5_1[N5-1];
	ftmp = Two*R;
	tmp5_1[N5-1] *= ftmp;
	ftmp = (R*a5 + coeffP*alpha);
	tmp5_2[N5-1] *= ftmp;

	int p = 0;
	for(int i=0; i < N5-1; i+=2, p++) { 
	  
	  // first row: (chi[i] = tmp5_1_i + D^dag tmp5_2_i)
	  // chi[i] = -beta_p alpha Hw psi_i
	  //          + sqrt(q_p)*(2 + a5 D^dag) psi_i+1
	  //
	  //       =  -beta_p alpha Ddag g_5 psi_i
	  //         + 2 sqrt(q_p) psi_{i+1} 
	  //         + D^dag ( sqrt(q_p)*a5 * psi_{i+1} )

	  //   tmp5_1 [i] = 2 sqrt(q_p) psi_{i+1}
	  //   tmp5_2 [i] = -beta_p alpha g_5 psi_i 
	  //               + a5 sqrt(q_p) psi[i+1]
	  //
	  two_q_sqrt = Two * q_sqrt[p];
	  q_sqrt_afive = a5* q_sqrt[p];

	  tmp5_1[i] = two_q_sqrt*psi[i+1];

	  tmp4 = Gamma(G5)*psi[i];
	  ftmp = alpha;
	  tmp5_2[i] = q_sqrt_afive*psi[i+1] - ftmp*tmp4;
	 

	  // second row of block:
	  //
	  // chi[i+1] = sqrt(q_p) ( 2 + a5 D^dag) psi_i
	  //          + alpha/beta_p H psi_{i+1} 
	  //          - sqrt(p_p/beta_p)( 2 + a5 D^dag) psi_N5-1
	  //
	  // = 2 sqrt(q_p) psi_i - 2 sqrt(p_p/beta_p)psi_N5-1
	  // + D^dag{ a5 sqrt(q_p) psi_i - a5 sqrt(p_p/beta_p) psi_N5-1
	  //         + alpha/beta_p gamma_5 psi_{i+1}
	  //
	  // tmp5_1[i+1] = 2 sqrt(q_p) psi_i - 2 sqrt(p_p/beta_p) psi[N5-1]
	  //
	  // tmp5_2[i+1] = a5 sqrt(q_p) psi_i - a5 sqrt(p_p/beta_p) psi_N5-1
	  //              + alpha/beta_p gamma_5 psi_{i+1}
	  
	  ftmp = Two*p_by_beta_sqrt[p];
	  tmp5_1[i+1] = two_q_sqrt * psi[i] - ftmp * psi[N5-1];
	  
	  ftmp = a5*p_by_beta_sqrt[p];
	  tmp5_2[i+1] = q_sqrt_afive*psi[i] - ftmp*psi[N5-1];
	  tmp4 = Gamma(G5)*psi[i+1];
	  ftmp = alpha*beta[p];
	  tmp5_2[i+1] += ftmp*tmp4;

	  // Coupling term: sqrt(p_p/beta_p)(2 + a5 Ddag) psi[i+1]
	  //
	  //    = 2 sqrt(p_p/beta_p) psi[i+1]
	  //    + Ddag  ( a5 sqrt(p_p/beta_p) ) psi_{i+1}
	  //
	  // tmp5_1[N5-1] = 2 sqrt(p_p/beta_p) psi[i+1]
	  // tmp5_2[N5-1] = a5 sqrt(p_p/beta_p) psi_{i+1}
	  ftmp = Two*p_by_beta_sqrt[p];
	  tmp5_1[N5-1] += ftmp*psi[i+1];

	  ftmp = a5*p_by_beta_sqrt[p];
	  tmp5_2[N5-1] += ftmp *psi[i+1];
	}
	
	// Vector Dirac op on tmp5_2
	for(int i=0; i < N5; i++) { 
	  Dw(chi[i], tmp5_2[i], MINUS);
	  chi[i]+=tmp5_1[i];
	}
	
      }
      break;
    default: 
      QDPIO::cerr << "Unknown value for PLUS /MINUS: " << isign << endl;
      QDP_abort(1);
    };

    getFermBC().modifyF(chi);

    END_CODE();
  }


  //! Derivative
  void 
  UnprecOvExtLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign) const
  {
    START_CODE();
    QDPIO::cout << "Not yet implemented " << endl;
    QDP_abort(1);

    getFermBC().zero(ds_u);

    END_CODE();
  }
}; // End Namespace Chroma

