// $Id: prec_dwf_linop_array_w.cc,v 1.10 2005-01-03 04:10:48 edwards Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"

namespace Chroma 
{ 
  // Check Conventions... Currently I (Kostas) am using Blum et.al.


  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   * \param N5_           extent of 5D  (Read)
   */
  void 
  EvenOddPrecDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				  const Real& WilsonMass_, const Real& m_q_, int N5_)
  {
    START_CODE();

    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1.0;
    N5  = N5_;

    D.create(u_);

    InvTwoKappa = 5.0 - WilsonMass ; 
    //InvTwoKappa =  WilsonMass - 5.0;
    TwoKappa = 1.0 / InvTwoKappa ;
    Kappa = TwoKappa/2.0 ;
  
    invDfactor =1.0/(1.0  + m_q/pow(InvTwoKappa,N5)) ;

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
  EvenOddPrecDWLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign,
				     const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    switch ( isign ) {
    
    case PLUS:
    {
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_- * psi[s+1] + P_+ * psi[s-1]
	chi[s][rb[cb]] = InvTwoKappa*psi[s] - 
	  0.5*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_- * psi[1] + mf* P_+ * psi[N5-1]
      chi[0][rb[cb]] = InvTwoKappa*psi[0] - 
	0.5*( psi[1]   - m_q*psi[N5m1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] +mf* P_- * psi[0]  -  P_+ * psi[N5-2]
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] - 
	0.5*( psi[N5m2] - m_q *psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );
    }
    break ;

    case MINUS:
    {    
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] - P_+ * psi[s+1] - P_- * psi[s-1]
	chi[s][rb[cb]] = InvTwoKappa*psi[s] - 
	  0.5*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s+1] - psi[s-1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_+ * psi[1] + mf* P_- * psi[N5-1]
      chi[0][rb[cb]] = InvTwoKappa*psi[0] - 
	0.5*( psi[1]   - m_q*psi[N5m1] + GammaConst<Ns,Ns*Ns-1>()*( psi[1]+m_q*psi[N5m1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] + mf* P_+ * psi[0]  -  P_- * psi[N5-2]
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] - 
	0.5*( psi[N5m2] - m_q *psi[0] - GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );
    }
    break ;
    }

    END_CODE();
  }


  //! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi     Pseudofermion field              (Read)
   * \param isign   Flag ( PLUS | MINUS )            (Read)
   * \param cb      checkerboard ( 0 | 1 )           (Read)
   */
  void 
  EvenOddPrecDWLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
					const multi1d<LatticeFermion>& psi, 
					enum PlusMinus isign,
					const int cb) const
  {
    START_CODE();

    chi.resize(N5);

    switch ( isign ) {

    case PLUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_q*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * (chi[s] - GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += Kappa*(chi[s-1] + GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += Kappa*(chi[s+1] - GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      fact = 0.5*m_q*TwoKappa;
      tt[rb[cb]] = fact*(chi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break ;
    
    case MINUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_q*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] -= fact * (chi[s] + GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] += Kappa*(chi[s-1] - GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += Kappa*(chi[s+1] + GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      tt[rb[cb]] = (0.5*m_q*TwoKappa)*(chi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break ;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }


  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  EvenOddPrecDWLinOpArray::Dminus(LatticeFermion& chi,
				  const LatticeFermion& psi,
				  enum PlusMinus isign,
				  int s5) const
  {
    QDPIO::cerr << "Dminus not implemented" << endl;
    QDP_abort(1);
  }

  //! Apply the the even-odd block onto a source vector
  void 
  EvenOddPrecDWLinOpArray::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
					     enum PlusMinus isign) const
  {
    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    for(int s(0);s<N5;s++)
    {
      D.deriv(ds_tmp,chi[s],psi[s],isign,0);
      for(int mu(0);mu<Nd;mu++)
	ds_u[mu] += Real(-0.5)*ds_tmp[mu];
    }
  }
  

  //! Apply the the odd-even block onto a source vector
  void 
  EvenOddPrecDWLinOpArray::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
					     enum PlusMinus isign) const
  {
    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    for(int s(0);s<N5;s++)
      {
	D.deriv(ds_tmp,chi[s],psi[s],isign,1);
	for(int mu(0);mu<Nd;mu++)
	  ds_u[mu] += Real(-0.5)*ds_tmp[mu];
      }
  }




#if 0
  // NOT WORKING YET - USE SLOW DEFAULT

  // THIS IS AN OPTIMIZED VERSION OF THE DERIVATIVE
  /*! ds_u = -(1/4)*(1/(Nd+m))*[\dot(D)_oe*D_eo + D_oe*\dot(D_eo)]*psi */
  void 
  EvenOddPrecWilsonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
				const multi1d<LatticeFermion>& chi, 
				const multi1d<LatticeFermion>& psi, 
				enum PlusMinus isign) const
  {
    START_CODE();

    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

    ds_u.resize(Nd);

    // ds_u = -(1/4)*(1/(Nd+m))*[\dot(D)_oe*D_eo + D_oe*\dot(D_eo)]*psi

    // First term
    //   ds_u = chi_o^dag * \dot(D)_oe * (D_eo*psi_o)
    LatticeFermion tmp;
    D.apply(tmp, psi, isign, 0);
    D.deriv(ds_u, chi, tmp, isign, 1);

    // Second term
    //   ds_tmp = (D_eo^dag*chi_o)^dag * \dot(D)_eo * psi_o
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    D.apply(tmp, chi, msign, 0);
    D.deriv(ds_tmp, tmp, psi, isign, 0);

    for(int mu = 0; mu < Nd; ++mu)
    {
      ds_u[mu] += ds_tmp[mu];
      ds_u[mu] *= Real(-0.25)*invfact;
    }

    END_CODE();
  }
#endif

}; // End Namespace Chroma

