// $Id: prec_dwf_linop_array_w.cc,v 1.15 2005-03-15 17:23:41 bjoo Exp $
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
   * \param aniso         aniso params  (Read)
   */
  EvenOddPrecDWLinOpArray::EvenOddPrecDWLinOpArray(
    const multi1d<LatticeColorMatrix>& u_, 
    const Real& WilsonMass_, const Real& m_q_, int N5_,
    const AnisoParam_t& aniso)
  {
    START_CODE();

    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1;
    N5  = N5_;

    multi1d<LatticeColorMatrix> u = u_;
    Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
  
    if (aniso.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != aniso.t_dir)
	  u[mu] *= ff;
      }
    }
    D.create(u);   // construct using possibly aniso glue

    InvTwoKappa = 1 + a5*(1 + (Nd-1)*ff - WilsonMass); 
    //InvTwoKappa =  WilsonMass - 5.0;
    TwoKappa = 1.0 / InvTwoKappa;
    Kappa = TwoKappa/2.0;
  
    invDfactor =1.0/(1.0  + m_q/pow(InvTwoKappa,N5));

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

    if( chi.size() != N5 ) chi.resize(N5);

    switch ( isign ) {
    
    case PLUS:
    {

      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_- * psi[s+1] - P_+ * psi[s-1]
	//chi[s][rb[cb]] = InvTwoKappa*psi[s] - 
	//   0.5*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]) ) ;

	// Recoded using chiralProject
	chi[s][rb[cb]] = InvTwoKappa*psi[s] - chiralProjectPlus(psi[s-1]);
	chi[s][rb[cb]] -= chiralProjectMinus(psi[s+1]);
      }

      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_- * psi[1] + mf* P_+ * psi[N5-1]
      // chi[0][rb[cb]] = InvTwoKappa*psi[0] - 
      //	0.5*( psi[1]   - m_q*psi[N5m1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]) ) ;

      // Recoded using chiralProject
      chi[0][rb[cb]] = InvTwoKappa*psi[0] - chiralProjectMinus(psi[1]);
      chi[0][rb[cb]] += m_q* chiralProjectPlus(psi[N5m1]);
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] +mf* P_- * psi[0]  -  P_+ * psi[N5-2]
      //chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] - 
      //	0.5*( psi[N5m2] - m_q *psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );

      // Recoded using chiralProject
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] - chiralProjectPlus(psi[N5m2]);
      chi[N5m1][rb[cb]] += m_q*chiralProjectMinus(psi[0]);
    }
    break ;

    case MINUS:
    {    


      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_+ * psi[s+1] - P_- * psi[s-1]
	//	chi[s][rb[cb]] = InvTwoKappa*psi[s] - 
	//  0.5*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s+1] - psi[s-1]) ) ;
	

	// Recoded using chiralProject
	chi[s][rb[cb]] = InvTwoKappa*psi[s] - chiralProjectPlus(psi[s+1]);
	chi[s][rb[cb]] -= chiralProjectMinus(psi[s-1]);

      }
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_+ * psi[1] + mf* P_- * psi[N5-1]
      // chi[0][rb[cb]] = InvTwoKappa*psi[0] - 
      //	0.5*( psi[1]   - m_q*psi[N5m1] + GammaConst<Ns,Ns*Ns-1>()*( psi[1]+m_q*psi[N5m1]) ) ;

      // Recoded using chiralProject
      chi[0][rb[cb]] = InvTwoKappa*psi[0] -  chiralProjectPlus(psi[1]);
      chi[0][rb[cb]] += m_q*chiralProjectMinus(psi[N5m1]);


      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] + mf* P_+ * psi[0]  -  P_- * psi[N5-2]
      // chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] - 
      // 0.5*( psi[N5m2] - m_q *psi[0] - GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );
      
      // Recoded using chiralProject
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1]-chiralProjectMinus(psi[N5m2]);
      chi[N5m1][rb[cb]] += m_q * chiralProjectPlus(psi[0]);
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

    if( chi.size() != N5 ) chi.resize(N5);

    switch ( isign ) {

    case PLUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      //Real fact(0.5*m_q*TwoKappa) ;
      // Recoded using ChiralProject which absorbs the factor of 0.5
      Real fact = m_q*TwoKappa;
      for(int s(0);s<N5-1;s++){
	//	chi[N5-1][rb[cb]] -= fact * (chi[s] - GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
       
	// Recoded using chiralProject 
	chi[N5-1][rb[cb]] -= fact*chiralProjectMinus(chi[s]);

	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++) {
	
	//chi[s][rb[cb]] += Kappa*(chi[s-1] + GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
	// Recoded using chiral project
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s-1]);

      }      

      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      
      for(int s(N5-2);s>-1;s--) {
	// chi[s][rb[cb]] += Kappa*(chi[s+1] - GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
	// Recoded using chiral Project
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s+1]);
      }

      //Finally the inverse of Rm 
      LatticeFermion tt;

      // fact = 0.5*m_q*TwoKappa;
      // Recoded using chiralProject which absorbs factor of 0.5
      fact = m_q*TwoKappa;

      // tt[rb[cb]] = fact*(chi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      tt[rb[cb]] = fact*chiralProjectPlus(chi[N5-1]);

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
      // Real fact(0.5*m_q*TwoKappa) ;
      // Recoded using chiralProject which absorbs factor of 0.5
      Real fact=m_q*TwoKappa;

      for(int s(0);s<N5-1;s++){
	// chi[N5-1][rb[cb]] -= fact * (chi[s] + GammaConst<Ns,Ns*Ns-1>()*chi[s])  ;
	// Recoded using chiralProject
	chi[N5-1][rb[cb]] -= fact*chiralProjectPlus(chi[s]);

	fact *= TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      
      for(int s(1);s<N5;s++) {
	// chi[s][rb[cb]] += Kappa*(chi[s-1] - GammaConst<Ns,Ns*Ns-1>()*chi[s-1]) ;
	// Recoded using chiralProject
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s-1]);

      }
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--) {
	//chi[s][rb[cb]] += Kappa*(chi[s+1] + GammaConst<Ns,Ns*Ns-1>()*chi[s+1]) ;
	// Recoded using chiralProject
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s+1]);;

      }

      //Finally the inverse of Rm 
      LatticeFermion tt;
      
      // Recoded using chiralProject
      tt[rb[cb]] = (m_q*TwoKappa)*chiralProjectMinus(chi[N5-1]);

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



  // THIS IS AN OPTIMIZED VERSION OF THE DERIVATIVE
  void 
  EvenOddPrecDWLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u,
				 const multi1d<LatticeFermion>& chi, 
				 const multi1d<LatticeFermion>& psi, 
				 enum PlusMinus isign) const
  {
    START_CODE();

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

    END_CODE();
  }

}; // End Namespace Chroma

