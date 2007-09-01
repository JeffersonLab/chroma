// $Id: eoprec_ovdwf_linop_array_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned Overlap-DWF (Borici) linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/eoprec_ovdwf_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   * \param N5_           extent of 5D  (Read)
   */
  void 
  EvenOddPrecOvDWLinOpArray::create(Handle< FermState<T,P,Q> > fs,
				    const Real& WilsonMass_, const Real& m_q_, int N5_)
  {
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1.0;
    N5  = N5_;


    // Extract a 4D fermstate to use for the Wilson dslash
    D.create(fs);

    c5InvTwoKappa = 1.0 - (Nd-WilsonMass) ;
    c5TwoKappa = 1.0 /  c5InvTwoKappa ;

    b5InvTwoKappa = 1.0 + (Nd-WilsonMass) ;
    b5TwoKappa = 1.0 /  b5InvTwoKappa ;

    TwoKappa = c5InvTwoKappa/b5InvTwoKappa ; 
    Kappa = TwoKappa/2.0 ;
  
    invDfactor =1.0/(1.0  + m_q*pow(TwoKappa,N5)) ;
  
  }


  //! Apply the even-even (odd-odd) coupling piece of the Borici operator
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
  EvenOddPrecOvDWLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
				       const multi1d<LatticeFermion>& psi, 
				       enum PlusMinus isign,
				       const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    switch ( isign ) {
    
    case PLUS:
    {
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_- * psi[s+1] + P_+ * psi[s-1]
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	  (0.5*c5InvTwoKappa)*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_- * psi[1] + mf* P_+ * psi[N5-1]
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
	(0.5*c5InvTwoKappa)*( psi[1]   - m_q*psi[N5m1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] +mf* P_- * psi[0]  -  P_+ * psi[N5-2]
      chi[N5m1][rb[cb]] = b5InvTwoKappa*psi[N5m1] - 
	(0.5*c5InvTwoKappa)*( psi[N5m2] - m_q *psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );
    }
    break ;

    case MINUS:
    {    
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] - P_+ * psi[s+1] - P_- * psi[s-1]
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	  (0.5*c5InvTwoKappa)*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s+1] - psi[s-1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] - P_+ * psi[1] + mf* P_- * psi[N5-1]
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
	(0.5*c5InvTwoKappa)*( psi[1]   - m_q*psi[N5m1] + GammaConst<Ns,Ns*Ns-1>()*( psi[1]+m_q*psi[N5m1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] + mf* P_+ * psi[0]  -  P_- * psi[N5-2]
      chi[N5m1][rb[cb]] = b5InvTwoKappa*psi[N5m1] - 
	(0.5*c5InvTwoKappa)*( psi[N5m2] - m_q *psi[0] - GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );
    }
    break ;
    }

    END_CODE();
  }


  //! Apply the inverse even-even (odd-odd) coupling piece of the Borici fermion operator
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
  EvenOddPrecOvDWLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
					  const multi1d<LatticeFermion>& psi, 
					  enum PlusMinus isign,
					  const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    // Copy and scale by TwoKappa (1/M0)
    for(int s(0);s<N5;s++)
      chi[s][rb[cb]] = b5TwoKappa * psi[s] ;


    switch ( isign ) {

    case PLUS:
    {
      
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
      LatticeFermion tt;            moveToFastMemoryHint(tt);
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
      LatticeFermion tt;                 moveToFastMemoryHint(tt);
      tt[rb[cb]] = (0.5*m_q*TwoKappa)*(chi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= tt  ;
	tt[rb[cb]] *= TwoKappa ;
      }
    }
    break ;
    }

    END_CODE();
  }


  //! Apply the even-odd (odd-even) coupling piece of the Borici operator
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
  EvenOddPrecOvDWLinOpArray::applyOffDiag(multi1d<LatticeFermion>& chi, 
					  const multi1d<LatticeFermion>& psi, 
					  enum PlusMinus isign,
					  const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);


    switch ( isign ) 
    {
    case PLUS:
    {
      LatticeFermion tmp;                        moveToFastMemoryHint(tmp);
      int otherCB = (cb + 1)%2 ;
      for(int s(1);s<N5-1;s++)
      {
	tmp[rb[otherCB]] = psi[s] + 
	  0.5*( psi[s+1] + psi[s-1] + 
		GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]));
	D.apply(chi[s],tmp,isign,cb);
	chi[s][rb[cb]] *= (-0.5);
      }
      
      int N5m1(N5-1);
      tmp[rb[otherCB]] = psi[0]  + 
	0.5*(psi[1] - m_q*psi[N5m1] 
	     - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]));
      D.apply(chi[0],tmp,isign,cb);
      chi[0][rb[cb]] *= (-0.5);
      
      int N5m2(N5-2);
      tmp[rb[otherCB]] = psi[N5m1] + 
	0.5*(psi[N5m2] - m_q *psi[0] 
	     + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]));

      D.apply(chi[N5m1],tmp,isign,cb);
      chi[N5m1][rb[cb]] *= (-0.5);
      
    }
    break ;

    case MINUS:
    {    
      multi1d<LatticeFermion> tmp(N5) ;     moveToFastMemoryHint(tmp);
      for(int s(0);s<N5;s++){
	D.apply(tmp[s],psi[s],isign,cb);
	tmp[s][rb[cb]] *= (-0.5) ;
      }
      for(int s(1);s<N5-1;s++){
	chi[s][rb[cb]] = tmp[s] + 
	  0.5*(tmp[s+1] + tmp[s-1] -
	       GammaConst<Ns,Ns*Ns-1>()*(tmp[s-1]-tmp[s+1]));
      }
      int N5m1(N5-1);
      chi[0][rb[cb]] = tmp[0]  + 
	0.5*(tmp[1] - m_q*tmp[N5m1] + 
	     GammaConst<Ns,Ns*Ns-1>()*(m_q*tmp[N5m1] + tmp[1]));
      
      
      int N5m2(N5-2);
      chi[N5m1][rb[cb]] = tmp[N5m1] + 
	0.5*(tmp[N5m2] - m_q *tmp[0] -
	     GammaConst<Ns,Ns*Ns-1>()*(tmp[N5m2] + m_q * tmp[0]));
      
    }
    break ;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }

  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  EvenOddPrecOvDWLinOpArray::Dminus(LatticeFermion& chi,
				    const LatticeFermion& psi,
				    enum PlusMinus isign,
				    int s5) const
  {
    LatticeFermion tt ;          moveToFastMemoryHint(tt);
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    chi = c5InvTwoKappa*psi +0.5*tt ;  //really -(-.5)D
  }
  
} // End Namespace Chroma

