// -*- C++ -*-
// $Id: eo4d_dwf_factors_array_w.cc,v 1.2 2003-11-23 06:00:50 edwards Exp $
/*! \file
 *  \brief 4D Even Odd preconditioned domain-wall fermion linear operator
 *   factors 
 */

#error "OBSOLETE: DO NOT USE"


// Check Convensions... Currently I am using Blum et.al.

#include "chromabase.h"
//#include "actions/ferm/linop/eo4d_dwf_factors_array_w.h"
#include "actions/ferm/linop/dslash_w.h"

void QDPEvenOdd4dDWDslashArray::create(const multi1d<LatticeColorMatrix>& u_, 
				    const int N5_)
{
  a5  = 1.0;
  N5  = N5_;
  D.create(u_);
}

#ifdef I_IMPLEMENTED_IT_IN_THE_HEADER
//! Apply the even-odd coupling piece  of the  domain-wall fermion operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void QDPEvenOdd4dDWDslashArray::operator()(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign) const
{
  for(int s(0);s<N5;s++){
    D(chi[s],psi[s],isign);
    chi[s] *= 0.5 ;
  }
}
#endif


void QDPEvenOdd4dDWDiagArray::create(const Real& WMass, const Real& mf, 
				     const int N5_)
{
  a5  = 1.0;
  N5  = N5_;
  WilsonMass = WMass ;
  InvTwoKappa = 5.0 - WilsonMass ;
  TwoKappa = 1.0 / InvTwoKappa ;
  Kappa = TwoKappa/2.0 ;
  m_f = mf ;
  
  invDfactor =1.0/(1.0 + m_f/pow(InvTwoKappa,N5)) ;

}

//! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void QDPEvenOdd4dDWDiagArray::operator()(multi1d<LatticeFermion>& chi, 
					 const multi1d<LatticeFermion>& psi, 
					 enum PlusMinus isign) const
{
  switch ( isign ) {

  case PLUS:
    {
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_- * psi[s+1] + P_+ * psi[s-1]
	chi[s] = InvTwoKappa*psi[s] + 
	  0.5*( psi[s+1] + psi[s-1] + Gamma(15)*(psi[s-1] - psi[s+1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] + P_- * psi[1] - mf* P_+ * psi[N5-1]
      chi[0] = InvTwoKappa*psi[0] + 
	0.5*( psi[1]   - m_f*psi[N5m1] - Gamma(15)*(m_f*psi[N5m1] + psi[1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] -mf* P_- * psi[0]  +  P_+ * psi[N5-2]
      chi[N5m1] = InvTwoKappa*psi[N5m1] + 
	0.5*( psi[N5m2] - m_f *psi[0] + Gamma(15)*(psi[N5m2] + m_f * psi[0]) );
    }
    break ;
    
  case MINUS:
    {
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_+ * psi[s+1] + P_- * psi[s-1]
	chi[s] = InvTwoKappa*psi[s] + 
	  0.5*( psi[s+1] + psi[s-1] + Gamma(15)*(psi[s+1] - psi[s-1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] + P_+ * psi[1] - mf* P_- * psi[N5-1]
      chi[0] = InvTwoKappa*psi[0] + 
	0.5*( psi[1]   - m_f*psi[N5m1] + Gamma(15)*( psi[1]+m_f*psi[N5m1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] -mf* P_+ * psi[0]  +  P_- * psi[N5-2]
      chi[N5m1] = InvTwoKappa*psi[N5m1] + 
	0.5*( psi[N5m2] - m_f *psi[0] - Gamma(15)*(psi[N5m2] + m_f * psi[0]) );
    }
    break ;
  }
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
void QDPEvenOdd4dDWDiagArray::apply(multi1d<LatticeFermion>& chi, 
				    const multi1d<LatticeFermion>& psi, 
				    enum PlusMinus isign,
				    const int cb) const
{
  switch ( isign ) {
    
  case PLUS:
    {
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_- * psi[s+1] + P_+ * psi[s-1]
	chi[s][rb[cb]] = InvTwoKappa*psi[s] + 
	  0.5*( psi[s+1] + psi[s-1] + Gamma(15)*(psi[s-1] - psi[s+1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] + P_- * psi[1] - mf* P_+ * psi[N5-1]
      chi[0][rb[cb]] = InvTwoKappa*psi[0] + 
	0.5*( psi[1]   - m_f*psi[N5m1] - Gamma(15)*(m_f*psi[N5m1] + psi[1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] -mf* P_- * psi[0]  +  P_+ * psi[N5-2]
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] + 
	0.5*( psi[N5m2] - m_f *psi[0] + Gamma(15)*(psi[N5m2] + m_f * psi[0]) );
    }
    break ;

  case MINUS:
    {    
      for(int s(1);s<N5-1;s++) // 1/2k psi[s] + P_+ * psi[s+1] + P_- * psi[s-1]
	chi[s][rb[cb]] = InvTwoKappa*psi[s] + 
	  0.5*( psi[s+1] + psi[s-1] + Gamma(15)*(psi[s+1] - psi[s-1]) ) ;
      
      int N5m1(N5-1) ;
      //s=0 -- 1/2k psi[0] + P_+ * psi[1] - mf* P_- * psi[N5-1]
      chi[0][rb[cb]] = InvTwoKappa*psi[0] + 
	0.5*( psi[1]   - m_f*psi[N5m1] + Gamma(15)*( psi[1]+m_f*psi[N5m1]) ) ;
      
      int N5m2(N5-2);
      //s=N5-1 -- 1/2k psi[N5-1] -mf* P_+ * psi[0]  +  P_- * psi[N5-2]
      chi[N5m1][rb[cb]] = InvTwoKappa*psi[N5m1] + 
	0.5*( psi[N5m2] - m_f *psi[0] - Gamma(15)*(psi[N5m2] + m_f * psi[0]) );
    }
    break ;
  }
}

//! Apply the inverse even-even (odd-odd) coupling piece of the domain-wall fermion operator
/*!
 * \ingroup linop
 *
 * The operator acts on the entire lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 * \param cb      checkerboard ( 0 | 1 )               (Read)
 */
void QDPEvenOdd4dDWDiagArray::applyInv(multi1d<LatticeFermion>& chi, 
				       const multi1d<LatticeFermion>& psi, 
				       enum PlusMinus isign,
				       const int cb) const
{

  switch ( isign ) {

  case PLUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_f*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] += fact * (chi[s] - Gamma(15)*chi[s])  ;
	fact *= -TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] -= Kappa*(chi[s-1] + Gamma(15)*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] -= Kappa*(chi[s+1] - Gamma(15)*chi[s+1]) ;
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      tt[rb[cb]] = (0.5*m_f*TwoKappa)*(chi[N5-1] + Gamma(15)*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] += tt  ;
	tt[rb[cb]] *= -TwoKappa ;
      }
    }
    break ;
    
  case MINUS:
    {
      // Copy and scale by TwoKappa (1/M0)
      for(int s(0);s<N5;s++)
	chi[s][rb[cb]] = TwoKappa * psi[s] ;
      
      // First apply the inverse of Lm 
      Real fact(0.5*m_f*TwoKappa) ;
      for(int s(0);s<N5-1;s++){
	chi[N5-1][rb[cb]] += fact * (chi[s] + Gamma(15)*chi[s])  ;
	fact *= -TwoKappa ;
      }
      
      //Now apply the inverse of L. Forward elimination 
      for(int s(1);s<N5;s++)
	chi[s][rb[cb]] -= Kappa*(chi[s-1] - Gamma(15)*chi[s-1]) ;
      
      //The inverse of D  now
      chi[N5-1][rb[cb]] *= invDfactor ;
      // That was easy....
      
      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] -= Kappa*(chi[s+1] + Gamma(15)*chi[s+1]) ;
      
      //Finally the inverse of Rm 
      LatticeFermion tt;
      tt[rb[cb]] = (0.5*m_f*TwoKappa)*(chi[N5-1] - Gamma(15)*chi[N5-1]);
      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] += tt  ;
	tt[rb[cb]] *= -TwoKappa ;
      }
    }
    break ;
  }

  //Done! That was not that bad after all....
  
}

