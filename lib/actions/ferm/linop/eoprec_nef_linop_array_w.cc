// $Id: eoprec_nef_linop_array_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"

using namespace QDP::Hints;


namespace Chroma 
{ 

  //! Creation routine
  /*! \ingroup fermact
   *
   * \param fs            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param b5_           NEF parameter (Read)
   * \param c5_           NEF parameter (Read)
   * \param m_q_          quark mass    (Read)
   * \param N5_           extent of 5D  (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::create(Handle< FermState<T,P,Q> > fs,
				     const Real& WilsonMass_, const Real &b5_, 
				     const Real &c5_, const Real& m_q_, int N5_)
  {
    START_CODE();
  
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    b5  = b5_;
    c5  = c5_;
    N5  = N5_;
  
    D.create(fs, N5);
  
  
    c5InvTwoKappa = 1.0 - c5*(Nd-WilsonMass) ;
    c5TwoKappa = 1.0 / c5InvTwoKappa ;
  
    b5InvTwoKappa = 1.0 + b5*(Nd-WilsonMass) ;
    b5TwoKappa = 1.0 / b5InvTwoKappa ;
  
    //InvTwoKappa = b5InvTwoKappa/c5InvTwoKappa ; 
    TwoKappa =  c5InvTwoKappa/b5InvTwoKappa ;
    Kappa = TwoKappa/2.0 ;
  
    invDfactor =1.0/(1.0  + m_q*pow(TwoKappa,N5)) ;


    END_CODE();
  }


  //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice. 
   * Total flopcount: 6 *N5*Nc*Ns flops/site
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
					const multi1d<LatticeFermion>& psi, 
					enum PlusMinus isign,
					const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    // Real c5Fact(0.5*c5InvTwoKappa) ; // The 0.5 is for the P+ and P-

    Real c5InvTwoKappamf = m_q*c5InvTwoKappa;
    switch ( isign ) {
    
    case PLUS:
    {
      // Total flopcount:6 N5 Nc Ns flops/site
 
      // This loop (N5-2) * 6 Nc Ns flops/site
      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_- * psi[s+1] - P_+ * psi[s-1]


	//	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	//  c5Fact*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s-1] - psi[s+1]) ) ;

	// Recoded using chiralProject and BLASology
	// 4Nc Ns flops/site
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - c5InvTwoKappa*chiralProjectPlus(psi[s-1]);
	// 2Nc Ns flops/site
	chi[s][rb[cb]] -= c5InvTwoKappa*chiralProjectMinus(psi[s+1]);

      }



      //s=0 -- 1/2k psi[0] - P_- * psi[1] + mf* P_+ * psi[N5-1]
      //      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
      //	c5Fact*( psi[1]   - m_q*psi[N5m1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5m1] + psi[1]) ) ;

      // Recoded using chiralProject with BLAS-ology. c5Fact-s factor of 
      // 1/2 absorbed by projectors

      // These two lines 6Nc Ns flops
      // 4Nc Ns flops/site
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - c5InvTwoKappa*chiralProjectMinus(psi[1]);
      // 2Nc Ns flops/site
      chi[0][rb[cb]] += c5InvTwoKappamf*chiralProjectPlus(psi[N5-1]);


      //s=N5-1 -- 1/2k psi[N5-1] +mf* P_- * psi[0]  -  P_+ * psi[N5-2]
      // chi[N5m1][rb[cb]] = b5InvTwoKappa*psi[N5m1] - 
      //	c5Fact*( psi[N5m2] - m_q *psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5m2] + m_q * psi[0]) );

      // Recoded with chiral projector and BLAS ology
      // These two lines: 6Nc Ns flops /site
      // 4Nc Ns flops/site
      chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - c5InvTwoKappa*chiralProjectPlus(psi[N5-2]);

      // 2Nc Ns flops/site
      chi[N5-1][rb[cb]] += c5InvTwoKappamf*chiralProjectMinus(psi[0]);
    }
    break ;

    case MINUS:
    {    
      for(int s(1);s<N5-1;s++) { // 1/2k psi[s] - P_+ * psi[s+1] - P_- * psi[s-1]
	// chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - 
	//  c5Fact*( psi[s+1] + psi[s-1] + GammaConst<Ns,Ns*Ns-1>()*(psi[s+1] - psi[s-1]) ) ;

	// Recoded using chiral projector and BLASology
	chi[s][rb[cb]] = b5InvTwoKappa*psi[s] - c5InvTwoKappa*chiralProjectPlus(psi[s+1]);
	chi[s][rb[cb]] -= c5InvTwoKappa*chiralProjectMinus(psi[s-1]);
      }

      //s=0 -- 1/2k psi[0] - P_+ * psi[1] + mf* P_- * psi[N5-1]
      //      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - 
      //	c5Fact*( psi[1]   - m_q*psi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*( psi[1]+m_q*psi[N5-1]) ) ;
      chi[0][rb[cb]] = b5InvTwoKappa*psi[0] - c5InvTwoKappa*chiralProjectPlus(psi[1]);
      chi[0][rb[cb]] += c5InvTwoKappamf * chiralProjectMinus(psi[N5-1]);
      


      //s=N5-1 -- 1/2k psi[N5-1] + mf* P_+ * psi[0]  -  P_- * psi[N5-2]
      // chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - 
      //	c5Fact*( psi[N5-2] - m_q *psi[0] - GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q * psi[0]) );

      // Recoded using Chiral Projector and BLAS ology
      chi[N5-1][rb[cb]] = b5InvTwoKappa*psi[N5-1] - c5InvTwoKappa*chiralProjectMinus(psi[N5-2]);
      chi[N5-1][rb[cb]] += c5InvTwoKappamf*chiralProjectPlus(psi[0]);
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
   * Total flopcount: (10*N5 - 8)*Nc*Ns flops
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign,
					   const int cb) const
  {
    START_CODE();
 
    if( chi.size() != N5 ) chi.resize(N5);
   

    switch ( isign ) {

    case PLUS:
    {
      
      // Copying and scaling, application of Lm^{-1} and L^{-1} and D^{-1}
      // Coalesced into a fused loop.

      Real fact = m_q*TwoKappa*b5TwoKappa*invDfactor;
      Real invDTwoKappa = invDfactor*TwoKappa;
      Real invDb5TwoKappa = invDfactor*b5TwoKappa;

      // 2Nc Ns flops/site
      chi[0][rb[cb]] = b5TwoKappa*psi[0];
      
      // 4Nc Ns flops/site
      chi[N5-1][rb[cb]] = invDb5TwoKappa*psi[N5-1] 
	- fact * chiralProjectMinus(psi[0]);
      
      fact *= TwoKappa;

      // (N5-2)*6NcNs flops/site
      for(int s = 1; s < N5-1; s++) {
	
	// 2Nc Ns flops/site
	chi[s][rb[cb]] = b5TwoKappa * psi[s] ;
	// 2Nc Ns flops/site
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s-1]);
	// 2Nc Ns flops/site
	chi[N5-1][rb[cb]] -= fact * chiralProjectMinus(psi[s]);
	fact *= TwoKappa ;

      }
      
      // 2NcNs flops/site
      chi[N5-1][rb[cb]] += invDTwoKappa*chiralProjectPlus(chi[N5-2]);

      
      //The inverse of R. Back substitution...... Getting there! 
      // (N5-1)*2Nc*Ns flops /site total
      for(int s(N5-2);s>-1;s--)
	// 2Nc Ns2 flops/site
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s+1]);

      //Finally the inverse of Rm 
      fact = m_q*TwoKappa;

      // (N5-1)*2NcNs flops/site total
      for(int s(0);s<N5-1;s++){
	// 2Nc Ns flops
	chi[s][rb[cb]] -= fact*chiralProjectPlus(chi[N5-1]);
	fact *= TwoKappa;
      }
    }
    break ;
    
    case MINUS:
    {

      // Copy and scale by TwoKappa (1/M0)
      // First apply the inverse of Lm 
      // Now apply the inverse of L. Forward elimination
      // The inverse of D  now
      // All fused.

      Real fact = m_q*TwoKappa*b5TwoKappa*invDfactor ;
      Real invDTwoKappa = invDfactor*TwoKappa;
      Real invDb5TwoKappa = invDfactor*b5TwoKappa;

      chi[0][rb[cb]] = b5TwoKappa*psi[0];
      chi[N5-1][rb[cb]] = invDb5TwoKappa*psi[N5-1];
      chi[N5-1][rb[cb]] -= fact*chiralProjectPlus(psi[0]);
      fact *= TwoKappa;

      for(int s(1);s<N5-1;s++) {
	// 2Nc Ns flops/sit
	chi[s][rb[cb]] = b5TwoKappa * psi[s] ;
	chi[s][rb[cb]] += TwoKappa*chiralProjectMinus(chi[s-1]);
	chi[N5-1][rb[cb]] -= fact * chiralProjectPlus(psi[s]);
	fact *= TwoKappa ;
      }

      chi[N5-1][rb[cb]] += invDTwoKappa*chiralProjectMinus(chi[N5-2]);

      // That was easy....

      //The inverse of R. Back substitution...... Getting there! 
      for(int s(N5-2);s>-1;s--)
	chi[s][rb[cb]] += TwoKappa*chiralProjectPlus(chi[s+1]);
      
      //Finally the inverse of Rm 
      fact = m_q*TwoKappa;

      for(int s(0);s<N5-1;s++){
	chi[s][rb[cb]] -= fact*chiralProjectMinus(chi[N5-1]);
	fact *= TwoKappa ;
      }
    }
    break ;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }

  //! Apply the even-odd (odd-even) coupling piece of the NEF operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   * Total flopcount 6 N5 Nc Ns + N5 * checkerboarded Dslash flops per site
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecNEFDWLinOpArray::applyOffDiag(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign,
					   const int cb) const 
  {
    START_CODE();

    Real fb5 = -Real(0.5)*b5 ;

    // Recoding with chiral projectors, a former factor of 0.5 is absorbed
    // into the projector
    Real fc5 = -Real(0.5)*c5 ;
    Real fc5mf = fc5*m_q;

    if( chi.size() != N5 ) chi.resize(N5);
  
    switch ( isign ) 
    {
    case PLUS:
    {
      multi1d<LatticeFermion> tmp(N5); moveToFastMemoryHint(tmp);
      int otherCB = (cb + 1)%2 ;


      // (N5 -2 )*6 *Nc*Ns flops 
      for(int s = 1; s < N5-1; s++){
	//	tmp[s][rb[otherCB]] = fb5*psi[s] + 
	//   fc5*(psi[s+1] + psi[s-1] +
	//       GammaConst<Ns,Ns*Ns-1>()*(psi[s-1]-psi[s+1]));
	//
	// Recoded with chiral projectors and BLAS ology.
	// 4Nc Ns flops
	tmp[s][rb[otherCB]] = fb5*psi[s] + fc5*chiralProjectPlus(psi[s-1]);

	// 2Nc Ns flops
	tmp[s][rb[otherCB]]+= fc5*chiralProjectMinus(psi[s+1]);
      }
      
      
      // tmp[0][rb[otherCB]] = fb5*psi[0]  + 
      //	fc5*(psi[1] - m_q*psi[N5-1] 
      //     - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5-1] + psi[1]));
      //
      // Recoded with chiralProjec and BLAS ology
      // 6Nc Ns flops/site
      //
      // 4Nc Ns flops/site
      tmp[0][rb[otherCB]] = fb5*psi[0] +fc5*chiralProjectMinus(psi[1]);
      // 2Nc Ns flops/site
      tmp[0][rb[otherCB]] -= fc5mf*chiralProjectPlus(psi[N5-1]);

      // tmp[N5-1][rb[otherCB]] = fb5*psi[N5-1] + 
      //	fc5*( psi[N5-2] - m_q *psi[0] +
      //           GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q * psi[0]));

      // 6Nc Ns flops:
      //    4Nc Ns flops
      tmp[N5-1][rb[otherCB]] = fb5*psi[N5-1] + fc5*chiralProjectPlus(psi[N5-2]);
      //   2Nc Ns flops
      tmp[N5-1][rb[otherCB]] -= fc5mf*chiralProjectMinus(psi[0]);

      

      // Replace this with a vector Dslash in time -- done
      D.apply(chi,tmp, isign, cb);
      
    }
    break ;
    
    case MINUS:
    { 
      multi1d<LatticeFermion> tmp(N5) ; moveToFastMemoryHint(tmp);

      // Replace this with a vector Dslash in time -- done
      D.apply(tmp,psi,isign,cb);
      

      for(int s(1);s<N5-1;s++){
	//	chi[s][rb[cb]] = fb5*tmp[s] + 
	//  fc5*(tmp[s+1] + tmp[s-1] -
	//       GammaConst<Ns,Ns*Ns-1>()*(tmp[s-1]-tmp[s+1]));
	//  
	//  Recoded using chiralProject and BLAS ology
	chi[s][rb[cb]] = fb5*tmp[s] + fc5*chiralProjectPlus(tmp[s+1]);
	chi[s][rb[cb]] += fc5*chiralProjectMinus(tmp[s-1]);
      }

      // chi[0][rb[cb]] = fb5*tmp[0]  + 
      //	fc5*(tmp[1] - m_q*tmp[N5-1] + 
      // GammaConst<Ns,Ns*Ns-1>()*(m_q*tmp[N5-1] + tmp[1]));
      //
      // Recoded using chiralProject and BLAS ology
      chi[0][rb[cb]] = fb5*tmp[0] + fc5*chiralProjectPlus(tmp[1]);
      chi[0][rb[cb]] -= fc5mf*chiralProjectMinus(tmp[N5-1]);

      
      // chi[N5-1][rb[cb]] = fb5*tmp[N5-1] + 
      //	fc5*( tmp[N5-2] - m_q *tmp[0] -
      //          GammaConst<Ns,Ns*Ns-1>()*(tmp[N5-2] + m_q * tmp[0]));
      chi[N5-1][rb[cb]] = fb5*tmp[N5-1] + fc5*chiralProjectMinus(tmp[N5-2]);
      chi[N5-1][rb[cb]] -= fc5mf*chiralProjectPlus(tmp[0]);

    }
    break ;
    }

    //Done! That was not that bad after all....
    //See, I told you so...
    END_CODE();
  }



  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  EvenOddPrecNEFDWLinOpArray::Dminus(LatticeFermion& chi,
				     const LatticeFermion& psi,
				     enum PlusMinus isign,
				     int s5) const
  {
    LatticeFermion tt ;      moveToFastMemoryHint(tt);
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    chi = c5InvTwoKappa*psi + (0.5*c5)*tt ;//It really is -(-0.5*c5)D 
  }


} // End Namespace Chroma

