// $Id: eoprec_nef_general_linop_array_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/eoprec_nef_general_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param b5_           NEF parameter array (Read)
   * \param c5_           NEF parameter array (Read)
   * \param m_q_          quark mass    (Read)
   * \param N5_           extent of 5D  (Read)
   */
  void 
  EvenOddPrecGenNEFDWLinOpArray::create(Handle< FermState<T,P,Q> > fs,
					const Real& WilsonMass_, 
					const Real& m_q_,
					const multi1d<Real>& b5_, 
					const multi1d<Real>& c5_,
					int N5_)
  {
    START_CODE();
  
    WilsonMass = WilsonMass_;
    m_q = m_q_;

    // Sanity checking:
    if( b5_.size() != N5_ ) { 
      QDPIO::cerr << "b5 array size and N5 are inconsistent" << endl;
      QDPIO::cerr << "b5_array.size() = " << b5_.size() << endl;
      QDPIO::cerr << "N5_ = " << N5_ << endl << flush;
      QDP_abort(1);
    }
    N5  = N5_;
    D.create(fs, N5);

    b5.resize(N5);
    c5.resize(N5);
    for(int i=0; i < N5; i++) { 
      b5[i] = b5_[i];
      c5[i] = c5_[i];
    }

    f_plus.resize(N5);
    f_minus.resize(N5);
    for(int i=0; i < N5; i++) { 
      f_plus[i] = b5[i]*( Real(Nd) - WilsonMass ) + 1;
      f_minus[i]= c5[i]*( Real(Nd) - WilsonMass ) - 1;
    }


    l.resize(N5-1);
    r.resize(N5-1);
  
    l[0] = -m_q*f_minus[N5-1]/f_plus[0];
    r[0] = -m_q*f_minus[0]/f_plus[0];

    for(int i=1; i < N5-1; i++) { 
      l[i] = -(f_minus[i-1]/f_plus[i])*l[i-1];
      r[i] = -(f_minus[i]/f_plus[i])*r[i-1];
    }

    a.resize(N5-1);
    b.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 
      a[i] = f_minus[i+1]/f_plus[i];
      b[i] = f_minus[i]/f_plus[i];
    }

    d.resize(N5);
    for(int i=0; i < N5; i++) { 
      d[i] = f_plus[i];
    }
    // Last bits of d can be computed 2 ways which should be equal
    Real tmp1 = f_minus[N5-2]*l[N5-2];
    Real tmp2 = f_minus[N5-1]*r[N5-2];

    // Sanity checking:
    // QDPIO::cout << "  The following 2 should be equal: ";
    // QDPIO::cout << tmp1 << " and " << tmp2 << endl;

  
    // Subtrace ONLY ONE of them onto d[N5-1]
    d[N5-1] -=  tmp1;
 
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
  EvenOddPrecGenNEFDWLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
					   const multi1d<LatticeFermion>& psi, 
					   enum PlusMinus isign,
					   int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    switch ( isign ) {
    
    case PLUS:
    {
      Real fact;


      // fplus[0]*psi[0] + fminus[0]P_- psi[1] - m fminus[0] P_+ psi[N5-1]
      // = fplus[0]*psi[0] + (fminus/2)[ (1-g5)psi[1] - m(1+g5) psi[N5-1] ]
      // = fplus[0]*psi[0] + (fminus/2)[ psi[1] - m psi[N5-1] 
      //                                  - g5( psi[1] + m psi[N5-1] ) ]
      // = fplus[0]*psi[0] + (fminus/2)[ tmp_1 - g5 tmp2 ]
      //
      // with tmp1 = psi[1] - m psi[N5-1]
      //      tmp2 = psi[1] + m psi[N5-1]

      /* Old Code */
      /*
      LatticeFermion tmp1, tmp2;

      tmp1[rb[cb]] = psi[1] - m_q*psi[N5-1];
      tmp2[rb[cb]] = psi[1] + m_q*psi[N5-1];
      fact = Real(0.5)*f_minus[0];
      
      chi[0][rb[cb]] = f_plus[0]*psi[0] + 
	fact*( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );
      */

      // Recoded using chiral projectors
      fact = f_minus[0]*m_q;
      chi[0][rb[cb]] = f_plus[0]*psi[0];
      chi[0][rb[cb]] += f_minus[0]*chiralProjectMinus(psi[1]);
      chi[0][rb[cb]] -= fact*chiralProjectPlus(psi[N5-1]);


      // -m fminus[N5-1] P_- psi[0] + f_minus[N5-1] P_+ psi[N5-2] 
      //     + f_plus[N5-1] psi[N5-1]
      //
      // fminus[N5-1] ( -m P_- psi[0] + P_+ psi[N5-2] ) 
      //     + f_plus[N5-1] psi[N5-1]
      //
      // fminus[N5-1]/2 ( (1 + g5) psi[N5-2] - m(1-g5) psi[0] )
      //     + f_plus[N5-1] psi[N5-1]
      //
      // fminus[N5-1]/2 ( psi[N5-2] - m psi[0] + g5( psi[N5-2] +m psi[0] ) )
      //     + f_plus[N5-1] psi[N5-1]
      //
      // fminus[N5-1]/2 (  tmp1  + g5 ( tmp2 )) + f_plus[N5-1] psi[N5-1]
      //
      //  tmp1 = psi[N5-2] - m psi[0]
      //  tmp2 = psi[N5-2] + m psi[0]

      /* OLD CODE */
      /*
      fact = Real(0.5)*f_minus[N5-1];
      tmp1[rb[cb]] = psi[N5-2] - m_q * psi[0];
      tmp2[rb[cb]] = psi[N5-2] + m_q * psi[0];
      
      chi[N5-1][rb[cb]] = f_plus[N5-1]*psi[N5-1] 
	+ fact*( tmp1 + GammaConst<Ns,Ns*Ns-1>() * tmp2);
      */

      // Recoded using chiral projectors
      // -m fminus[N5-1] P_- psi[0] + f_minus[N5-1] P_+ psi[N5-2] 
      //     + f_plus[N5-1] psi[N5-1]
      //
      // fminus[N5-1] ( -m P_- psi[0] + P_+ psi[N5-2] ) 
      //     + f_plus[N5-1] psi[N5-1]
      fact = f_minus[N5-1]*m_q;
      chi[N5-1][rb[cb]] = f_plus[N5-1]*psi[N5-1];
      chi[N5-1][rb[cb]] += f_minus[N5-1]*chiralProjectPlus(psi[N5-2]);
      chi[N5-1][rb[cb]] -= fact*chiralProjectMinus(psi[0]);


      for(int s=1; s<N5-1; s++) {


	// fminus[s] P_+ psi[s-1] + fplus[s] psi[s] + fminus[s] P_- psi[s+1]
	//= fplus[s] psi[s] + fminus[s]( P_+ psi[s-1] + P_-psi[s+1])
        //= fplus[s] psi[s] + (fminus[s]/2)( psi[s-1] + psi[s+1] + 
	//                                   g_5*{ psi[s-1] - psi[s + 1] } )

	// OLD CODE 
        /*
	fact = Real(0.5)*f_minus[s];
	tmp1[rb[cb]] = psi[s-1] + psi[s+1];
	tmp2[rb[cb]] = psi[s-1] - psi[s+1];
	
	chi[s][rb[cb]] = f_plus[s]*psi[s] +
	  fact*( tmp1 + GammaConst<Ns,Ns*Ns-1>() * tmp2  );
	*/

	// Recoded using chiral projectors
	chi[s][rb[cb]] = f_plus[s]*psi[s];
	chi[s][rb[cb]] += f_minus[s]*chiralProjectMinus(psi[s+1]);
	chi[s][rb[cb]] += f_minus[s]*chiralProjectPlus(psi[s-1]);       

      }
      

    }
    break ;

    case MINUS:
    {    

      // Scarily different. Daggering in the 5th dimension
      // changes the f_minuses from constant along a row to
      // varying along a row...

      Real fact;

      // OLD CODE 
      /*
      LatticeFermion tmp1, tmp2;

      // Fact is constant now.
      fact = Real(0.5);

      
      tmp1[rb[cb]] = f_minus[1]*psi[1] - m_q*f_minus[N5-1]*psi[N5-1];
      tmp2[rb[cb]] = f_minus[1]*psi[1] + m_q*f_minus[N5-1]*psi[N5-1];

      chi[0][rb[cb]] = f_plus[0]*psi[0] + 
	fact*( tmp1 + GammaConst<Ns, Ns*Ns-1>()*tmp2 );
      */

      // Recoded using chiral projectors
      fact = m_q * f_minus[N5-1];
      chi[0][rb[cb]] = f_plus[0]*psi[0];
      chi[0][rb[cb]] += f_minus[1]*chiralProjectPlus(psi[1]);
      chi[0][rb[cb]] -= fact*chiralProjectMinus(psi[N5-1]);


      // OLD CODE 
      /*
      tmp1[rb[cb]] = f_minus[N5-2]*psi[N5-2] - m_q * f_minus[0]* psi[0];
      tmp2[rb[cb]] = f_minus[N5-2]*psi[N5-2] + m_q * f_minus[0]* psi[0];
      
      chi[N5-1][rb[cb]] = f_plus[N5-1]*psi[N5-1]
	+ fact*( tmp1 - GammaConst<Ns,Ns*Ns-1>()*tmp2 ) ;
      */

      // Recoded using Chiral Projectors
      fact = m_q * f_minus[0];
      chi[N5-1][rb[cb]]  = f_plus[N5-1]*psi[N5-1];
      chi[N5-1][rb[cb]] += f_minus[N5-2]*chiralProjectMinus(psi[N5-2]);
      chi[N5-1][rb[cb]] -= fact*chiralProjectPlus(psi[0]);


      for(int s=1; s<N5-1; s++) {

	// OLD CODE
	/*
	tmp1[rb[cb]] = f_minus[s-1]*psi[s-1] + f_minus[s+1]*psi[s+1];
	tmp2[rb[cb]] = f_minus[s-1]*psi[s-1] - f_minus[s+1]*psi[s+1];

	chi[s][rb[cb]] = f_plus[s]*psi[s] +
	  fact*( tmp1 - GammaConst<Ns,Ns*Ns-1>()*tmp2 );
	*/

	// Recoded using Chiral Projectors
	chi[s][rb[cb]] = f_plus[s]*psi[s];
	chi[s][rb[cb]] += f_minus[s-1]*chiralProjectMinus(psi[s-1]);
	chi[s][rb[cb]] += f_minus[s+1]*chiralProjectPlus(psi[s+1]);

      }

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
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)
   */
  void 
  EvenOddPrecGenNEFDWLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
					      const multi1d<LatticeFermion>& psi, 
					      enum PlusMinus isign,
					      int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    // I use two temporaries
    multi1d<LatticeFermion> z(N5);        moveToFastMemoryHint(z);
    multi1d<LatticeFermion> z_prime(N5);  moveToFastMemoryHint(z_prime);
    Real fact;

    switch ( isign ) {

    case PLUS:
    {

      // First apply the inverse of Lm :
      // Solve Lm z = chi
      z[N5-1][rb[cb]] = psi[N5-1];
      for(int s=0; s < N5-1; s++){
	z[s][rb[cb]] = psi[s];


	// OLD CODE
	// The factor of 1/2 is for the projection expression
	//	fact = Real(0.5)*l[s];
        //      z[N5-1][rb[cb]] -=  fact*(psi[s] - GammaConst<Ns,Ns*Ns-1>()*psi[s])  ;
       
	// Recoded with projector
	z[N5-1][rb[cb]] -= l[s]*chiralProjectMinus(psi[s]);

      }
      
      //Now apply the inverse of L. Forward elimination 
      //
      // L z' = z
      z_prime[0][rb[cb]] = z[0];
      for(int s = 0; s < N5-1; s++) {

	// OLD CODE 
	// The factor of 1/2 is for the projection part
	//	fact = Real(0.5)*a[s];
	// z_prime[s+1][rb[cb]] = z[s+1] - fact*(z_prime[s] + GammaConst<Ns,Ns*Ns-1>()*z_prime[s]);

	// Recoded using chiralProjectors
	z_prime[s+1][rb[cb]] = z[s+1] - a[s]*chiralProjectPlus(z_prime[s]);
		     
      }
	
      // z = D^{-1} z'
      for(int s=0; s < N5; s++) { 
	fact = Real(1)/d[s];
	z[s][rb[cb]] = fact*z_prime[s];
      }

      // z' = U^{-1} z => U z' = z
      //The inverse of R. Back substitution...... Getting there! 
      z_prime[N5-1][rb[cb]] = z[N5-1];
      for(int s=N5-2; s >=0; s-- ) { 

	// OLD CODE
	// fact = Real(0.5)*b[s];

	// z_prime[s][rb[cb]] = z[s] - fact*(z_prime[s+1] - GammaConst<Ns,Ns*Ns-1>()*z_prime[s+1]);

	// Recoded Using ChiralProjectors
	z_prime[s][rb[cb]] = z[s] - b[s]*chiralProjectMinus(z_prime[s+1]);
      }

      //Finally the inverse of Rm 
      chi[N5-1][rb[cb]] = z_prime[N5-1];
      for(int s=0; s < N5-1; s++) { 

	// OLD CODE
	// fact = Real(0.5)*r[s]; 

	// chi[s][rb[cb]] = z_prime[s] - fact*(chi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
	
	// Recoded using ChiralProjectors
	chi[s][rb[cb]] = z_prime[s] - r[s]*chiralProjectPlus(chi[N5-1]);

      }
    }
    break ;
    
    case MINUS:
    {

      // First apply the inverse of Rm :
      // Solve (Rm)^{T} z = psi
      
      // This is the same as the Lm case above but with l P_ => r P+
      
      z[N5-1][rb[cb]] = psi[N5-1];
      for(int s=0; s < N5-1; s++){
	z[s][rb[cb]] = psi[s];

	// OLD CODE 
	// The factor of 1/2 is for the projection expression
	// fact = Real(0.5)*r[s];
        // z[N5-1][rb[cb]] -=  fact*(psi[s] + GammaConst<Ns,Ns*Ns-1>()*psi[s])  ;
	// Recoded Using Chiral Projectors
	z[N5-1][rb[cb]] -= r[s]*chiralProjectPlus(psi[s]);
      }
      
      //Now apply the inverse of U^{T}. Forward elimination 
      //
      // U^T z' = z . This is the same as L z' = z above but with 
      // a P+ => b P-

      z_prime[0][rb[cb]] = z[0];
      for(int s = 0; s < N5-1; s++) {

	// OLD CODE
	// The factor of 1/2 is for the projection part
	// fact = Real(0.5)*b[s];
	// z_prime[s+1][rb[cb]] = z[s+1] - fact*(z_prime[s] - GammaConst<Ns,Ns*Ns-1>()*z_prime[s]);

	// Recoded using Chiral Projectors
	z_prime[s+1][rb[cb]] = z[s+1] - b[s]*chiralProjectMinus(z_prime[s]);
      }
	
      // z = D^{-1} z'
      for(int s=0; s < N5; s++) { 
	fact = Real(1)/d[s];
	z[s][rb[cb]] = fact*z_prime[s];
      }

      // z' = (L^T)^{-1} z ie L^T z' = z
      //
      // This is the same as the U z' =z case above
      // but with b P_ => a P+

      z_prime[N5-1][rb[cb]] = z[N5-1];
      for(int s=N5-2; s >=0; s-- ) { 

	// OLD CODE
	// fact = Real(0.5)*a[s];
	//z_prime[s][rb[cb]] = z[s] - fact*(z_prime[s+1] + GammaConst<Ns,Ns*Ns-1>()*z_prime[s+1]);

	// Recoded using Chiral Projectors
	z_prime[s][rb[cb]] = z[s] - a[s]*chiralProjectPlus(z_prime[s+1]);
      }

      //Finally the inverse of Lm^T 
      //Same as the inverse of R above but with r P+ => l P-
      chi[N5-1][rb[cb]] = z_prime[N5-1];
      for(int s=0; s < N5-1; s++) { 

	// OLD CODE
	//	fact = Real(0.5)*l[s]; 
	// chi[s][rb[cb]] = z_prime[s] - fact*(chi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);

	
	// Recoded using Chiral Projectors
	chi[s][rb[cb]] = z_prime[s] - l[s]*chiralProjectMinus(chi[N5-1]);
      }
     
    }
    break ;
    }

    END_CODE();
  }

  //! Apply the even-odd (odd-even) coupling piece of the NEF operator
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
  EvenOddPrecGenNEFDWLinOpArray::applyOffDiag(multi1d<LatticeFermion>& chi, 
					      const multi1d<LatticeFermion>& psi, 
					      enum PlusMinus isign,
					      int cb) const 
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    switch ( isign ) 
    {
    case PLUS:
    {
      multi1d<LatticeFermion> tmp(N5);  moveToFastMemoryHint(tmp);
      Real fact1;
      Real fact2;
      Real fact2mf;

      int otherCB = (cb + 1)%2 ;

      //  OLD CODE 
      /*
      LatticeFermion tmp1;
      LatticeFermion tmp2;


      fact1 = -Real(0.5)*b5[0];
      fact2 = -Real(0.25)*c5[0];
      tmp1[rb[otherCB]] = psi[1] - m_q*psi[N5-1];
      tmp2[rb[otherCB]] = psi[1] + m_q*psi[N5-1];
      
      tmp[0][rb[otherCB]] = fact1*psi[0] 
	+ fact2*(tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2);
      */

      // Recoded using Chiral Projectors
      fact1 = -Real(0.5)*b5[0];
      fact2 = -Real(0.5)*c5[0];
      fact2mf = m_q*fact2;

      tmp[0][rb[otherCB]] = fact1*psi[0];
      tmp[0][rb[otherCB]] += fact2*chiralProjectMinus(psi[1]);
      tmp[0][rb[otherCB]] -= fact2mf*chiralProjectPlus(psi[N5-1]);


      // OLD CODE 
      /*
      fact1 = -Real(0.5)*b5[N5-1];
      fact2 = -Real(0.25)*c5[N5-1];

      tmp1[rb[otherCB]]=psi[N5-2]-m_q*psi[0];
      tmp2[rb[otherCB]]=psi[N5-2]+m_q*psi[0];

      tmp[N5-1][rb[otherCB]] = fact1*psi[N5-1]
	+ fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2 );
      */

      // Recoded using Chiral Projectors
      fact1 = -Real(0.5)*b5[N5-1];
      fact2 = -Real(0.5)*c5[N5-1];
      fact2mf = m_q*fact2;

      tmp[N5-1][rb[otherCB]] = fact1*psi[N5-1];
      tmp[N5-1][rb[otherCB]] += fact2*chiralProjectPlus(psi[N5-2]);
      tmp[N5-1][rb[otherCB]] -= fact2mf*chiralProjectMinus(psi[0]);

      for(int s=1; s < N5-1; s++) { 

	// OLD CODE 
	/*
	  fact1 = -Real(0.5)*b5[s];
	  fact2 = -Real(0.25)*c5[s];
	  
	  tmp1[rb[otherCB]]=psi[s-1] + psi[s+1];
	  tmp2[rb[otherCB]]=psi[s-1] - psi[s+1];
	  
	  tmp[s][rb[otherCB]]= fact1*psi[s] 
	  + fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2) ;
	*/

	// Recoded using Chiral Projectors
	fact1 = -Real(0.5)*b5[s];
	fact2 = -Real(0.5)*c5[s];
	tmp[s][rb[otherCB]] = fact1*psi[s];
	tmp[s][rb[otherCB]] += fact2*chiralProjectPlus(psi[s-1]);
	tmp[s][rb[otherCB]] += fact2*chiralProjectMinus(psi[s+1]);
      }

      // Replace with vector dslash later -- done
      D.apply(chi, tmp, isign, cb);
      

    }
    break ;
    
    case MINUS:
    { 
      multi1d<LatticeFermion> tmp_d(N5) ; moveToFastMemoryHint(tmp_d);

      // Replace with a vector dslash later -- done
      D.apply(tmp_d, psi, isign, cb);
      

      // OLD CODE
      /*
      LatticeFermion tmp1;
      LatticeFermion tmp2;


      Real one_quarter = Real(0.25);

      Real factb= -Real(0.5)*b5[0];
      Real fact  = m_q*c5[N5-1];
      tmp1[rb[cb]] = c5[1]*tmp_d[1] - fact*tmp_d[N5-1];
      tmp2[rb[cb]] = c5[1]*tmp_d[1] + fact*tmp_d[N5-1];

      chi[0][rb[cb]] = factb*tmp_d[0]
	- one_quarter*( tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2 );
      */

      // Recoded with Chiral Projector
      Real factb = -Real(0.5)*b5[0];
      Real factc1 = -Real(0.5)*c5[1];
      Real factc2 = -Real(0.5)*m_q*c5[N5-1];

      chi[0][rb[cb]] = factb*tmp_d[0];
      chi[0][rb[cb]] += factc1*chiralProjectPlus(tmp_d[1]);
      chi[0][rb[cb]] -= factc2*chiralProjectMinus(tmp_d[N5-1]);


      // OLD CODE 
      /*
      factb = -Real(0.5)*b5[N5-1];
      fact = m_q * c5[0];
      
      tmp1[rb[cb]] = c5[N5-2]*tmp_d[N5-2] - fact * tmp_d[0];
      tmp2[rb[cb]] = c5[N5-2]*tmp_d[N5-2] + fact * tmp_d[0];

      chi[N5-1][rb[cb]] = factb * tmp_d[N5-1]
	- one_quarter * ( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );
      */

      // Recoded with Chiral Projector
      factb = -Real(0.5)*b5[N5-1];
      factc1 = -Real(0.5)*c5[N5-2];
      factc2 = -Real(0.5)*m_q*c5[0];

      chi[N5-1][rb[cb]] = factb*tmp_d[N5-1];
      chi[N5-1][rb[cb]] += factc1*chiralProjectMinus(tmp_d[N5-2]);
      chi[N5-1][rb[cb]] -= factc2*chiralProjectPlus(tmp_d[0]);


      for(int s=1; s<N5-1; s++){
	factb = -Real(0.5) * b5[s];
	factc1 = -Real(0.5) * c5[s-1];
	factc2 = -Real(0.5) * c5[s+1];

	// OLD CODE 
	//	tmp1[rb[cb]] = c5[s-1]*tmp_d[s-1] + c5[s+1]*tmp_d[s+1];
	//     tmp2[rb[cb]] = c5[s-1]*tmp_d[s-1] - c5[s+1]*tmp_d[s+1];

	// chi[s][rb[cb]] = factb*tmp_d[s] 
	//  - one_quarter*( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );

	// Recoded with Chiral Projector
	chi[s][rb[cb]] = factb*tmp_d[s];
	chi[s][rb[cb]] += factc1*chiralProjectMinus(tmp_d[s-1]);
	chi[s][rb[cb]] += factc2*chiralProjectPlus(tmp_d[s+1]);

      }
    }
    break ;
    }


    END_CODE();
  }


  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  EvenOddPrecGenNEFDWLinOpArray::Dminus(LatticeFermion& chi,
					const LatticeFermion& psi,
					enum PlusMinus isign,
					int s5) const
  {
    LatticeFermion tt ;         moveToFastMemoryHint(tt);
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    Real fact = Real(0.5)*c5[s5];
    chi = f_minus[s5]*psi - fact*tt ;
  }



  // Derivative of even-odd linop component
  /* 
   * This is a copy of the above applyOffDiag with the D.apply(...) replaced
   * by  D.deriv(ds_tmp,...) like calls.
   */
  void 
  EvenOddPrecGenNEFDWLinOpArray::applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u,
						   const multi1d<LatticeFermion>& chi, 
						   const multi1d<LatticeFermion>& psi, 
						   enum PlusMinus isign,
						   int cb) const 
  {
    START_CODE();

    ds_u.resize(Nd);
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
						   
    ds_u = zero;

    switch ( isign ) {
    case PLUS:
      {
	LatticeFermion tmp;            moveToFastMemoryHint(tmp);
	Real fact1;
	Real fact2;
	Real fact2mf;
	int otherCB = (cb + 1)%2 ;

	// OLD CODE
	/*
	LatticeFermion tmp1;
	LatticeFermion tmp2;

	
	fact1 = -Real(0.5)*b5[0];
	fact2 = -Real(0.25)*c5[0];
	tmp1[rb[otherCB]] = psi[1] - m_q*psi[N5-1];
	tmp2[rb[otherCB]] = psi[1] + m_q*psi[N5-1];
	
	tmp[rb[otherCB]] = fact1*psi[0] 
	  + fact2*(tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2);
	*/

	// Recoded with chiralProject
	fact1 = -Real(0.5)*b5[0];
	fact2 = -Real(0.5)*c5[0];
	fact2mf = m_q*fact2;

	tmp[rb[otherCB]] = fact1*psi[0];
	tmp[rb[otherCB]] += fact2*chiralProjectMinus(psi[1]);
	tmp[rb[otherCB]] -= fact2mf*chiralProjectPlus(psi[N5-1]);

	D.deriv(ds_u, chi[0], tmp, isign, cb);
	
	/*
	fact1 = -Real(0.5)*b5[N5-1];
	fact2 = -Real(0.25)*c5[N5-1];
	
	tmp1[rb[otherCB]]=psi[N5-2]-m_q*psi[0];
	tmp2[rb[otherCB]]=psi[N5-2]+m_q*psi[0];
	
	tmp[rb[otherCB]] = fact1*psi[N5-1]
	  + fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2 );
	*/

	// Recoded using Chiral Projectors
	fact1 = -Real(0.5)*b5[N5-1];
	fact2 = -Real(0.5)*c5[N5-1];
	fact2mf = m_q*fact2;
	
	tmp[rb[otherCB]] = fact1*psi[N5-1];
	tmp[rb[otherCB]] += fact2*chiralProjectPlus(psi[N5-2]);
	tmp[rb[otherCB]] -= fact2mf*chiralProjectMinus(psi[0]);

	D.deriv(ds_tmp, chi[N5-1], tmp, isign, cb);
	ds_u += ds_tmp;
	
	
	for(int s=1; s < N5-1; s++) { 

	  // OLD CODE
	  /*
	  fact1 = -Real(0.5)*b5[s];
	  fact2 = -Real(0.25)*c5[s];
	  
	  tmp1[rb[otherCB]]=psi[s-1] + psi[s+1];
	  tmp2[rb[otherCB]]=psi[s-1] - psi[s+1];
	  
	  tmp[rb[otherCB]]= fact1*psi[s] 
	    + fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2) ;
	  */

	  	// Recoded using Chiral Projectors
	  fact1 = -Real(0.5)*b5[s];
	  fact2 = -Real(0.5)*c5[s];
	  tmp[rb[otherCB]] = fact1*psi[s];
	  tmp[rb[otherCB]] += fact2*chiralProjectPlus(psi[s-1]);
	  tmp[rb[otherCB]] += fact2*chiralProjectMinus(psi[s+1]);
	  D.deriv(ds_tmp, chi[s], tmp, isign, cb);
	  ds_u += ds_tmp;
	}
      }
      break ;
      
    case MINUS:
      { 
	
	// New Minus case by Balint
	// There are 3 cases which are summed over.
	// Best way to compute them is to go back to basics.
	
	// s=0:
	// 
	//  F = -(b_5[0]/2) chi(0,cb)^dag Ddot(dagger, cb) psi(0, 1-cb)
	//      -(c_5[1]/2) chi(0,cb)^dag P_{+} Ddot(dagger, cb) psi(1, 1-cb)
	// +( m c_5[N-1]/2) chi(0,cb)^dag P_{-} Ddot(dagger, cb) psi(N-1, 1-cb)
	
	// General 0 < s < N-1:
	//
	//  F = -(b_5[s]/2) chi(s,cb)^dag Ddot(dagger, cb) psi(s, 1-cb) 
	//    -(c_5[s-1]/2) chi(s,cb)^dag P_{-} Ddot(dagger, cb) psi(s-1, 1-cb)
	//    -(c_5[s+1]/2) chi(s,cb)^dag P_{+} Ddot(dagger, cb) psi(s+1, 1-cb)
	
	// and for s = N-1
	//
	// F = -(b_5[N-1]/2) chi(N-1,cb)^dag Ddot(dagger,cb) psi(N-1, 1-cb)
	//     -(c_5[N-2]/2) chi(N-1,cb)^dag P_{-} Ddot(dagger,cb) psi(N-2, 1-cb)
	//   +( m c_5[0] /2) chi(N-1,cb)^dag P_{+} Ddot(dagger,cb) psi(0, 1-cb)
	//
	//
	// If we compute X, P_{+}X, and P_{-}X up front then we need to 
	// evaluate 3 force terms for each i. If we multiply out the projections	// and do gamma_5's etc it brings the number of terms up to 6.
	// 
	// We also fold the coefficients into the right hand (chi) vectors
	// as needed

	// Storage for P_{+} X
	multi1d<LatticeFermion> chi_plus(N5);  moveToFastMemoryHint(chi_plus);
	multi1d<LatticeFermion> chi_minus(N5); moveToFastMemoryHint(chi_minus);
	
	for(int s=0; s < N5; s++) {
	  // OLD CODE
	  /*
	  chi_plus[s][rb[cb]]  = chi[s] + GammaConst<Ns,Ns*Ns-1>()*chi[s];
	  chi_plus[s][rb[cb]] *= Real(0.5);
	  
	  chi_minus[s][rb[cb]] = chi[s] - GammaConst<Ns,Ns*Ns-1>()*chi[s];
	  chi_minus[s][rb[cb]]*= Real(0.5);
	  */

	  // Recoded using chiralProject
	  chi_plus[s][rb[cb]] = chiralProjectPlus(chi[s]);
	  chi_minus[s][rb[cb]] = chiralProjectMinus(chi[s]);

	  // Zero out unwanted checkerboard... May remove this later
	  chi_plus[s][rb[1-cb]] = zero;
	  chi_minus[s][rb[1-cb]] = zero;
	}
	
	
	ds_u = zero;
	LatticeFermion chi_tmp = zero;       moveToFastMemoryHint(chi_tmp);


	Real ftmp;

	// First term:
	// -(b5[0]/2) chi[0]^dag Ddot^dag psi[0]
	ftmp=Real(0.5)*b5[0];
	chi_tmp[rb[cb]] = ftmp*chi[0];
	D.deriv(ds_tmp, chi_tmp, psi[0], isign, cb);
	ds_u -= ds_tmp;
	

	// -(c5[1]/2) chi[0]^dag P+  Ddot^dag psi[1]
	ftmp=Real(0.5)*c5[1];
	chi_tmp[rb[cb]] = ftmp*chi_plus[0];
	D.deriv(ds_tmp, chi_tmp, psi[1], isign, cb);
	ds_u -= ds_tmp;
	
	// +(m c5[N-1]) chi[0]^dag P_- Ddot^dag psi[N-1]
	ftmp=Real(0.5)*c5[N5-1]*m_q;
	chi_tmp[rb[cb]] = ftmp*chi_minus[0];
	D.deriv(ds_tmp, chi_tmp, psi[N5-1], isign, cb);
	ds_u += ds_tmp;
	
	
	// Middle bits
	for(int s=1; s < N5-1; s++) { 
	
	  //    -(b5[s]/2) chi[s] Ddot^dag psi[s]	
	  ftmp = Real(0.5)*b5[s];
	  chi_tmp[rb[cb]] = ftmp*chi[s];
	  D.deriv(ds_tmp, chi_tmp, psi[s], isign, cb);
	  ds_u -= ds_tmp;
	  
	  //    -(c5[s-1]/2) chi[s]^dag P- Ddot^dag psi[s-1]
	  ftmp = Real(0.5)*c5[s-1];
	  chi_tmp[rb[cb]] = ftmp*chi_minus[s];
	  D.deriv(ds_tmp, chi_tmp, psi[s-1], isign,cb);
	  ds_u -= ds_tmp;
	  
	  //   -(c5[s+1]/2) chi[s]^dag P+ Ddot^dag psi[s+1]
	  ftmp = Real(0.5)*c5[s+1];
	  chi_tmp[rb[cb]] = ftmp*chi_plus[s];
	  D.deriv(ds_tmp, chi_tmp, psi[s+1], isign, cb);
	  ds_u -= ds_tmp;
	  
	}
	
	
	// Lat bit
	//  - (b5[N5-1]/2) chi[N5-1] Ddot^dagger psi[N5-1]
	ftmp = Real(0.5)*b5[N5-1];
	chi_tmp[rb[cb]] = ftmp * chi[N5-1];
	D.deriv(ds_tmp, chi_tmp, psi[N5-1], isign, cb);
	ds_u -= ds_tmp;
	
	//  -(c5[N5-2]/2) chi[N5-1] P_ Ddot^dagger psi[N5-2]
	ftmp = Real(0.5)*c5[N5-2];
	chi_tmp[rb[cb]] = ftmp*chi_minus[N5-1];
	D.deriv(ds_tmp, chi_tmp, psi[N5-2], isign, cb);
	ds_u -= ds_tmp;
	
	//  +(m c5[0]/2) chi[N5-1] P+ Ddot^dagger psi[0]
	ftmp = Real(0.5)*c5[0]*m_q;
	chi_tmp[rb[cb]] = ftmp * chi_plus[N5-1];
	D.deriv(ds_tmp, chi_tmp, psi[0], isign, cb);
	ds_u += ds_tmp;
	
      }
      break ;
    }


    END_CODE();
  }

} // End Namespace Chroma

