// $Id: prec_nef_general_linop_array_w.cc,v 1.3 2004-12-12 21:22:16 edwards Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_nef_general_linop_array_w.h"


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
EvenOddPrecGenNEFDWLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
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
  D.create(u_);

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
				      const int cb) const
{
  START_CODE();

  chi.resize(N5);

  switch ( isign ) {
    
  case PLUS:
    {
      Real fact;
      LatticeFermion tmp1, tmp2;


      // fplus[0]*psi[0] + fminus[0]P_- psi[1] - m fminus[0] P_+ psi[N5-1]
      // = fplus[0]*psi[0] + (fminus/2)[ (1-g5)psi[1] - m(1+g5) psi[N5-1] ]
      // = fplus[0]*psi[0] + (fminus/2)[ psi[1] - m psi[N5-1] 
      //                                  - g5( psi[1] + m psi[N5-1] ) ]
      // = fplus[0]*psi[0] + (fminus/2)[ tmp_1 - g5 tmp2 ]
      //
      // with tmp1 = psi[1] - m psi[N5-1]
      //      tmp2 = psi[1] + m psi[N5-1]

      tmp1[rb[cb]] = psi[1] - m_q*psi[N5-1];
      tmp2[rb[cb]] = psi[1] + m_q*psi[N5-1];
      fact = Real(0.5)*f_minus[0];
      
      chi[0][rb[cb]] = f_plus[0]*psi[0] + 
	fact*( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );


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
      
      fact = Real(0.5)*f_minus[N5-1];
      tmp1[rb[cb]] = psi[N5-2] - m_q * psi[0];
      tmp2[rb[cb]] = psi[N5-2] + m_q * psi[0];
      
      chi[N5-1][rb[cb]] = f_plus[N5-1]*psi[N5-1] 
	+ fact*( tmp1 + GammaConst<Ns,Ns*Ns-1>() * tmp2);
      
      for(int s=1; s<N5-1; s++) {
	
	// fminus[s] P_+ psi[s-1] + fplus[s] psi[s] + fminus[s] P_- psi[s+1]
	//= fplus[s] psi[s] + fminus[s]( P_+ psi[s-1] + P_-psi[s+1])
        //= fplus[s] psi[s] + (fminus[s]/2)( psi[s-1] + psi[s+1] + 
	//                                   g_5*{ psi[s-1] - psi[s + 1] } )
	fact = Real(0.5)*f_minus[s];
	tmp1[rb[cb]] = psi[s-1] + psi[s+1];
	tmp2[rb[cb]] = psi[s-1] - psi[s+1];
	
	chi[s][rb[cb]] = f_plus[s]*psi[s] +
	  fact*( tmp1 + GammaConst<Ns,Ns*Ns-1>() * tmp2  );
      }
      

    }
    break ;

  case MINUS:
    {    

      // Scarily different. Daggering in the 5th dimension
      // changes the f_minuses from constant along a row to
      // varying along a row...

      Real fact;
      LatticeFermion tmp1, tmp2;

      // Fact is constant now.
      fact = Real(0.5);

      
      tmp1[rb[cb]] = f_minus[1]*psi[1] - m_q*f_minus[N5-1]*psi[N5-1];
      tmp2[rb[cb]] = f_minus[1]*psi[1] + m_q*f_minus[N5-1]*psi[N5-1];

      chi[0][rb[cb]] = f_plus[0]*psi[0] + 
	     fact*( tmp1 + GammaConst<Ns, Ns*Ns-1>()*tmp2 );


      tmp1[rb[cb]] = f_minus[N5-2]*psi[N5-2] - m_q * f_minus[0]* psi[0];
      tmp2[rb[cb]] = f_minus[N5-2]*psi[N5-2] + m_q * f_minus[0]* psi[0];
      
      chi[N5-1][rb[cb]] = f_plus[N5-1]*psi[N5-1]
	+ fact*( tmp1 - GammaConst<Ns,Ns*Ns-1>()*tmp2 ) ;
	

      for(int s=1; s<N5-1; s++) {

	tmp1[rb[cb]] = f_minus[s-1]*psi[s-1] + f_minus[s+1]*psi[s+1];
	tmp2[rb[cb]] = f_minus[s-1]*psi[s-1] - f_minus[s+1]*psi[s+1];

	chi[s][rb[cb]] = f_plus[s]*psi[s] +
	  fact*( tmp1 - GammaConst<Ns,Ns*Ns-1>()*tmp2 );
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
				      const int cb) const
{
  START_CODE();

  chi.resize(N5);

  // I use two temporaries
  multi1d<LatticeFermion> z(N5);
  multi1d<LatticeFermion> z_prime(N5);
  Real fact;

  switch ( isign ) {

  case PLUS:
    {

      // First apply the inverse of Lm :
      // Solve Lm z = chi
      z[N5-1][rb[cb]] = psi[N5-1];
      for(int s=0; s < N5-1; s++){
	z[s][rb[cb]] = psi[s];

	// The factor of 1/2 is for the projection expression
	fact = Real(0.5)*l[s];
        z[N5-1][rb[cb]] -=  fact*(psi[s] - GammaConst<Ns,Ns*Ns-1>()*psi[s])  ;
       
      }
      
      //Now apply the inverse of L. Forward elimination 
      //
      // L z' = z
      z_prime[0][rb[cb]] = z[0];
      for(int s = 0; s < N5-1; s++) {
	// The factor of 1/2 is for the projection part
	fact = Real(0.5)*a[s];
	z_prime[s+1][rb[cb]] = z[s+1] - fact*(z_prime[s] + GammaConst<Ns,Ns*Ns-1>()*z_prime[s]);
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
	fact = Real(0.5)*b[s];

	z_prime[s][rb[cb]] = z[s] - fact*(z_prime[s+1] - GammaConst<Ns,Ns*Ns-1>()*z_prime[s+1]);
      }

      //Finally the inverse of Rm 
      chi[N5-1][rb[cb]] = z_prime[N5-1];
      for(int s=0; s < N5-1; s++) { 
	fact = Real(0.5)*r[s]; 

	chi[s][rb[cb]] = z_prime[s] - fact*(chi[N5-1] + GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
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

	// The factor of 1/2 is for the projection expression
	fact = Real(0.5)*r[s];
        z[N5-1][rb[cb]] -=  fact*(psi[s] + GammaConst<Ns,Ns*Ns-1>()*psi[s])  ;
       
      }
      
      //Now apply the inverse of U^{T}. Forward elimination 
      //
      // U^T z' = z . This is the same as L z' = z above but with 
      // a P+ => b P-

      z_prime[0][rb[cb]] = z[0];
      for(int s = 0; s < N5-1; s++) {
	// The factor of 1/2 is for the projection part
	fact = Real(0.5)*b[s];
	z_prime[s+1][rb[cb]] = z[s+1] - fact*(z_prime[s] - GammaConst<Ns,Ns*Ns-1>()*z_prime[s]);
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
	fact = Real(0.5)*a[s];
	z_prime[s][rb[cb]] = z[s] - fact*(z_prime[s+1] + GammaConst<Ns,Ns*Ns-1>()*z_prime[s+1]);
      }

      //Finally the inverse of Lm^T 
      //Same as the inverse of R above but with r P+ => l P-
      chi[N5-1][rb[cb]] = z_prime[N5-1];
      for(int s=0; s < N5-1; s++) { 
	fact = Real(0.5)*l[s]; 

	chi[s][rb[cb]] = z_prime[s] - fact*(chi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*chi[N5-1]);
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
					 const int cb) const 
{
  START_CODE();

  chi.resize(N5);

  switch ( isign ) 
  {
  case PLUS:
    {
      LatticeFermion tmp;
      LatticeFermion tmp1;
      LatticeFermion tmp2;
      Real fact1;
      Real fact2;
      int otherCB = (cb + 1)%2 ;
 
      fact1 = -Real(0.5)*b5[0];
      fact2 = -Real(0.25)*c5[0];
      tmp1[rb[otherCB]] = psi[1] - m_q*psi[N5-1];
      tmp2[rb[otherCB]] = psi[1] + m_q*psi[N5-1];
      
      tmp[rb[otherCB]] = fact1*psi[0] 
	+ fact2*(tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2);

      D.apply(chi[0], tmp, isign, cb);
      
      
      fact1 = -Real(0.5)*b5[N5-1];
      fact2 = -Real(0.25)*c5[N5-1];

      tmp1[rb[otherCB]]=psi[N5-2]-m_q*psi[0];
      tmp2[rb[otherCB]]=psi[N5-2]+m_q*psi[0];

      tmp[rb[otherCB]] = fact1*psi[N5-1]
	+ fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2 );

      D.apply(chi[N5-1], tmp, isign, cb);
      

      for(int s=1; s < N5-1; s++) { 
	fact1 = -Real(0.5)*b5[s];
	fact2 = -Real(0.25)*c5[s];

	tmp1[rb[otherCB]]=psi[s-1] + psi[s+1];
	tmp2[rb[otherCB]]=psi[s-1] - psi[s+1];
	
	tmp[rb[otherCB]]= fact1*psi[s] 
	  + fact2*(tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2) ;
	
	D.apply(chi[s], tmp, isign, cb);
      }
    }
    break ;
    
  case MINUS:
    { 
      multi1d<LatticeFermion> tmp_d(N5) ;

      for(int s=0; s<N5; s++){
	D.apply(tmp_d[s], psi[s], isign, cb);
      }

      LatticeFermion tmp1;
      LatticeFermion tmp2;


      Real one_quarter = Real(0.25);

      Real factb= -Real(0.5)*b5[0];
      Real fact  = m_q*c5[N5-1];
      tmp1[rb[cb]] = c5[1]*tmp_d[1] - fact*tmp_d[N5-1];
      tmp2[rb[cb]] = c5[1]*tmp_d[1] + fact*tmp_d[N5-1];

      chi[0][rb[cb]] = factb*tmp_d[0]
	- one_quarter*( tmp1 + GammaConst<Ns,Ns*Ns-1>()*tmp2 );



      factb = -Real(0.5)*b5[N5-1];
      fact = m_q * c5[0];
      
      tmp1[rb[cb]] = c5[N5-2]*tmp_d[N5-2] - fact * tmp_d[0];
      tmp2[rb[cb]] = c5[N5-2]*tmp_d[N5-2] + fact * tmp_d[0];

      chi[N5-1][rb[cb]] = factb * tmp_d[N5-1]
	- one_quarter * ( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );

      for(int s=1; s<N5-1; s++){
	factb = -Real(0.5) * b5[s];
	
	tmp1[rb[cb]] = c5[s-1]*tmp_d[s-1] + c5[s+1]*tmp_d[s+1];
	tmp2[rb[cb]] = c5[s-1]*tmp_d[s-1] - c5[s+1]*tmp_d[s+1];

	chi[s][rb[cb]] = factb*tmp_d[s] 
	  - one_quarter*( tmp1 - GammaConst<Ns, Ns*Ns-1>()*tmp2 );
	
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
  LatticeFermion tt ;
  D.apply(tt,psi,isign,0);
  D.apply(tt,psi,isign,1);
  Real fact = Real(0.5)*c5[s5];
  chi = f_minus[s5]*psi - fact*tt ;
}

}; // End Namespace Chroma

