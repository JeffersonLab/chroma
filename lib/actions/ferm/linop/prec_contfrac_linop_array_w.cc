// $Id: prec_contfrac_linop_array_w.cc,v 1.1 2004-09-07 14:52:30 bjoo Exp $
/*! \file
 *  \brief  4D-style even-odd preconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_contfrac_linop_array_w.h"

// Check Conventions... Currently I (Kostas) am using Blum et.al.

//! Creation routine
/*! \ingroup fermact
 *
 * \param state_         Gauge Field ConenctState            (Read)
 * \param N5_            Extent of 5D                        (Read)
 * \param m_q_           Quark Mass                          (Read)
 * \param WilsonMass_    Wilson Mass (should be -ve)         (Read)
 * \param k_             Coefficients of Continued Fraction  (Read)
 * \param c_             Coefficients of Equiv Transform     (Read)
 */
void 
EvenOddPrecContFracLinOpArray::create(Handle< ConnectState > state_,
				      const int N5_,
				      const Real& m_q_,
				      const Real& WilsonMass_,
				      const multi1d<Real>& k_,
				      const multi1d<Real>& c_)  
{
  START_CODE();

  N5 = N5_;
  m_q = m_q_;
  WilsonMass = WilsonMass_;
  D.create(u_);

  // Auxiliary factor
  M_factor = (Real(Nd) + WilsonMass);

  // Sanity checking
  // k_.size() should be the same as N5
  if( k_.size() != N5 ) { 
    QDPIO::cerr << "EvenOddPrecContFracLinOpArray: k_.size()!=N5 "
		<< " k_.size()=" << k_.size() << " N5_=" << N5 << endl;
    QDP_abort(1);
  }

  // c.size() should be one less than N5 
  if( c_.size() != (N5-1) ) { 
    QDPIO::cerr << "EvenOddPrecContFracLinOpArray: c_.size() != (N5-1) "
		<< " c_.size()= " << c_.size() << " N5-1 " << N5-1 << endl;
    QDP_abort(1);
  }

  // If OK, copy the coefficients
  for(int i=0; i < N5; i++) { 
    k[i] = k_[i];
  }
  for(int i=0; i < N5-1; i++) { 
    c[i] = c_[i];
  }

  // Set up diagonal coeffs of tridiag matrix
  a.resize(N5);

  a[0] = ((Real(1) + m_q)/(Real(1) - m_q)) + k[0]*M_factor;

  for(int sgnH=-1, int i=1; i < N5; i++, sgnH *= -1) { 
    // a_i = (-1)^{i}*(c_{i-1}^2)*k[i]*(Nd + M_w)
    a[i] = Real(sgnH)*c[i-1]*c[i-1]*k[i]*M_factor;
  }

  // Set up off diagonal coeffs of tridiag matrix
  b.resize(N5-1);
  b[0] = c[0];
  for(int i=1; i < N5-1; i++) { 
    b[i] = c[i]*c[i-1];
  }

  // Compute the LDU decomposition coeffs.
  // First the coeffs of the D (most important)
  d.resize(N5);
  u_l.resize(N5-1);

  d[N5-1] = a[N5-1];
  for(int i=N5-2; i >=0; i--) {
    d[i] = a[i] - (b[i]*b[i]/d[i+1]);
    u_l[i] = b[i]/d[i+1];
  }

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
EvenOddPrecContFracLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
				   const multi1d<LatticeFermion>& psi, 
				   enum PlusMinus isign,
				   const int cb) const
{
  START_CODE();

  // We don't care about the isign because our operator is Hermitian
  // Apply matrix
  //   [ A_N5-1 B_N5-1    0     ...   ...  ]  [ psi_0    ]
  //   [ B_N5-2 A_N5-2 B_N5-1   ...   ...  ]  [ psi_1    ]
  //   [  0      ...    ...     ...   ...  ]  [ psi_2    ]
  //   [  ...    ...    B_2  A_2  B_1   0  ]  [ ...      ]
  //   [  ...    ...    0    B_1  A_1  B_0 ]  [ psi_N5-2 ]
  //   [  ...    ...    ...   0   B_0  A_0 ]  [ psi_N5-1 ]

  // With A_i = gamma_5 a_i = a_i gamma_5
  // and  B_i = b_i I

  LatticeFermion tmp;
  int G5=Ns*Ns-1;

  // chi[N5-1] = b[0] * psi[N5-2] + a[0] gamma_5 psi[N5-1]
  tmp[rb[cb]]  = Gamma(G5)*psi[ N5-1 ];
  chi[ N5-1 ][ rb[cb] ] = b[0]*psi[ N5-2 ]
  chi[ N5-1 ][ rb[cb] ] += a[0]*tmp;

  for(int s = N5-2, int index=1; s >= 0; s--, index++ ) { 

    // chi[s] = b[i]psi[s-1] + a[i]gamma_5 psi[s] + b[i-1] psi[s+1]
    chi[ s ][rb[cb]] = b[index]*psi[s-1];

    tmp[rb[cb]] = Gamma(G5)*psi[s];
    chi[ s ][rb[cb]] += a[index]*tmp;
    
    chi[ s ][rb[cb]] += b[index-1]*psi[s+1];
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
EvenOddPrecContFracLinOpArray::applyDiagInv(multi1d<LatticeFermion>& chi, 
				      const multi1d<LatticeFermion>& psi, 
				      enum PlusMinus isign,
				      const int cb) const
{
  START_CODE();

  multi1d<LatticeFermion> y(N5);
  LatticeFermion tmp;
  Real coeff;

  const int G5 = Ns*Ns-1;

  // Forward substitute L y = psi
  y[0][rb[cb]] = psi[0];

  for(int s = 1, index=N5-2; s < N5; s++, index--) { 
    tmp[rb[cb]] = Gamma(G5)*y[s-1];
    y[s][rb[cb]] = psi[s] - u_l[ index ]*tmp;
  }

  // Invert diagonal piece -- store back in y
  for(int s = N5-1,  index=0; s >= 0; s--, index++) { 
    tmp[rb[cb]] = Gamma(G5)*y[i];
    coeff = Real(1)/d[index];
    y[s][rb[cb]] = coeff*tmp;
  }

  // Backsubstitute U chi = y
  chi[N5-1][rb[cb]] = y[N5-1];

  for(int s = N5-2, index=0; s >= 0; s--; index++) {
    tmp = Gamma(G5)*y[s+1];
    chi[s][rb[cb]] = y[s] - u_l[index]*tmp;
  }

  //Done! That was not that bad after all....
  //See, I told you so...
  END_CODE();
}
