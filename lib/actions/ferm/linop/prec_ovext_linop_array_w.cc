/* $Id: prec_ovext_linop_array_w.cc,v 1.1 2005-04-22 13:27:42 bjoo Exp $
/*! \file
*  \brief EvenOddPreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
*/

#include "chromabase.h"
#include "actions/ferm/linop/prec_ovext_linop_array_w.h"

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
  EvenOddPrecOvExtLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
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

    START_CODE();

    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Ds(new WilsonDslash(u_));

    Dslash  = Ds;  // Copy Handle -- M now owns dslash

    Npoles = Npoles_;
    N5 = 2*Npoles_ + 1;

    Real R = (Real(1) + Mass_)/(Real(1)-Mass_);
    Real alpha = b5_ + c5_;
    Real a5 = b5_ - c5_;

    Q = Nd - OverMass_;
    A = -alpha*Q;
    E = Real(2)*R  + (R*a5 + alpha*coeffP_)*Q;

    Aprime = Real(0.5)*alpha;
    Eprime = -Real(0.5)*(R *a5 + coeffP_*alpha);

    B.resize(Npoles);
    C.resize(Npoles);
    D.resize(Npoles);
    Bprime.resize(Npoles);
    Cprime.resize(Npoles);
    Dprime.resize(Npoles);

    for(int p=0; p < Npoles; p++) { 
      B[p] = sqrt( rootQ_[p] * beta_[p] )* ( Real(2) + a5*Q );
      D[p] = sqrt( resP_[p] * beta_[p] ) * ( Real(2) + a5*Q );
      C[p] = alpha*beta_[p]*Q;

      Bprime[p] = -Real(0.5)*a5*sqrt( rootQ_[p] * beta_[p] );
      Dprime[p] = -Real(0.5)*a5*sqrt( resP_[p] * beta_[p] );
      Cprime[p] = -Real(0.5)*alpha*beta_[p];
    }

    multi1d<Real> detBlock(Npoles);
    for(int p=0; p < Npoles; p++) { 
      detBlock[p] = A*C[p] - B[p]*B[p];
    }

    Atilde.resize(Npoles);
    Btilde.resize(Npoles);
    Ctilde.resize(Npoles);

    for(int p=0; p < Npoles; p++) { 
      Atilde[p] = A/detBlock[p];
      Btilde[p] = B[p]/detBlock[p];
      Ctilde[p] = C[p]/detBlock[p];
    }

    S = E;
    for(int p=0; p < Npoles; p++) { 
      S += D[p]*D[p]*Atilde[p];
    }

    D_bd_inv.resize(2*Npoles);

    int pole=0;
    for(int i=0; i < 2*Npoles; i+=2, pole++) { 
      D_bd_inv[i] = D[pole]*Btilde[pole];
      D_bd_inv[i+1] = D[pole]*Atilde[pole];
    }

    END_CODE();
  }

  //! Apply the even-even (odd-odd) coupling piece of the domain-wall fermion operator
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )          (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   */
  void 
  EvenOddPrecOvExtLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
					const multi1d<LatticeFermion>& psi, 
					enum PlusMinus isign,
					const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);
    int G5 = Ns*Ns - 1;

    switch( isign ) { 
    case PLUS: 
      {
	LatticeFermion tmp4;
	tmp4[rb[cb]] = Gamma(G5)*psi[N5-1];

	// Lowest corner
	chi[N5-1][rb[cb]] = E*tmp4;
	
	int p=0; 
	for(int i=0; i < N5-1; i+=2, p++) { 
	  tmp4[rb[cb]] = Gamma(G5)*psi[i];
	  chi[i][rb[cb]] = A*tmp4 + B[p]*psi[i+1];

	  tmp4[rb[cb]] = Gamma(G5)*psi[i+1];
	  chi[i+1][rb[cb]] = B[p]*psi[i] + C[p]*tmp4;
	  chi[i+1][rb[cb]] += D[p]*psi[N5-1];

	  chi[N5-1][rb[cb]] -= D[p]*psi[i+1];
	}
      }
      break;
    case MINUS:
      {
	LatticeFermion tmp4;
	tmp4[rb[cb]] = Gamma(G5)*psi[N5-1];

	// Lowest corner
	chi[N5-1][rb[cb]] = E*tmp4;
	
	int p=0; 
	for(int i=0; i < N5-1; i+=2, p++) { 
	  tmp4[rb[cb]] = Gamma(G5)*psi[i];
	  chi[i][rb[cb]] = A*tmp4 + B[p]*psi[i+1];

	  tmp4[rb[cb]] = Gamma(G5)*psi[i+1];
	  chi[i+1][rb[cb]] = B[p]*psi[i] + C[p]*tmp4;
	  chi[i+1][rb[cb]] -= D[p]*psi[N5-1];

	  chi[N5-1][rb[cb]] += D[p]*psi[i+1];
	}
      }
      break;
          default: 
      QDPIO::cerr << "Unknown value for PLUS /MINUS: " << isign << endl;
      QDP_abort(1);
    };

    END_CODE();
  }

  void 
  EvenOddPrecOvExtLinOpArray::applyOffDiag(multi1d<LatticeFermion>& chi, 
					const multi1d<LatticeFermion>& psi, 
					enum PlusMinus isign,
					const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);
    int G5 = Ns*Ns - 1;

    switch( isign ) { 
    case PLUS: 
      {
	multi1d<LatticeFermion> Dpsi(N5);
	for(int i=0; i < N5; i++) { 
	  Dslash->apply(Dpsi[i], psi[i], PLUS, cb);
	}

	LatticeFermion tmp4;
	tmp4[rb[cb]] = Gamma(G5)*Dpsi[N5-1];

	// Lowest corner
	chi[N5-1][rb[cb]] = Eprime*tmp4;
	
	int p=0; 
	for(int i=0; i < N5-1; i+=2, p++) { 

	  tmp4[rb[cb]] = Gamma(G5)*Dpsi[i];
	  chi[i][rb[cb]] = Aprime*tmp4 + Bprime[p]*Dpsi[i+1];

	  tmp4[rb[cb]] = Gamma(G5)*Dpsi[i+1];
	  chi[i+1][rb[cb]] = Bprime[p]*Dpsi[i] + Cprime[p]*tmp4;
	  chi[i+1][rb[cb]] += Dprime[p]*Dpsi[N5-1];

	  chi[N5-1][rb[cb]] -= Dprime[p]*Dpsi[i+1];
	}
      }
      break;
    case MINUS:
      {
	multi1d<LatticeFermion> tmp(N5);
	int otherCB = (cb + 1)%2;


	LatticeFermion tmp4;
	tmp4[rb[otherCB]] = Gamma(G5)*psi[N5-1];

	// Lowest corner
	tmp[N5-1][rb[otherCB]] = Eprime*tmp4;
	
	int p=0; 
	for(int i=0; i < N5-1; i+=2, p++) { 
	  tmp4[rb[otherCB]] = Gamma(G5)*psi[i];
	  tmp[i][rb[otherCB]] = Aprime*tmp4 + Bprime[p]*psi[i+1];

	  tmp4[rb[otherCB]] = Gamma(G5)*psi[i+1];
	  tmp[i+1][rb[otherCB]] = Bprime[p]*psi[i] + Cprime[p]*tmp4;
	  tmp[i+1][rb[otherCB]] -= Dprime[p]*psi[N5-1];

	  tmp[N5-1][rb[cb]] += Dprime[p]*psi[i+1];
	}

	for(int i=0; i < N5; i++) { 
	  Dslash->apply(chi[i], tmp[i], MINUS, cb);
	}
      }
      break;
          default: 
      QDPIO::cerr << "Unknown value for PLUS /MINUS: " << isign << endl;
      QDP_abort(1);
    };

    END_CODE();
  }


  void
  EvenOddPrecOvExtLinOpArray::applyDiagInv(multi1d<LatticeFermion>& psi,
					 const multi1d<LatticeFermion>& chi,
					 enum PlusMinus isign,
					 const int cb) const 
  {
    START_CODE();

    multi1d<LatticeFermion> tmp5(N5);
    LatticeFermion tmp4;
    Real ftmp;

    if( psi.size() != N5 ) { psi.resize(N5); }

    int G5 = Ns*Ns-1;

    // Apply tmp5 = L^{-1} chi
    // L only has components along the last row, so this is a wholesale
    // copy.
    // The last row looks like:
    //
    // [ (+/-) D_p Btilde_p, (-/+) D_p Atilde_p g_5, ... , 1 ]
    //=[ (+/-) D_bd_inv[0], (-/+) D_bd_inv[1] g_5 ,.. , 1 ]
    //
    // (+/-) <=> isign=(PLUS, MINUS)
    //
    // And we subtract it ie:
    //  tmp5[N5-1] = chi[N5-1] - (+/-) term1 - (-/+) term2 ...
    //             =           - sgn term1 + sgn term2
    //
    //  with sgn = (+1/-1) <=> isign = PLUS/MINUS

    // Start off.
    tmp5[N5-1][rb[cb]] = chi[N5-1];

    Real sign;
    switch (isign) {
    case PLUS:
      sign = Real(1);
      break;
    case MINUS:
       sign = Real(-1);
       break;
    }

    for(int i=0; i < 2*Npoles; i+=2) { 
      // Copy
      tmp5[i][rb[cb]] = chi[i];
      tmp5[i+1][rb[cb]] = chi[i+1];

      // Fixup last component
      ftmp = sign*D_bd_inv[i];
      tmp5[N5-1][rb[cb]] -= ftmp*chi[i];
      tmp4[rb[cb]] = Gamma(G5)*chi[i+1];
      ftmp = sign*D_bd_inv[i+1];
      tmp5[N5-1][rb[cb]] += ftmp*tmp4;
    }

    // Apply Inverse of Block Diagonal Matrix to tmp5 
    //
    // [ A g_5  B_0                                   ]
    // [ B_0    C_0 g_5                               ]
    // [                A g_5  B_1                    ]
    // [                B_0    C_1 g_5                ]
    // [                               ...            ]
    // [                                        S g_5 ]
    //
    // which is
    // [ Ctilde_0 g_5  - Btilde_0                                       ]
    // [ -Btilde_0   Atilde_0 g_5                                       ] 
    // [                          Ctilde_1 g5  -Btilde_1                ]
    // [                          -Btilde_1    Atilde_1 g5              ]
    // [                                                   ...          ]
    // [                                                       (1/S)g_5 ]
    //
    // with Atilde_p = (1/d_p) A
    //      Btilde_p = (1/d_p) B_p,  Ctilde_p = (1/d_p)C_p
    //
    //      S = E + sum_p D^2_p A tilde_p
    // 
    //      d_p = det [ A g_5  B_p     ]
    //                [ B_p    C_p g_5 ]
    // 
    //          = A C_p - B_p^2
    //
    //  Atilde_p, Btilde_p, Ctilde_p and S are precomputed in the constructor.
    //
    // This matrix is Hermitian.


    int p=0;
    for(int i=0; i < 2*Npoles; i+=2, p++) {
      // Do all the blocks
      tmp4[rb[cb]] = Gamma(G5)*tmp5[i];
      psi[i][rb[cb]] = Ctilde[p]*tmp4 - Btilde[p]*tmp5[i+1];
      tmp4[rb[cb]] = Gamma(G5)*tmp5[i+1];
      psi[i+1][rb[cb]] = Atilde[p]*tmp4 - Btilde[p]*tmp5[i];
    }
    // Do the last bit
    tmp4[rb[cb]] = Gamma(G5)*tmp5[N5-1];
    ftmp = Real(1)/S;
    psi[N5-1][rb[cb]] = ftmp*tmp4;


      
    // Now apply: [ 1 .....                (-/+) Btilde_0 D_0    ] ^ {-1}
    //            [ 0 1                    (+/-) Atilde_0 D_0 G5 ]
    //            [ 0 0 1                  (-/+) Btilde_1 D_1    ]
    //            [ 0 0 0 1                (+/-) Atilde_1 D_1 G5 ]
    //            [ ...                               ...        ]
    //            [                                       1      ]
    // 
    // to tmp5_2 with backsubstitution. The elements in the rightmost
    // row are stored in D_bd_inv, but the sign is flipped compared to 
    // the forward substitution case. In either case we subtract
    // the term from the rightmost col ie:
    //  psi[ i ] = tmp5_2[i]    - (-/+) term * psi[N5-1];
    //           = tmp5_2[i]  + sign * term * psi[N5-1];
    //
    //  psi[i+1] = tmp5_2[i+1]  - (+/-) term * psi[N5-1];
    //           = tmp5_2[i+1]  - sign * term * psi[N5-1];
    //
    //  with sign = 1/-1 <=> isign = PLUS/MINUS
    // (ie signs flip in opposite order to forward subsitution)


    // Backsubstitution:
    // First piece already done

    // Only need this so precompute it out here.
    tmp4 = Gamma(G5)*psi[N5-1];

    for(int i=0; i < 2*Npoles; i+=2) {
      ftmp = sign*D_bd_inv[i];
      psi[i][rb[cb]]  += ftmp*psi[N5-1];

      ftmp = sign*D_bd_inv[i+1];
      psi[i+1][rb[cb]] -= ftmp*tmp4;
    }


    // And we are done. That was not so bad now was it?
    END_CODE();
  }

  void 
  EvenOddPrecOvExtLinOpArray::applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u, 
						const multi1d<LatticeFermion>& chi, 
						const multi1d<LatticeFermion>& psi, 
						enum PlusMinus isign,
						int cb) const
  {
    START_CODE();

    QDPIO::cout << "Not yet implemented " << endl;
    QDP_abort(1);    

    END_CODE();
  }
  //! Derivative
  void 
  EvenOddPrecOvExtLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign) const
  {
    START_CODE();

    QDPIO::cout << "Not yet implemented " << endl;
    QDP_abort(1);

    END_CODE();
  }
}; // End Namespace Chroma

