/* $Id: eoprec_ovext_linop_array_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
*  \brief EvenOddPreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
*/

#include "chromabase.h"
#include "actions/ferm/linop/eoprec_ovext_linop_array_w.h"

namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup fermact
   *
   * \param fs            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   */
  void 
  EvenOddPrecOvExtLinOpArray::create(Handle< FermState<T,P,Q> > fs,
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


    Npoles = Npoles_;
    N5 = 2*Npoles_ + 1;

    Dslash.create(fs,N5);

    Real R = (Real(1) + Mass_)/(Real(1)-Mass_);
    Real alpha = b5_ + c5_;
    Real a5 = b5_ - c5_;

    Real QQ = Nd - OverMass_;
    A = -alpha*QQ;
    E = Real(2)*R  + (R*a5 + alpha*coeffP_)*QQ;

    Aprime = Real(0.5)*alpha;
    Eprime = -Real(0.5)*(R *a5 + coeffP_*alpha);

    B.resize(Npoles);
    C.resize(Npoles);
    D.resize(Npoles);
    Bprime.resize(Npoles);
    Cprime.resize(Npoles);
    Dprime.resize(Npoles);

    for(int p=0; p < Npoles; p++) { 
      B[p] = sqrt( rootQ_[p] * beta_[p] )* ( Real(2) + a5*QQ );
      D[p] = sqrt( resP_[p] * beta_[p] ) * ( Real(2) + a5*QQ );
      C[p] = alpha*beta_[p]*QQ;

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
   *
   * Flopcount = 10*Nc*Ns*(N5-1) + 2*Nc*Ns
   *           = 10*Nc*Ns*N5 - 8*Nc*Ns  = (10N5 - 8)*Nc*Ns
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
	// Lowest corner
	// 2*Nc*Ns flops/cbsite	
	chi[N5-1][rb[cb]] = E* ( GammaConst<Ns,Ns*Ns-1>()*psi[N5-1] );

	
	int p=0; 
	// ((N5-1)/2)*20*Nc*Ns flops/cbsite = 10*Nc*Ns*(N5-1) flops/cbsite
	for(int i=0; i < N5-1; i+=2, p++) { 
	  // 6*Nc*Ns flops/cbsite
	  chi[i][rb[cb]] = B[p]*psi[i+1] + A*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

	  // 6*Nc*Ns flops/cbsite
	  chi[i+1][rb[cb]] = B[p]*psi[i] + C[p]*(GammaConst<Ns,Ns*Ns-1>()*psi[i+1]);
	  // 4*Nc*Ns flops/cbsite
	  chi[i+1][rb[cb]] += D[p]*psi[N5-1];

	  // 4*Nc*Ns flops/cbsite
	  chi[N5-1][rb[cb]] -= D[p]*psi[i+1];
	}
      }
      break;
    case MINUS:
      {
	// Lowest corner
	// 2*Nc*Ns cbsite flops
	chi[N5-1][rb[cb]] = E*(GammaConst<Ns,Ns*Ns-1>()*psi[N5-1]);
	
	int p=0; 
	// 10 Nc*Ns*(N5-1) cbsite flops for loop
	for(int i=0; i < N5-1; i+=2, p++) { 
	  // 6*Nc*Ns cbsite flops
	  chi[i][rb[cb]] = B[p]*psi[i+1] + A*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

	  // 6*Nc*Ns cbsite flops
	  chi[i+1][rb[cb]] = B[p]*psi[i] + C[p]*(GammaConst<Ns,Ns*Ns-1>()*psi[i+1]);
	  // 4*Nc*Ns cbsite flops
	  chi[i+1][rb[cb]] -= D[p]*psi[N5-1];
	  
	  // 4*Nc*Ns cbsite flops
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


  /* Flopcount = N5*1320 + (10*N5-8)*Nc*Ns cbsite flops */
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
	multi1d<LatticeFermion> Dpsi(N5); moveToFastMemoryHint(Dpsi);

	// N5 * Dslash flops = N5 * 1320 cbsite flops
	Dslash.apply(Dpsi,psi,PLUS,cb);

     
	// Lowest corner
	// 2*Nc*Ns*cbsite flops
	chi[N5-1][rb[cb]] = Eprime*(GammaConst<Ns,Ns*Ns-1>()*Dpsi[N5-1]);
	
	int p=0; 
	// Total flops for loop: 10*(N5-1)*Nc*Ns cbsite flops
	for(int i=0; i < N5-1; i+=2, p++) { 

	  // 6*Nc*Ns*cbsite flops
	  chi[i][rb[cb]] = Bprime[p]*Dpsi[i+1] + Aprime*(GammaConst<Ns,Ns*Ns-1>()*Dpsi[i]);

	  // 6*Nc*Ns*cbsite flops
	  chi[i+1][rb[cb]] = Bprime[p]*Dpsi[i] + Cprime[p]*(GammaConst<Ns,Ns*Ns-1>()*Dpsi[i+1]);

	  // 4*Nc*Ns*cbsite flops
	  chi[i+1][rb[cb]] += Dprime[p]*Dpsi[N5-1];

	  // 4*Nc*Ns*cbsite flops
	  chi[N5-1][rb[cb]] -= Dprime[p]*Dpsi[i+1];
	}
      }
      break;
    case MINUS:
      {
	multi1d<LatticeFermion> tmp(N5); moveToFastMemoryHint(tmp);

	int otherCB = (cb + 1)%2;


	// Lowest corner
	// 2*Nc*Ns cbsite flops
	tmp[N5-1][rb[otherCB]] = Eprime*(GammaConst<Ns,Ns*Ns-1>()*psi[N5-1]);
	
	int p=0; 
	for(int i=0; i < N5-1; i+=2, p++) { 
	  // 6*Nc*Ns cbsite flops
	  tmp[i][rb[otherCB]] = Bprime[p]*psi[i+1] + Aprime*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);
	  
	  // 6*Nc*Ns cbsite flops
	  tmp[i+1][rb[otherCB]] = Bprime[p]*psi[i] + Cprime[p]*(GammaConst<Ns,Ns*Ns-1>()*psi[i+1]);

	  // 4*Nc*Ns cbsite flops
	  tmp[i+1][rb[otherCB]] -= Dprime[p]*psi[N5-1];
	  
	  // 4*Nc*Ns cbsite flops
	  tmp[N5-1][rb[otherCB]] += Dprime[p]*psi[i+1];
	}

	// N5 *1320 cbsite flops
	Dslash.apply(chi,tmp,MINUS,cb);
      }
      break;
          default: 
      QDPIO::cerr << "Unknown value for PLUS /MINUS: " << isign << endl;
      QDP_abort(1);
    };

    END_CODE();
  }

  /* Flopcount:  Npoles*10*Nc*Ns cbsite flops
               + Npoles*12*Nc*Ns cbsite flops
	       + 2*Nc*Ns cbsite flops
	       + Npoles*8*Nc*Ns cbsite flops
	       = Npoles*30*Nc*Ns + 2Nc*Ns cbsite flops
	       = (30Npoles + 2)*Nc*Ns cbsite flops
	       = (30(N5-1)/2 + 2)*Nc*Ns cbsite flops
	       = (15N5-13)Nc*Nc cbsite flops
  */

  void
  EvenOddPrecOvExtLinOpArray::applyDiagInv(multi1d<LatticeFermion>& psi,
					 const multi1d<LatticeFermion>& chi,
					 enum PlusMinus isign,
					 const int cb) const 
  {
    START_CODE();

    multi1d<LatticeFermion> tmp5(N5);   moveToFastMemoryHint(tmp5);

  
    LatticeFermion tmp4;                moveToFastMemoryHint(tmp4);
  
    Real ftmp;
    Real ftmp2;

    if( psi.size() != N5 ) { psi.resize(N5); }

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


    // Npoles * 10*Nc*Ns flops
    for(int i=0; i < 2*Npoles; i+=2) { 
      // Copy
      tmp5[i][rb[cb]] = chi[i];
      tmp5[i+1][rb[cb]] = chi[i+1];


      // Fixup last component
      // 6Nc*Ns cbsite flops
      ftmp = sign*D_bd_inv[i];
      ftmp2 = sign*D_bd_inv[i+1];

      tmp4[rb[cb]] =ftmp *chi[i] - ftmp2*(GammaConst<Ns,Ns*Ns-1>()*chi[i+1]);
      
      // 4Nc*Ns cbsite flops
      tmp5[N5-1][rb[cb]] -= tmp4;
      
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
    // Npoles*12*Nc*Ns flops
    for(int i=0; i < 2*Npoles; i+=2, p++) {
      // Do all the blocks
      // 6Nc*Ns flops
      ftmp=-Btilde[p];
      psi[i][rb[cb]] =  ftmp*tmp5[i+1] + Ctilde[p]*(GammaConst<Ns,Ns*Ns-1>()*tmp5[i]);

      // 6Nc*Ns flops
      psi[i+1][rb[cb]] = ftmp*tmp5[i] + Atilde[p]*(GammaConst<Ns,Ns*Ns-1>()*tmp5[i+1]);
    }

    // Do the last bit
    ftmp = Real(1)/S;

    // 2*Nc*Ns flops
    psi[N5-1][rb[cb]] = ftmp*(GammaConst<Ns,Ns*Ns-1>()*tmp5[N5-1]);


      
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
    
    tmp4[rb[cb]] = Real(1)*(GammaConst<Ns,Ns*Ns-1>()*psi[N5-1]);
    
    // Npoles * 8Nc*Ns flops
    for(int i=0; i < 2*Npoles; i+=2) {
      ftmp = sign*D_bd_inv[i];
      ftmp2 = sign*D_bd_inv[i+1];

      // 4*Nc*Ns flops
      psi[i][rb[cb]]  += ftmp*psi[N5-1];

      // 4*Nc*Ns flops
      psi[i+1][rb[cb]] -= ftmp2*tmp4;
	
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

} // End Namespace Chroma

