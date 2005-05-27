// $Id: prec_ovlap_contfrac5d_linop_array_opt_w.cc,v 1.1 2005-05-27 18:39:38 edwards Exp $
/*! \file
 *  \brief Optimized version of 5D continued frac. linop
 */


#include "chromabase.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/prec_ovlap_contfrac5d_linop_array_opt_w.h"


namespace Chroma 
{ 
  OptEvenOddPrecOvlapContFrac5DLinOpArray::OptEvenOddPrecOvlapContFrac5DLinOpArray(
    Handle<const ConnectState> state,
    const Real& _m_q,
    const Real& _OverMass,
    int _N5,
    const Real& _scale_fac,
    const multi1d<Real>& _alpha,
    const multi1d<Real>& _beta,
    const bool _isLastZeroP ) :
    m_q(_m_q), OverMass(_OverMass), N5(_N5), scale_fac(_scale_fac), 
    alpha(_alpha), beta(_beta), isLastZeroP(_isLastZeroP)
  {
    START_CODE();

    Handle< const DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > Ds(new WilsonDslash(state->getLinks()));
    Dslash  = Ds;  // Copy Handle -- M now owns dslash

    // The mass ratio
    Real mass = ( Real(1) + m_q ) / (Real(1) - m_q);

    // Now compute some coefficients.
    // First the beta_tilde_i
    // Basically this is beta[n]*Hsign*scale_fac
    // Now N5 is always odd. So the first Hsign is +
    // and the last one should also be
    // Hence at the end of this loop Hsign should be flipped from +->-
    beta_tilde.resize(N5);
    int Hsign = 1;
    for(int i=0; i < N5; i++) 
    { 
      // Flip Hsign
      beta_tilde[i] = beta[i]*Hsign*scale_fac; 

      /*
	QDPIO::cout << "beta["<<i<<"]=" << beta[i]
	<< "  Hsign=" << Hsign
	<< "  scale_fac=" << scale_fac 
	<< "  beta_tilde["<<i<<"]=" << beta_tilde[i] << endl;

      */
      Hsign = -Hsign;
    }

    // Sanity Check
    if ( Hsign > 0 ) {
      QDPIO::cerr << "Something is wrong. At the end of this loop"
		  << " Hsign should be -ve" << endl;
    }

    // Now the a_i's and b_i's
    a.resize(N5);
    for(int i=0; i < N5-1; i++) { 
      a[i] = beta_tilde[i]*(Nd - OverMass);
    }
    a[N5-1] = mass + (beta_tilde[N5-1]*(Nd - OverMass));

    /*
      QDPIO::cout << "Nd - OverMass = " << Nd- OverMass << endl;
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "a["<<i<<"]= " << a[i] << endl;
      }
    */
    // Now the d-s
    multi1d<Real> d(N5);
    invd.resize(N5);

    d[0] = a[0];
    invd[0] = Real(1)/d[0];

    for(int i=1; i < N5; i++) { 
      d[i] = a[i] - (alpha[i-1]*alpha[i-1])/d[i-1];
      invd[i] = Real(1)/d[i];
    }


    /*
      for(int i=0; i < N5; i++) { 
      QDPIO::cout << "d["<<i<<"]=" << d[i] << endl;
      }
    */

    // Now the u-s
    u.resize(N5-1);
    for(int i=0; i < N5-1; i++) { 
      u[i] = alpha[i]/d[i];
    }

    off_diag_coeff.resize(N5);
    for(int i=0; i < N5; i++) { 
      off_diag_coeff[i] = -Real(0.5)*beta_tilde[i];
    }
    /*
      for(int i=0; i < N5-1; i++) { 
      QDPIO::cout << "u["<<i<<"] = " << u[i] << endl;
      }
    */
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
   *
   *
   *  Flopcount: N5*6NcNs + (N5-2)*4NcNs = NcNs( 6N5 +4(N5-2) ) = (10N5-8) Nc Ns / cb_site
   */
  void 
  OptEvenOddPrecOvlapContFrac5DLinOpArray::applyDiag(multi1d<LatticeFermion>& chi, 
						     const multi1d<LatticeFermion>& psi, 
						     enum PlusMinus isign,
						     const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);
   
    // We don't care about the isign because our operator is Hermitian
    // Apply matrix
    //   [ A_0  B_0   0     ...                       ]  [ psi_0    ]
    //   [ B_0  A_1  B_1                   ...   ...  ]  [ psi_1    ]
    //   [  0   ...  ...     ...                 ...  ]  [ psi_2    ]
    //   [  ...    ...    0    B_N5-3  A_N5-2  B_N5-2 ]  [ psi_N5-2 ]
    //   [  ...    ...    ...   0      B_N5-2  A_N5-1 ]  [ psi_N5-1 ]

    // With A_i = gamma_5 a_i = a_i gamma_5
    // and  B_i = b_i I = alpha_i I

    LatticeFermion tmp;
//    int G5=Ns*Ns-1;

    // First 0ne 
    // Operation: chi[0][rb[cb]] = a[0] G5 psi[0] + alpha[0]*psi[1]
    //
    //  Useful flops: 6Nc Ns / site 
    // tmp[rb[cb]] = Gamma(G5)*psi[0];
    // chi[0][rb[cb]] = a[0]*tmp;
    
    if( N5 > 1 ) { 
      chi[0][rb[cb]] = alpha[0]*psi[1] + a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }
    else {
      chi[0][rb[cb]] = a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
    }

    // All the rest
    for(int i=1; i < N5; i++) { 

      // Operation: 
      //   N5 - 1 times:
      //    chi[i]  = alpha[i-1]*psi[i-1] + a[i] Gamma_5 *psi[i]
      //   N5 - 2 times:
      //    chi[i] += alpha[i]*psi[i+1];
      //  Useful flops = (N5-1) * 6NcNs + (N5-2)*4Nc*Ns

      /*
      // B_{i-1} psi_[i-1]
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1];

      // A_{i} psi[i] = a_{i} g_5 psi[i]
      tmp[rb[cb]] = Gamma(G5)*psi[i];
      chi[i][rb[cb]] += a[i]*tmp;
      */
      chi[i][rb[cb]] = alpha[i-1]*psi[i-1] + a[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-1) {
	chi[i][rb[cb]] += alpha[i]*psi[i+1];
      }
    }

    END_CODE();
  }


  //! Apply the inverse even-even (odd-odd)
  /*!
   * \ingroup linop
   *
   * Here we apply the LDU decomposition
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   * \param cb      checkerboard ( 0 | 1 )               (Read)

   *
   * Total flopcount: (N5-1)*4NcNs + 2NcNs + (N5-1)*6NcNs
   *                 = (N5-1)*10NcNs + 2NcNs
   *                 = (10N5-8) Nc Ns per_cb_site
   */
  void 
  OptEvenOddPrecOvlapContFrac5DLinOpArray::applyDiagInv(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
    if( chi.size() != N5 )  chi.resize(N5);

    multi1d<Fermion> y(N5);

//    LatticeFermion tmp;
//    const int G5 = Ns*Ns-1;

    for(int site=rb[cb].start(); site <= rb[cb].end(); ++site)
    {
      // Solve L y = psi
      {
	//y[0] = psi[0];
	REAL *In  = (REAL *) &(psi[0].elem(site).elem(0).elem(0).real());
	REAL* Out = (REAL *) &(y[0].elem().elem(0).elem(0).real());

	register int index_x = 0;
	register int index_z = 0;

	Out[index_z++] = In[index_x++]; // 00r
	Out[index_z++] = In[index_x++]; // 00i
	Out[index_z++] = In[index_x++]; // 01r
	Out[index_z++] = In[index_x++]; // 01i
	Out[index_z++] = In[index_x++]; // 02r
	Out[index_z++] = In[index_x++]; // 02i

	Out[index_z++] = In[index_x++]; // 10r
	Out[index_z++] = In[index_x++]; // 10i
	Out[index_z++] = In[index_x++]; // 11r
	Out[index_z++] = In[index_x++]; // 11i
	Out[index_z++] = In[index_x++]; // 12r
	Out[index_z++] = In[index_x++]; // 12i

	Out[index_z++] = In[index_x++]; // 20r
	Out[index_z++] = In[index_x++]; // 20i
	Out[index_z++] = In[index_x++]; // 21r
	Out[index_z++] = In[index_x++]; // 21i
	Out[index_z++] = In[index_x++]; // 22r
	Out[index_z++] = In[index_x++]; // 22i

	Out[index_z++] = In[index_x++]; // 30r
	Out[index_z++] = In[index_x++]; // 30i
	Out[index_z++] = In[index_x++]; // 31r
	Out[index_z++] = In[index_x++]; // 31i
	Out[index_z++] = In[index_x++]; // 32r
	Out[index_z++] = In[index_x++]; // 32i
      }

      // (N5-1)*4NcNs
      for(int i = 1; i < N5; i++) 
      {
	// tmp[rb[cb]] = Gamma(G5)*y[i-1];
	// y[i][rb[cb]] = psi[i] - u[i-1]*tmp;
	//y[i][rb[cb]] = psi[i] - u[i-1]*(GammaConst<Ns,Ns*Ns-1>()*y[i-1]);

	REAL  a =  u[i-1].elem().elem().elem().elem();
	REAL *InScale = (REAL *) &(y[i-1].elem().elem(0).elem(0).real());
	REAL *Add     = (REAL *) &(psi[i].elem(site).elem(0).elem(0).real());
	REAL* Out     = (REAL *) &(y[i].elem().elem(0).elem(0).real());
      
	// (Vector) out = (Vector) Add + (Scalar) (*scalep) * (Vector) P{+} InScale 
	// void xmayz_g5(REAL *Out,REAL *scalep,REAL *Add, REAL *InScale,int n_4vec)
	register int index_x = 0;
	register int index_y = 0;
	register int index_z = 0;
  
	// Spin Component 0 (AYPX)
	REAL x0r = Add[index_x++];
	REAL y0r = InScale[index_y++];
	REAL z0r = x0r - a*y0r;
	Out[index_z++] = z0r;
  
	REAL x0i = Add[index_x++];  
	REAL y0i = InScale[index_y++];
	REAL z0i = x0i - a*y0i;
	Out[index_z++] = z0i;

	REAL x1r = Add[index_x++];    
	REAL y1r = InScale[index_y++];
	REAL z1r = x1r - a*y1r;
	Out[index_z++] = z1r;

	REAL x1i = Add[index_x++];    
	REAL y1i = InScale[index_y++];
	REAL z1i = x1i - a*y1i;
	Out[index_z++] = z1i;
    
	REAL x2r = Add[index_x++];
	REAL y2r = InScale[index_y++];     
	REAL z2r = x2r - a*y2r;
	Out[index_z++] = z2r;
   
	REAL x2i = Add[index_x++]; 
	REAL y2i = InScale[index_y++];
	REAL z2i = x2i - a*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 1 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r - a*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i - a*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r - a*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i - a*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r - a*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i - a*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 2 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r + a*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i + a*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r + a*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i + a*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r + a*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i + a*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 3 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r + a*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i + a*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r + a*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i + a*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r + a*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i + a*y2i;  
	Out[index_z++] = z2i;
      } 

      // Invert diagonal piece  y <- D^{-1} y
      // N5 times: y = (1/d_i) gamma_5 y[i] 
      // Useful flops: N5 * 2NcNs flops / site
      // Rolled this into the bottom loop

      // Backsubstitute U chi = y

      // 2NcNs

      {
	//chi[N5-1][rb[cb]] = invd[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*y[N5-1]);
	REAL a =  invd[N5-1].elem().elem().elem().elem();
	REAL *In  = (REAL *) &(y[N5-1].elem().elem(0).elem(0).real());
	REAL *Out = (REAL *) &(chi[N5-1].elem(site).elem(0).elem(0).real());

	// (Vector) out = (Scalar) (*scalep) * (Vector) In
	//void scal_g5(REAL *Out, REAL *scalep, REAL *In, int n_4vec)
	register int index_x = 0;
	register int index_z = 0;

	// Spin Component 0
	REAL x0r = In[index_x++];
	REAL z0r = a*x0r;
	Out[index_z++] = z0r;
    
	REAL x0i = In[index_x++];
	REAL z0i = a*x0i;
	Out[index_z++] = z0i;
    
	REAL x1r = In[index_x++];
	REAL z1r = a*x1r;
	Out[index_z++] = z1r;
    
	REAL x1i = In[index_x++];
	REAL z1i = a*x1i;
	Out[index_z++] = z1i;
    
	REAL x2r = In[index_x++];     
	REAL z2r = a*x2r;
	Out[index_z++] = z2r;
    
	REAL x2i = In[index_x++];
	REAL z2i = a*x2i;
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = In[index_x++];
	z0r = a*x0r;
	Out[index_z++] = z0r;
    
	x0i = In[index_x++];
	z0i = a*x0i;
	Out[index_z++] = z0i;
    
	x1r = In[index_x++];
	z1r = a*x1r;
	Out[index_z++] = z1r;
    
	x1i = In[index_x++];
	z1i = a*x1i;
	Out[index_z++] = z1i;
    
	x2r = In[index_x++];     
	z2r = a*x2r;
	Out[index_z++] = z2r;
    
	x2i = In[index_x++];
	z2i = a*x2i;
	Out[index_z++] = z2i;

	// Spin Component 2
	x0r = In[index_x++];
	z0r = a*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = a*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = a*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = a*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = a*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = a*x2i;
	Out[index_z++] =-z2i;

	// Spin Component 3
	x0r = In[index_x++];
	z0r = a*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = a*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = a*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = a*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = a*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = a*x2i;
	Out[index_z++] =-z2i;
      }

      // N5-1 * 6NcNs
      for(int i = N5-2; i >= 0; i--) 
      {
	// tmp[rb[cb]] = Gamma(G5)*chi[i+1]
	// chi[i][rb[cb]] = y[i] - u[i]*tmp;
	// y[i][rb[cb]] = invd[i]*(GammaConst<Ns,Ns*Ns-1>()*y[i]);
	// chi[i][rb[cb]] = y[i] - u[i]*(GammaConst<Ns,Ns*Ns-1>()*chi[i+1]);

	//chi[i][rb[cb]] = GammaConst<Ns,Ns*Ns-1>()*(invd[i]*y[i]-u[i]*chi[i+1]); 
	REAL a = invd[i].elem().elem().elem().elem();
	REAL b = u[i].elem().elem().elem().elem();
	REAL *InScale = (REAL *) &(y[i].elem().elem(0).elem(0).real());
	REAL *Add     = (REAL *) &(chi[i+1].elem(site).elem(0).elem(0).real());
	REAL* Out     = (REAL *) &(chi[i].elem(site).elem(0).elem(0).real());

	// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*g5*vector)Add)
	//g5_axmbyz(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
	register int index_x = 0;
	register int index_y = 0;
	register int index_z = 0;

	// Spin Component 0 (AXPY3)
	REAL x0r = InScale[index_x++];
	REAL y0r = Add[index_y++];
	REAL z0r = a*x0r ;
	z0r -= b*y0r;
	Out[index_z++] = z0r;
    
	REAL x0i = InScale[index_x++];
	REAL y0i = Add[index_y++];
	REAL z0i = a*x0i;
	z0i -= b*y0i;
	Out[index_z++] = z0i;
    
	REAL x1r = InScale[index_x++];
	REAL y1r = Add[index_y++];
	REAL z1r = a*x1r ;
	z1r -= b*y1r;
	Out[index_z++] = z1r;
    
	REAL x1i = InScale[index_x++];
	REAL y1i = Add[index_y++];
	REAL z1i = a*x1i;
	z1i -= b*y1i;
	Out[index_z++] = z1i;
    
	REAL x2r = InScale[index_x++];     
	REAL y2r = Add[index_y++];
	REAL z2r = a*x2r ;
	z2r -= b*y2r;
	Out[index_z++] = z2r;
    
	REAL x2i = InScale[index_x++];
	REAL y2i = Add[index_y++];
	REAL z2i = a*x2i ;
	z2i -=  b*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = a*x0r;
	z0r -= b*y0r;
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = a*x0i;
	z0i -= b*y0i;
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = a*x1r ;
	z1r -= b*y1r;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = a*x1i;
	z1i -= b*y1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = a*x2r;
	z2r -=  b*y2r;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i = a*x2i;
	z2i -= b*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 2 (AXPY3)
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = b*y0r;
	z0r -= a*x0r ;
   
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = b*y0i;
	z0i -= a*x0i;   
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = b*y1r;
	z1r -= a*x1r ;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = b*y1i;
	z1i -= a*x1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = b*y2r;
	z2r -= a*x2r ;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i =  b*y2i;  
	z2i -= a*x2i ;
	Out[index_z++] = z2i;

	// Spin Component 3 (AXPY3)
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = b*y0r;
	z0r -= a*x0r ;
    
	Out[index_z++] = z0r;

	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = b*y0i;
	z0i -= a*x0i;   
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = b*y1r;
	z1r -= a*x1r ;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = b*y1i;
	z1i -= a*x1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = b*y2r;
	z2r -= a*x2r ;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i =  b*y2i;  
	z2i -= a*x2i ;
	Out[index_z++] = z2i;
      }
    }
  }


  //! Apply the off diagonal block
  /*!
   * \param chi     result     	                   (Modify)
   * \param psi     source     	                   (Read)
   * \param isign   Flag ( PLUS | MINUS )   	   (Read)
   * \param cb      checkerboard ( 0 | 1 )         (Read)
   *
   * Total flopcount: (N5-1)*(Dslash_cb_flops + 2NcNs) per cb_site (isLastZeroP==true)
   *                   N5*(Dslash_cb_flops+2NcNs) per cb_site      (isLastZeroP==false)
   *
   * 
   */
  void OptEvenOddPrecOvlapContFrac5DLinOpArray::applyOffDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi,
    enum PlusMinus isign,
    const int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    LatticeFermion tmp;
    
//    int G5 = Ns*Ns-1;

    // Optimisation... do up to the penultimate block...

    // (N5-1)( Dslash + 2NcNs) flops/site

    for(int i=0; i < N5-1; i++) { 
      // CB is CB of TARGET
      // gamma_5 Dslash is hermitian so I can ignore isign

      // Apply g5 Dslash
      Dslash->apply(tmp, psi[i], PLUS, cb);
  
      // chi[i][rb[cb]] = Gamma(G5)*tmp;

      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[i][rb[cb]] = off_diag_coeff[i]*(GammaConst<Ns,Ns*Ns-1>()*tmp);
    }

 
    // Only do last block if beta_tilde[i] is not zero
    // ( Dslash + 2NcNs) flops per site if done.
    if( !isLastZeroP ) {
      Dslash->apply(tmp, psi[N5-1], PLUS, cb);
      // chi[N5-1][rb[cb]] = Gamma(G5)*tmp;

      // Chi_i is now -(1/2) beta_tilde_i Dslash 
      chi[N5-1][rb[cb]] = off_diag_coeff[N5-1]*(GammaConst<Ns,Ns*Ns-1>()*tmp);
    }
    else { 
      chi[N5-1][rb[cb]] = zero;
    }
  
    END_CODE();
  }


  // Derivative of even-odd linop component
  /* 
   * This is a copy of the above applyOffDiag with the D.apply(...) replaced
   * by  D.deriv(ds_tmp,...) like calls.
   */
  void 
  OptEvenOddPrecOvlapContFrac5DLinOpArray::applyDerivOffDiag(multi1d<LatticeColorMatrix>& ds_u,
							     const multi1d<LatticeFermion>& chi, 
							     const multi1d<LatticeFermion>& psi, 
							     enum PlusMinus isign,
							     int cb) const 
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
						   
    LatticeFermion tmp;
    Real coeff;
    int G5 = Ns*Ns-1;

    switch (isign)
    {
    case PLUS:
      // Optimisation... do up to the penultimate block...
      for(int i=0; i < N5; i++) 
      {
	if (i == N5-1 && isLastZeroP) continue;

	// CB is CB of TARGET
	// consider case of gamma_5 Dslash
	tmp[rb[cb]] = Gamma(G5)*chi[i];

	// Multiply coefficient
	coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[cb]] *= coeff;

	// Apply g5 Dslash
	Dslash->deriv(ds_tmp, tmp, psi[i], PLUS, cb);
	ds_u += ds_tmp;
      }
      break;

    case MINUS:
      // Optimisation... do up to the penultimate block...
      for(int i=0; i < N5; i++) 
      {
	if (i == N5-1 && isLastZeroP) continue;

	// CB is CB of TARGET
	// consider case of Dslash^dag gamma_5
	tmp[rb[1-cb]] = Gamma(G5)*psi[i];

	// Multiply coefficient
	coeff = -Real(0.5)*beta_tilde[i];

	// Chi_i is now -(1/2) beta_tilde_i Dslash 
	tmp[rb[1-cb]] *= coeff;

	// Apply g5 Dslash
	Dslash->deriv(ds_tmp, chi[i], tmp, MINUS, cb);
	ds_u += ds_tmp;
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }

    END_CODE();
  }




  // THIS IS AN OPTIMIZED VERSION OF THE DERIVATIVE
  void 
  OptEvenOddPrecOvlapContFrac5DLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u,
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

