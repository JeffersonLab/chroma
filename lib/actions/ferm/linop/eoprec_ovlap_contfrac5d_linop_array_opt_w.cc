// $Id: eoprec_ovlap_contfrac5d_linop_array_opt_w.cc,v 3.2 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief Optimized version of 5D continued frac. linop
 */


#include "chromabase.h"
#include "actions/ferm/linop/eoprec_ovlap_contfrac5d_linop_array_opt_w.h"


namespace Chroma 
{ 


#if 0
  // FOR SOME UNKNOWN REASON, THIS CODE IS MUCH SLOWER THAN THE QDP VERSION!!


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
  OptEvenOddPrecOvlapContFrac5DLinOpArray::applyDiag(
    multi1d<LatticeFermion>& chi, 
    const multi1d<LatticeFermion>& psi, 
    enum PlusMinus isign,
    const int cb) const
  {
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

    // First 0ne 
    // Operation: chi[0][rb[cb]] = a[0] G5 psi[0] + alpha[0]*psi[1]
    //
    //  Useful flops: 6Nc Ns / site 
    // tmp[rb[cb]] = Gamma(G5)*psi[0];
    // chi[0][rb[cb]] = a[0]*tmp;
    
    const int start = rb[cb].start();
    const int end   = rb[cb].end();

    if( N5 > 1 ) { 
      //chi[0][rb[cb]] = alpha[0]*psi[1] + a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
      
      REAL scalep   = alpha[0].elem().elem().elem().elem();
      REAL scalep2  = a[0].elem().elem().elem().elem();
      REAL *InScale = (REAL *) &(psi[1].elem(start).elem(0).elem(0).real());
      REAL *Add     = (REAL *) &(psi[0].elem(start).elem(0).elem(0).real());
      REAL *Out     = (REAL *) &(chi[0].elem(start).elem(0).elem(0).real());

      // (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*g5*vector)Add)
      //void axpbyz_g5(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
      register int index_x = 0;
      register int index_y = 0;
      register int index_z = 0;
	
      for(int site=start; site <= end; ++site)
      {
	// Spin Component 0 (AXPY3)
	REAL x0r = InScale[index_x++];
	REAL y0r = Add[index_y++];
	REAL z0r = scalep*x0r ;
	z0r += scalep2*y0r;
	Out[index_z++] = z0r;
    
	REAL x0i = InScale[index_x++];
	REAL y0i = Add[index_y++];
	REAL z0i = scalep*x0i;
	z0i += scalep2*y0i;
	Out[index_z++] = z0i;
    
	REAL x1r = InScale[index_x++];
	REAL y1r = Add[index_y++];
	REAL z1r = scalep*x1r ;
	z1r += scalep2*y1r;
	Out[index_z++] = z1r;
    
	REAL x1i = InScale[index_x++];
	REAL y1i = Add[index_y++];
	REAL z1i = scalep*x1i;
	z1i += scalep2*y1i;
	Out[index_z++] = z1i;
    
	REAL x2r = InScale[index_x++];     
	REAL y2r = Add[index_y++];
	REAL z2r = scalep*x2r ;
	z2r += scalep2*y2r;
	Out[index_z++] = z2r;
    
	REAL x2i = InScale[index_x++];
	REAL y2i = Add[index_y++];
	REAL z2i = scalep*x2i ;
	z2i +=  scalep2*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep*x0r;
	z0r += scalep2*y0r;
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep*x0i;
	z0i += scalep2*y0i;
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep*x1r ;
	z1r += scalep2*y1r;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep*x1i;
	z1i += scalep2*y1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep*x2r;
	z2r +=  scalep2*y2r;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i = scalep*x2i;
	z2i += scalep2*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 2 (AXPY3)
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep*x0r ;
	z0r -= scalep2*y0r;
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep*x0i;
	z0i -= scalep2*y0i;
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep*x1r ;
	z1r -= scalep2*y1r;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep*x1i;
	z1i -= scalep2*y1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep*x2r ;
	z2r -= scalep2*y2r;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i = scalep*x2i ;
	z2i -=  scalep2*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 3
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep*x0r;
	z0r -= scalep2*y0r;
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep*x0i;
	z0i -= scalep2*y0i;
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep*x1r ;
	z1r -= scalep2*y1r;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep*x1i;
	z1i -= scalep2*y1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep*x2r;
	z2r -=  scalep2*y2r;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i = scalep*x2i;
	z2i -= scalep2*y2i;  
	Out[index_z++] = z2i;
      }
    }
    else 
    {
      //chi[0][rb[cb]] = a[0]*(GammaConst<Ns,Ns*Ns-1>()*psi[0]);
      REAL scalep =  a[0].elem().elem().elem().elem();
      REAL *In    = (REAL *) &(psi[0].elem(start).elem(0).elem(0).real());
      REAL *Out   = (REAL *) &(chi[0].elem(start).elem(0).elem(0).real());

      // (Vector) out = (Scalar) (*scalep) * (Vector) In
      //void scal_g5(REAL *Out, REAL *scalep, REAL *In, int n_4vec)
      register int index_x = 0;
      register int index_z = 0;

      for(int site=start; site <= end; ++site)
      {
	// Spin Component 0
	REAL x0r = In[index_x++];
	REAL z0r = scalep*x0r;
	Out[index_z++] = z0r;
    
	REAL x0i = In[index_x++];
	REAL z0i = scalep*x0i;
	Out[index_z++] = z0i;
    
	REAL x1r = In[index_x++];
	REAL z1r = scalep*x1r;
	Out[index_z++] = z1r;
    
	REAL x1i = In[index_x++];
	REAL z1i = scalep*x1i;
	Out[index_z++] = z1i;
    
	REAL x2r = In[index_x++];     
	REAL z2r = scalep*x2r;
	Out[index_z++] = z2r;
    
	REAL x2i = In[index_x++];
	REAL z2i = scalep*x2i;
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] = z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] = z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] = z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] = z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] = z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
	Out[index_z++] = z2i;

	// Spin Component 2
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
	Out[index_z++] =-z2i;

	// Spin Component 3
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
	Out[index_z++] =-z2i;
      }
    }

    // All the rest
    for(int i=1; i < N5; i++) 
    {
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
      {
	//chi[i][rb[cb]] = alpha[i-1]*psi[i-1] + a[i]*(GammaConst<Ns,Ns*Ns-1>()*psi[i]);

	REAL scalep   = alpha[i-1].elem().elem().elem().elem();
	REAL scalep2  = a[i].elem().elem().elem().elem();
	REAL *InScale = (REAL *) &(psi[i-1].elem(start).elem(0).elem(0).real());
	REAL *Add     = (REAL *) &(psi[i].elem(start).elem(0).elem(0).real());
	REAL *Out     = (REAL *) &(chi[i].elem(start).elem(0).elem(0).real());

	// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (scalep2)*g5*vector)Add)
	//void axpbyz_g5(REAL *Out,REAL *scalep,REAL *InScale, REAL *scalep2, REAL *Add,int n_4vec)
	register int index_x = 0;
	register int index_y = 0;
	register int index_z = 0;

	for(int site=start; site <= end; ++site)
	{
	  // Spin Component 0 (AXPY3)
	  REAL x0r = InScale[index_x++];
	  REAL y0r = Add[index_y++];
	  REAL z0r = scalep*x0r ;
	  z0r += scalep2*y0r;
	  Out[index_z++] = z0r;
    
	  REAL x0i = InScale[index_x++];
	  REAL y0i = Add[index_y++];
	  REAL z0i = scalep*x0i;
	  z0i += scalep2*y0i;
	  Out[index_z++] = z0i;
    
	  REAL x1r = InScale[index_x++];
	  REAL y1r = Add[index_y++];
	  REAL z1r = scalep*x1r ;
	  z1r += scalep2*y1r;
	  Out[index_z++] = z1r;
    
	  REAL x1i = InScale[index_x++];
	  REAL y1i = Add[index_y++];
	  REAL z1i = scalep*x1i;
	  z1i += scalep2*y1i;
	  Out[index_z++] = z1i;
    
	  REAL x2r = InScale[index_x++];     
	  REAL y2r = Add[index_y++];
	  REAL z2r = scalep*x2r ;
	  z2r += scalep2*y2r;
	  Out[index_z++] = z2r;
    
	  REAL x2i = InScale[index_x++];
	  REAL y2i = Add[index_y++];
	  REAL z2i = scalep*x2i ;
	  z2i +=  scalep2*y2i;  
	  Out[index_z++] = z2i;

	  // Spin Component 1
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r;
	  z0r += scalep2*y0r;
	  Out[index_z++] = z0r;
    
	  x0i = InScale[index_x++];
	  y0i = Add[index_y++];
	  z0i = scalep*x0i;
	  z0i += scalep2*y0i;
	  Out[index_z++] = z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r ;
	  z1r += scalep2*y1r;
	  Out[index_z++] = z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i;
	  z1i += scalep2*y1i;
	  Out[index_z++] = z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r;
	  z2r +=  scalep2*y2r;
	  Out[index_z++] = z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i;
	  z2i += scalep2*y2i;  
	  Out[index_z++] = z2i;

	  // Spin Component 2 (AXPY3)
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r ;
	  z0r -= scalep2*y0r;
	  Out[index_z++] = z0r;
    
	  x0i = InScale[index_x++];
	  y0i = Add[index_y++];
	  z0i = scalep*x0i;
	  z0i -= scalep2*y0i;
	  Out[index_z++] = z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r ;
	  z1r -= scalep2*y1r;
	  Out[index_z++] = z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i;
	  z1i -= scalep2*y1i;
	  Out[index_z++] = z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r ;
	  z2r -= scalep2*y2r;
	  Out[index_z++] = z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i ;
	  z2i -=  scalep2*y2i;  
	  Out[index_z++] = z2i;

	  // Spin Component 3
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r;
	  z0r -= scalep2*y0r;
	  Out[index_z++] = z0r;
    
	  x0i = InScale[index_x++];
	  y0i = Add[index_y++];
	  z0i = scalep*x0i;
	  z0i -= scalep2*y0i;
	  Out[index_z++] = z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r ;
	  z1r -= scalep2*y1r;
	  Out[index_z++] = z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i;
	  z1i -= scalep2*y1i;
	  Out[index_z++] = z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r;
	  z2r -=  scalep2*y2r;
	  Out[index_z++] = z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i;
	  z2i -= scalep2*y2i;  
	  Out[index_z++] = z2i;
	}
      }

      // When i hits N5-1, we don't have the B_N5-1 term
      if(i < N5-1) 
      {
	//chi[i][rb[cb]] += alpha[i]*psi[i+1];

	REAL scalep   = alpha[i].elem().elem().elem().elem();
	REAL* InScale = (REAL *)&(psi[i+1].elem(start).elem(0).elem(0).real());
	REAL* Add     = (REAL *)&(chi[i].elem(start).elem(0).elem(0).real());
	REAL* Out     = (REAL *)&(chi[i].elem(start).elem(0).elem(0).real());

	//vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
	// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (Vector) Add
	//void vaxpy3(REAL *Out,REAL *scalep,REAL *InScale, REAL *Add,int n_3vec)
	register int index_x = 0;
	register int index_y = 0;
	register int index_z = 0;

	for(int site=start; site <= end; ++site)
	{
	  // Spin Component 0 (AXPY3)
	  REAL x0r = InScale[index_x++];
	  REAL y0r = Add[index_y++];
	  REAL z0r = scalep*x0r + y0r;
	  Out[index_z++] =(REAL) z0r;
    
	  REAL x0i = InScale[index_x++];
 	  REAL y0i = Add[index_y++];
	  REAL z0i = scalep*x0i + y0i;
	  Out[index_z++] =(REAL) z0i;
    
	  REAL x1r = InScale[index_x++];
	  REAL y1r = Add[index_y++];
	  REAL z1r = scalep*x1r + y1r;
	  Out[index_z++] = (REAL)z1r;
    
	  REAL x1i = InScale[index_x++];
	  REAL y1i = Add[index_y++];
	  REAL z1i = scalep*x1i + y1i;
	  Out[index_z++] = (REAL)z1i;
    
	  REAL x2r = InScale[index_x++];     
	  REAL y2r = Add[index_y++];
	  REAL z2r = scalep*x2r + y2r;
	  Out[index_z++] = (REAL)z2r;
    
	  REAL x2i = InScale[index_x++];
	  REAL y2i = Add[index_y++];
	  REAL z2i = scalep*x2i + y2i;  
	  Out[index_z++] = (REAL)z2i;

	  // Spin Component 1 (AXPY3)
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r + y0r;
	  Out[index_z++] =(REAL) z0r;
    
	  x0i = InScale[index_x++];
 	  y0i = Add[index_y++];
	  z0i = scalep*x0i + y0i;
	  Out[index_z++] =(REAL) z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r + y1r;
	  Out[index_z++] = (REAL)z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i + y1i;
	  Out[index_z++] = (REAL)z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r + y2r;
	  Out[index_z++] = (REAL)z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i + y2i;  
	  Out[index_z++] = (REAL)z2i;

	  // Spin Component 2 (AXPY3)
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r + y0r;
	  Out[index_z++] =(REAL) z0r;
    
	  x0i = InScale[index_x++];
 	  y0i = Add[index_y++];
	  z0i = scalep*x0i + y0i;
	  Out[index_z++] =(REAL) z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r + y1r;
	  Out[index_z++] = (REAL)z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i + y1i;
	  Out[index_z++] = (REAL)z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r + y2r;
	  Out[index_z++] = (REAL)z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i + y2i;  
	  Out[index_z++] = (REAL)z2i;

	  // Spin Component 3 (AXPY3)
	  x0r = InScale[index_x++];
	  y0r = Add[index_y++];
	  z0r = scalep*x0r + y0r;
	  Out[index_z++] =(REAL) z0r;
    
	  x0i = InScale[index_x++];
 	  y0i = Add[index_y++];
	  z0i = scalep*x0i + y0i;
	  Out[index_z++] =(REAL) z0i;
    
	  x1r = InScale[index_x++];
	  y1r = Add[index_y++];
	  z1r = scalep*x1r + y1r;
	  Out[index_z++] = (REAL)z1r;
    
	  x1i = InScale[index_x++];
	  y1i = Add[index_y++];
	  z1i = scalep*x1i + y1i;
	  Out[index_z++] = (REAL)z1i;
    
	  x2r = InScale[index_x++];     
	  y2r = Add[index_y++];
	  z2r = scalep*x2r + y2r;
	  Out[index_z++] = (REAL)z2r;
    
	  x2i = InScale[index_x++];
	  y2i = Add[index_y++];
	  z2i = scalep*x2i + y2i;  
	  Out[index_z++] = (REAL)z2i;
	}
      }
    }
  }
#endif


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

	REAL  scalep  =  u[i-1].elem().elem().elem().elem();
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
	REAL z0r = x0r - scalep*y0r;
	Out[index_z++] = z0r;
  
	REAL x0i = Add[index_x++];  
	REAL y0i = InScale[index_y++];
	REAL z0i = x0i - scalep*y0i;
	Out[index_z++] = z0i;

	REAL x1r = Add[index_x++];    
	REAL y1r = InScale[index_y++];
	REAL z1r = x1r - scalep*y1r;
	Out[index_z++] = z1r;

	REAL x1i = Add[index_x++];    
	REAL y1i = InScale[index_y++];
	REAL z1i = x1i - scalep*y1i;
	Out[index_z++] = z1i;
    
	REAL x2r = Add[index_x++];
	REAL y2r = InScale[index_y++];     
	REAL z2r = x2r - scalep*y2r;
	Out[index_z++] = z2r;
   
	REAL x2i = Add[index_x++]; 
	REAL y2i = InScale[index_y++];
	REAL z2i = x2i - scalep*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 1 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r - scalep*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i - scalep*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r - scalep*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i - scalep*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r - scalep*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i - scalep*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 2 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r + scalep*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i + scalep*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r + scalep*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i + scalep*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r + scalep*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i + scalep*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 3 (AYPX)
	x0r = Add[index_x++];
	y0r = InScale[index_y++];
	z0r = x0r + scalep*y0r;
	Out[index_z++] = z0r;
  
	x0i = Add[index_x++];  
	y0i = InScale[index_y++];
	z0i = x0i + scalep*y0i;
	Out[index_z++] = z0i;

	x1r = Add[index_x++];    
	y1r = InScale[index_y++];
	z1r = x1r + scalep*y1r;
	Out[index_z++] = z1r;

	x1i = Add[index_x++];    
	y1i = InScale[index_y++];
	z1i = x1i + scalep*y1i;
	Out[index_z++] = z1i;
    
	x2r = Add[index_x++];
	y2r = InScale[index_y++];     
	z2r = x2r + scalep*y2r;
	Out[index_z++] = z2r;
   
	x2i = Add[index_x++]; 
	y2i = InScale[index_y++];
	z2i = x2i + scalep*y2i;  
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
	REAL scalep =  invd[N5-1].elem().elem().elem().elem();
	REAL *In    = (REAL *) &(y[N5-1].elem().elem(0).elem(0).real());
	REAL *Out   = (REAL *) &(chi[N5-1].elem(site).elem(0).elem(0).real());

	// (Vector) out = (Scalar) (*scalep) * (Vector) In
	//void scal_g5(REAL *Out, REAL *scalep, REAL *In, int n_4vec)
	register int index_x = 0;
	register int index_z = 0;

	// Spin Component 0
	REAL x0r = In[index_x++];
	REAL z0r = scalep*x0r;
	Out[index_z++] = z0r;
    
	REAL x0i = In[index_x++];
	REAL z0i = scalep*x0i;
	Out[index_z++] = z0i;
    
	REAL x1r = In[index_x++];
	REAL z1r = scalep*x1r;
	Out[index_z++] = z1r;
    
	REAL x1i = In[index_x++];
	REAL z1i = scalep*x1i;
	Out[index_z++] = z1i;
    
	REAL x2r = In[index_x++];     
	REAL z2r = scalep*x2r;
	Out[index_z++] = z2r;
    
	REAL x2i = In[index_x++];
	REAL z2i = scalep*x2i;
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] = z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] = z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] = z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] = z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] = z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
	Out[index_z++] = z2i;

	// Spin Component 2
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
	Out[index_z++] =-z2i;

	// Spin Component 3
	x0r = In[index_x++];
	z0r = scalep*x0r;
	Out[index_z++] =- z0r;
    
	x0i = In[index_x++];
	z0i = scalep*x0i;
	Out[index_z++] =- z0i;
    
	x1r = In[index_x++];
	z1r = scalep*x1r;
	Out[index_z++] =-z1r;
    
	x1i = In[index_x++];
	z1i = scalep*x1i;
	Out[index_z++] =-z1i;
    
	x2r = In[index_x++];     
	z2r = scalep*x2r;
	Out[index_z++] =-z2r;
    
	x2i = In[index_x++];
	z2i = scalep*x2i;
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
	REAL scalep   = invd[i].elem().elem().elem().elem();
	REAL scalep2  = u[i].elem().elem().elem().elem();
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
	REAL z0r = scalep*x0r ;
	z0r -= scalep2*y0r;
	Out[index_z++] = z0r;
    
	REAL x0i = InScale[index_x++];
	REAL y0i = Add[index_y++];
	REAL z0i = scalep*x0i;
	z0i -= scalep2*y0i;
	Out[index_z++] = z0i;
    
	REAL x1r = InScale[index_x++];
	REAL y1r = Add[index_y++];
	REAL z1r = scalep*x1r ;
	z1r -= scalep2*y1r;
	Out[index_z++] = z1r;
    
	REAL x1i = InScale[index_x++];
	REAL y1i = Add[index_y++];
	REAL z1i = scalep*x1i;
	z1i -= scalep2*y1i;
	Out[index_z++] = z1i;
    
	REAL x2r = InScale[index_x++];     
	REAL y2r = Add[index_y++];
	REAL z2r = scalep*x2r ;
	z2r -= scalep2*y2r;
	Out[index_z++] = z2r;
    
	REAL x2i = InScale[index_x++];
	REAL y2i = Add[index_y++];
	REAL z2i = scalep*x2i ;
	z2i -=  scalep2*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 1
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep*x0r;
	z0r -= scalep2*y0r;
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep*x0i;
	z0i -= scalep2*y0i;
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep*x1r ;
	z1r -= scalep2*y1r;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep*x1i;
	z1i -= scalep2*y1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep*x2r;
	z2r -=  scalep2*y2r;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i = scalep*x2i;
	z2i -= scalep2*y2i;  
	Out[index_z++] = z2i;

	// Spin Component 2 (AXPY3)
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep2*y0r;
	z0r -= scalep*x0r ;
   
	Out[index_z++] = z0r;
    
	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep2*y0i;
	z0i -= scalep*x0i;   
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep2*y1r;
	z1r -= scalep*x1r ;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep2*y1i;
	z1i -= scalep*x1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep2*y2r;
	z2r -= scalep*x2r ;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i =  scalep2*y2i;  
	z2i -= scalep*x2i ;
	Out[index_z++] = z2i;

	// Spin Component 3 (AXPY3)
	x0r = InScale[index_x++];
	y0r = Add[index_y++];
	z0r = scalep2*y0r;
	z0r -= scalep*x0r ;
    
	Out[index_z++] = z0r;

	x0i = InScale[index_x++];
	y0i = Add[index_y++];
	z0i = scalep2*y0i;
	z0i -= scalep*x0i;   
	Out[index_z++] = z0i;
    
	x1r = InScale[index_x++];
	y1r = Add[index_y++];
	z1r = scalep2*y1r;
	z1r -= scalep*x1r ;
	Out[index_z++] = z1r;
    
	x1i = InScale[index_x++];
	y1i = Add[index_y++];
	z1i = scalep2*y1i;
	z1i -= scalep*x1i;
	Out[index_z++] = z1i;
    
	x2r = InScale[index_x++];     
	y2r = Add[index_y++];
	z2r = scalep2*y2r;
	z2r -= scalep*x2r ;
	Out[index_z++] = z2r;
    
	x2i = InScale[index_x++];
	y2i = Add[index_y++];
	z2i =  scalep2*y2i;  
	z2i -= scalep*x2i ;
	Out[index_z++] = z2i;
      }
    }
  }

} // End Namespace Chroma

