/*! \file
 *  \brief Clover term linear operator
 *
 *  This particular implementation is specific to a scalar-like
 *  architecture
 */

#include "chromabase.h"
#include "actions/ferm/linop/clover_term_ssed.h"
#include "meas/glue/mesfield.h"


namespace Chroma 
{ 
  //extern void ssed_clover_apply(REAL64* diag, REAL64* offd, REAL64* psiptr, REAL64* chiptr, int n_sites);

  // Empty constructor. Must use create later
  SSEDCloverTerm::SSEDCloverTerm() {
    bool retry = false;

    // Allocate the arrays, with proper alignment

    try { 
      tri_diag = (PrimitiveClovDiag *)QDP::Allocator::theQDPAllocator::Instance().allocate( Layout::sitesOnNode()*sizeof(PrimitiveClovDiag) , QDP::Allocator::FAST );
    }
    catch( std::bad_alloc ) { 
      retry = true;
    }

    if( retry ) { 
      try { 
	tri_diag = (PrimitiveClovDiag *)QDP::Allocator::theQDPAllocator::Instance().allocate( Layout::sitesOnNode()*sizeof(PrimitiveClovDiag) , QDP::Allocator::DEFAULT );
      }
      catch( std::bad_alloc ) { 
	QDPIO::cerr << "Failed to allocate the tri_diag" << std::endl << std::flush ;
	QDP_abort(1);
      }
    }


    retry = false;
    try { 
      tri_off_diag = (PrimitiveClovOffDiag *)QDP::Allocator::theQDPAllocator::Instance().allocate( Layout::sitesOnNode()*sizeof(PrimitiveClovOffDiag) , QDP::Allocator::FAST );
    }
    catch( std::bad_alloc ) { 
      retry = true;
    }

    if( retry ) { 
      try { 
	tri_off_diag = (PrimitiveClovOffDiag *)QDP::Allocator::theQDPAllocator::Instance().allocate( Layout::sitesOnNode()*sizeof(PrimitiveClovOffDiag) , QDP::Allocator::DEFAULT );
      }
      catch( std::bad_alloc ) { 
	QDPIO::cerr << "Failed to allocate the tri_off_diag" << std::endl << std::flush ;
	QDP_abort(1);
      }
    }


  }


  SSEDCloverTerm::~SSEDCloverTerm() {
    if ( tri_diag != 0x0 ) { 
      QDP::Allocator::theQDPAllocator::Instance().free(tri_diag);
    }
    
    if ( tri_off_diag != 0x0 ) { 	
      QDP::Allocator::theQDPAllocator::Instance().free(tri_off_diag);
    }	
  }

  //! Creation routine
  void SSEDCloverTerm::create(Handle< FermState<T,P,Q> > fs,
			       const CloverFermActParams& param_)
  {
    START_CODE();
    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "SSEDCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }

    {
      Real ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= Real(0.5) * ff;
      param.clovCoeffT *= Real(0.5);
    }

    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    Real diag_mass;
    {
      Real ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }

    /* Calculate F(mu,nu) */
    multi1d<LatticeColorMatrix> f;
    mesField(f, u);
    makeClov(f, diag_mass);
    
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++) {
      choles_done[i] = false;
    }

#if 0
    // Testing code
    LatticeFermion foo;
    LatticeFermion res1=zero;
    LatticeFermion res2=zero;
    gaussian(foo);
    apply(res1,foo,PLUS, 0);
    applySite(res2,foo,PLUS, 0);

    LatticeFermion diff=res1-res2;
    XMLFileWriter fred("fred");
    push(fred, "stuff");
    write(fred, "diff", diff);
    pop(fred);

    QDPIO::cout << "sqrt( norm2( diff))= " << sqrt(norm2(diff)) << std::endl << std::flush;
    QDP_abort(1);
#endif

    END_CODE();
  }

    // Now copy
  void SSEDCloverTerm::create(Handle< FermState<T,P,Q> > fs,
			     const CloverFermActParams& param_,
			     const SSEDCloverTerm& from_)
  {
    START_CODE();

    const SSEDCloverTerm& from = dynamic_cast<const SSEDCloverTerm&>(from_);
    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "SSEDCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }
    
    {
      Real ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= Real(0.5) * ff;
      param.clovCoeffT *= Real(0.5);
    }
    
    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    Real diag_mass;
    {
      Real ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }
    
    
    /* Calculate F(mu,nu) */
    //multi1d<LatticeColorMatrix> f;
    //mesField(f, u);
    //makeClov(f, diag_mass);
    
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++) {
      choles_done[i] = from.choles_done[i];
    }
    
    tr_log_diag_ = from.tr_log_diag_;
    
    // Should be allocated by constructor.
    // tri_diag.resize(from.tri_diag.size());
    // tri_off_diag.resize(from.tri_off_diag.size());

    // This I claim is so performance uncritical, that
    // for now I won't parallelize it...
    for(int site=0; site < Layout::sitesOnNode(); site++) {
      for(int block = 0; block < 2; block++) { 
	for(int comp=0; comp< 6; comp++) { 
	  tri_diag[site][block][comp] = from.tri_diag[site][block][comp];
	}
	for(int comp=0; comp < 15; comp++) { 
	  tri_off_diag[site][block][comp] = from.tri_off_diag[site][block][comp];
	}
      }
    }
			    
    END_CODE();
  }


  /*
   * MAKCLOV 
   *
   *  In this routine, MAKCLOV calculates

   *    1 - (1/4)*sigma(mu,nu) F(mu,nu)

   *  using F from mesfield

   *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]

   *  using basis of SPPROD and stores in a lower triangular matrix
   *  (no diagonal) plus real diagonal

   *  where
   *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
   *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
   *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
   *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

   *  and

   *         | sigF(1)   sigF(3)     0         0     |
   *  sigF = | sigF(5)  -sigF(1)     0         0     |
   *         |   0         0      -sigF(0)  -sigF(2) |
   *         |   0         0      -sigF(4)   sigF(0) |
   *  where
   *    sigF(i) is a color matrix

   *  sigF(0) = i*(ClovT*E_z + ClovR*B_z)
   *          = i*(ClovT*F(3,2) + ClovR*F(1,0))
   *  sigF(1) = i*(ClovT*E_z - ClovR*B_z)
   *          = i*(ClovT*F(3,2) - ClovR*F(1,0))
   *  sigF(2) = i*(E_+ + B_+)
   *  sigF(3) = i*(E_+ - B_+)
   *  sigF(4) = i*(E_- + B_-)
   *  sigF(5) = i*(E_- - B_-)
   *  i*E_+ = (i*ClovT*E_x - ClovT*E_y)
   *        = (i*ClovT*F(3,0) - ClovT*F(3,1))
   *  i*E_- = (i*ClovT*E_x + ClovT*E_y)
   *        = (i*ClovT*F(3,0) + ClovT*F(3,1))
   *  i*B_+ = (i*ClovR*B_x - ClovR*B_y)
   *        = (i*ClovR*F(2,1) + ClovR*F(2,0))
   *  i*B_- = (i*ClovR*B_x + ClovR*B_y)
   *        = (i*ClovR*F(2,1) - ClovR*F(2,0))

   *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

   *  NOTE: the above definitions assume that the time direction, t_dir,
   *        is 3. In general F(k,j) is multiplied with ClovT if either
   *        k=t_dir or j=t_dir, and with ClovR otherwise.

   *+++
   *  Here are some notes on the origin of this routine. NOTE, ClovCoeff or u0
   *  are not actually used in MAKCLOV.
   *
   *  The clover mass term is suppose to act on a std::vector like
   *
   *  chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi

   *  Definitions used here (NOTE: no "i")
   *   sigma(mu,nu) = gamma(mu)*gamma(nu) - gamma(nu)*gamma(mu)
   *                = 2*gamma(mu)*gamma(nu)   for mu != nu
   *       
   *   chi = sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu < nu
   *       = (1/2) * sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu != nu
   *       = (1/4) * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)*psi
   *
   *
   * chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi
   *     = psi - (ClovCoeff/u0^3) * kappa * chi
   *     == psi - kappa * chi
   *
   *  We have absorbed ClovCoeff/u0^3 into kappa. A u0 was previously absorbed into kappa
   *  for compatibility to ancient conventions. 
   *---

   * Arguments:
   *  \param f         field strength tensor F(cb,mu,nu)        (Read)
   *  \param diag_mass effective mass term                      (Read)
   */

  namespace SSEDCloverEnv { 
    struct MakeClovArgs  {
      LatticeColorMatrix& f0;
      LatticeColorMatrix& f1;
      LatticeColorMatrix& f2;
      LatticeColorMatrix& f3;
      LatticeColorMatrix& f4;
      LatticeColorMatrix& f5;
      PrimitiveClovDiag* tri_diag;
      PrimitiveClovOffDiag* tri_off_diag;
      const Real& diag_mass;
    };

    inline
    void makeClovSiteLoop(int lo, int hi, int myId, MakeClovArgs* a)
    {
      LatticeColorMatrix& f0=a->f0;
      LatticeColorMatrix& f1=a->f1;
      LatticeColorMatrix& f2=a->f2;
      LatticeColorMatrix& f3=a->f3;
      LatticeColorMatrix& f4=a->f4;
      LatticeColorMatrix& f5=a->f5;
      PrimitiveClovDiag* tri_diag=a->tri_diag;
      PrimitiveClovOffDiag* tri_off_diag=a->tri_off_diag;
      const Real& diag_mass = a->diag_mass;

      for(int site = lo; site < hi; ++site) {

	for(int jj = 0; jj < 2; jj++) {

	  for(int ii = 0; ii < 2*Nc; ii++) {

	    tri_diag[site][jj][ii] = diag_mass.elem().elem().elem();
	  }
	}
	
	
	/* The appropriate clover coeffients are already included in the
	   field strengths F(mu,nu)! */
	
	RComplex<REAL> E_minus;
	RComplex<REAL> B_minus;
	RComplex<REAL> ctmp_0;
	RComplex<REAL> ctmp_1;
	RScalar<REAL> rtmp_0;
	RScalar<REAL> rtmp_1;
	
	for(int i = 0; i < Nc; ++i) {

	  /*# diag_L(i,0) = 1 - i*diag(E_z - B_z) */
	  /*#             = 1 - i*diag(F(3,2) - F(1,0)) */
	  ctmp_0 = f5.elem(site).elem().elem(i,i);
	  ctmp_0 -= f0.elem(site).elem().elem(i,i);
	  rtmp_0 = imag(ctmp_0);
	  tri_diag[site][0][i] += rtmp_0;
	  
	  /*# diag_L(i+Nc,0) = 1 + i*diag(E_z - B_z) */
	  /*#                = 1 + i*diag(F(3,2) - F(1,0)) */
	  tri_diag[site][0][i+Nc] -= rtmp_0;
	  
	  /*# diag_L(i,1) = 1 + i*diag(E_z + B_z) */
	  /*#             = 1 + i*diag(F(3,2) + F(1,0)) */
	  ctmp_1 = f5.elem(site).elem().elem(i,i);
	  ctmp_1 += f0.elem(site).elem().elem(i,i);
	  rtmp_1 = imag(ctmp_1);
	  tri_diag[site][1][i] -= rtmp_1;
	  
	  /*# diag_L(i+Nc,1) = 1 - i*diag(E_z + B_z) */
	  /*#                = 1 - i*diag(F(3,2) + F(1,0)) */
	  tri_diag[site][1][i+Nc] += rtmp_1;
	}
	
	/*# Construct lower triangular portion */
	/*# Block diagonal terms */
	for(int i = 1; i < Nc; ++i) {

	  for(int j = 0; j < i; ++j) {
	    
	    int elem_ij  = i*(i-1)/2 + j;
	    int elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;
	    
	    /*# L(i,j,0) = -i*(E_z - B_z)[i,j] */
	    /*#          = -i*(F(3,2) - F(1,0)) */
	    ctmp_0 = f0.elem(site).elem().elem(i,j);
	    ctmp_0 -= f5.elem(site).elem().elem(i,j);
	    tri_off_diag[site][0][elem_ij] = timesI(ctmp_0);
	    
	    /*# L(i+Nc,j+Nc,0) = +i*(E_z - B_z)[i,j] */
	    /*#                = +i*(F(3,2) - F(1,0)) */
	    tri_off_diag[site][0][elem_tmp] = -tri_off_diag[site][0][elem_ij];
	    
	    /*# L(i,j,1) = i*(E_z + B_z)[i,j] */
	    /*#          = i*(F(3,2) + F(1,0)) */
	    ctmp_1 = f5.elem(site).elem().elem(i,j);
	    ctmp_1 += f0.elem(site).elem().elem(i,j);
	    tri_off_diag[site][1][elem_ij] = timesI(ctmp_1);
	    
	    /*# L(i+Nc,j+Nc,1) = -i*(E_z + B_z)[i,j] */
	    /*#                = -i*(F(3,2) + F(1,0)) */
	    tri_off_diag[site][1][elem_tmp] = -tri_off_diag[site][1][elem_ij];
	  }
	}
      
	/*# Off-diagonal */
	for(int i = 0; i < Nc; ++i) {

	  for(int j = 0; j < Nc; ++j) {

	    // Flipped index
	    // by swapping i <-> j. In the past i would run slow
	    // and now j runs slow
	    int elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;
	    
	    /*# i*E_- = (i*E_x + E_y) */
	    /*#       = (i*F(3,0) + F(3,1)) */
	    E_minus = timesI(f2.elem(site).elem().elem(i,j));
	    E_minus += f4.elem(site).elem().elem(i,j);
	    
	    /*# i*B_- = (i*B_x + B_y) */
	    /*#       = (i*F(2,1) - F(2,0)) */
	    B_minus = timesI(f3.elem(site).elem().elem(i,j));
	    B_minus -= f1.elem(site).elem().elem(i,j);
	    
	    /*# L(i+Nc,j,0) = -i*(E_- - B_-)  */
	    tri_off_diag[site][0][elem_ij] = B_minus - E_minus;
	    
	    /*# L(i+Nc,j,1) = +i*(E_- + B_-)  */
	    tri_off_diag[site][1][elem_ij] = E_minus + B_minus;
	  }
	}
      } // End Site Loop
    } // End Function
  }; // End Namespace


  void SSEDCloverTerm::makeClov(const multi1d<LatticeColorMatrix>& f, const Real& diag_mass)
  {
    START_CODE();

    LatticeColorMatrix f0;
    LatticeColorMatrix f1;
    LatticeColorMatrix f2;
    LatticeColorMatrix f3;
    LatticeColorMatrix f4;
    LatticeColorMatrix f5;
    
    const int nodeSites = QDP::Layout::sitesOnNode();

    START_CODE();
  
    if ( Nd != 4 )
      QDP_error_exit("expecting Nd == 4", Nd);
  
    if ( Ns != 4 )
      QDP_error_exit("expecting Ns == 4", Ns);
  
    f0 = f[0] * getCloverCoeff(0,1);
    f1 = f[1] * getCloverCoeff(0,2);
    f2 = f[2] * getCloverCoeff(0,3);
    f3 = f[3] * getCloverCoeff(1,2);
    f4 = f[4] * getCloverCoeff(1,3);
    f5 = f[5] * getCloverCoeff(2,3);    

    SSEDCloverEnv::MakeClovArgs a={f0,f1,f2,f3,f4,f5,tri_diag,tri_off_diag, diag_mass};
    dispatch_to_threads(nodeSites, a, SSEDCloverEnv::makeClovSiteLoop);

    END_CODE();
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  void SSEDCloverTerm::choles(int cb)
  {
    START_CODE();

    // When you are doing the cholesky - also fill out the trace_log_diag piece)
    //chlclovms(tr_log_diag_, cb);
    ldagdlinv(tr_log_diag_,cb);
    
    END_CODE();
  }


  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   *
   * \return logarithm of the determinant  
   */
  Double SSEDCloverTerm::cholesDet(int cb) const
  {
    START_CODE();

    if( choles_done[cb] == false ) 
    {
      QDPIO::cout << "Error: YOu have not done the Cholesky.on this operator on this subset" << std::endl;
      QDPIO::cout << "You sure you shouldn't be asking invclov?" << std::endl;
      QDP_abort(1);
    }
LatticeReal ff=tr_log_diag_;

    if( param.sub_zero_usedP ) { 
 	QDPIO::cout << "Subtracting "<< param.sub_zero<<std::endl;
	LatticeReal tmp;
	tmp[rb[cb]] = param.sub_zero;
	ff[rb[cb]] -= tmp;
    }
    END_CODE();

    return sum(ff, rb[cb]);
  }    


  namespace SSEDCloverEnv { 
    struct LDagDLInvArgs { 
      LatticeReal& tr_log_diag;
      PrimitiveClovDiag* tri_diag;
      PrimitiveClovOffDiag* tri_off_diag;
      int cb;
    };

    inline 
    void lDagDLInvSiteLoop(int lo, int hi, int myId, LDagDLInvArgs *a)
    {
      LatticeReal& tr_log_diag = a->tr_log_diag;
      PrimitiveClovDiag* tri_diag = a->tri_diag;
      PrimitiveClovOffDiag* tri_off_diag = a->tri_off_diag;
      int cb = a->cb;

      RScalar<REAL> zip=0;
      
      int N = 2*Nc;

      // Loop through the sites.
      for(int ssite=lo; ssite < hi; ++ssite)  {

	int site = rb[cb].siteTable()[ssite];
	
	int site_neg_logdet=0;
	// Loop through the blocks on the site.
	for(int block=0; block < 2; block++) { 
	  
	  // Triangular storage 
	  // Big arrays to get good alignment...
	  
	  RScalar<REAL> inv_d[8] QDP_ALIGN16;
	  RComplex<REAL> inv_offd[16] QDP_ALIGN16;
	  RComplex<REAL> v[8] QDP_ALIGN16;
	  RScalar<REAL>  diag_g[8] QDP_ALIGN16;
	  // Algorithm 4.1.2 LDL^\dagger Decomposition
	  // From Golub, van Loan 3rd ed, page 139
	  for(int i=0; i < N; i++) { 
	    inv_d[i] = tri_diag[site][block][i];
	  }
	  for(int i=0; i < 15; i++) { 
	    inv_offd[i]  =tri_off_diag[site][block][i];
	  }
	  
	  for(int j=0; j < N; ++j) { 
	    
	    // Compute v(0:j-1)
	    //
	    // for i=0:j-2
	    //   v(i) = A(j,i) A(i,i)
	    // end
	    
	    
	    for(int i=0; i < j; i++) { 
	      int elem_ji = j*(j-1)/2 + i;
	      
	      RComplex<REAL> A_ii = cmplx( inv_d[i], zip );
	      v[i] = A_ii*adj(inv_offd[elem_ji]);
	    }
	    
	    // v(j) = A(j,j) - A(j, 0:j-2) v(0:j-2)
	    //                 ^ This is done with a loop over k ie:
	    //
	    // v(j) = A(j,j) - sum_k A*(j,k) v(k)     k=0...j-2
	    //
	    //      = A(j,j) - sum_k A*(j,k) A(j,k) A(k,k)
	    //      = A(j,j) - sum_k | A(j,k) |^2 A(k,k)
	    
	    v[j] = cmplx(inv_d[j],zip);
	    
	    for(int k=0; k < j; k++) { 
	      int elem_jk = j*(j-1)/2 + k;
	      v[j] -= inv_offd[elem_jk]*v[k];
	    }
	    
	    
	    // At this point in time v[j] has to be real, since
	    // A(j,j) is from diag ie real and all | A(j,k) |^2 is real
	    // as is A(k,k)
	    
	    // A(j,j) is the diagonal element - so store it.
	    inv_d[j] = real( v[j] );
	    
	    // Last line of algorithm:
	    // A( j+1 : n, j) = ( A(j+1:n, j) - A(j+1:n, 1:j-1)v(1:k-1) ) / v(j)
	    //
	    // use k as first colon notation and l as second so
	    // 
	    // for k=j+1 < n-1
	    //      A(k,j) = A(k,j) ;
	    //      for l=0 < j-1
	    //         A(k,j) -= A(k, l) v(l)
	    //      end
	    //      A(k,j) /= v(j);
	    //
	    for(int k=j+1; k < N; k++) { 
	      int elem_kj = k*(k-1)/2 + j;
	      for(int l=0; l < j; l++) { 
		int elem_kl = k*(k-1)/2 + l;
		inv_offd[elem_kj] -= inv_offd[elem_kl] * v[l];
	      }
	      inv_offd[elem_kj] /= v[j];
	    }
	  }
	  
	  // Now fix up the inverse
	  RScalar<REAL> one;
	  one.elem() = (REAL)1;
	  
	  for(int i=0; i < N; i++) { 
	    diag_g[i] = one/inv_d[i];
	    
	    // Compute the trace log
	    // NB we are always doing trace log | A | 
	    // (because we are always working with actually A^\dagger A
	    //  even in one flavour case where we square root)
	    tr_log_diag.elem(site).elem().elem().elem() += log(fabs(inv_d[i].elem()));
	    // However, it is worth counting just the no of negative logdets
	    // on site
	    if( inv_d[i].elem() < 0 ) { 
	      site_neg_logdet++;
	    }
	  }
	  // Now we need to invert the L D L^\dagger 
	  // We can do this by solving:
	  //
	  //  L D L^\dagger M^{-1} = 1   
	  //
	  // This can be done by solving L D X = 1  (X = L^\dagger M^{-1})
	  //
	  // Then solving L^\dagger M^{-1} = X
	  //
	  // LD is lower diagonal and so X will also be lower diagonal.
	  // LD X = 1 can be solved by forward substitution.
	  //
	  // Likewise L^\dagger is strictly upper triagonal and so
	  // L^\dagger M^{-1} = X can be solved by forward substitution.
	  RComplex<REAL> sum;
	  for(int k = 0; k < N; ++k) {

	    for(int i = 0; i < k; ++i) {
	      zero_rep(v[i]);
	    }
	    
	    /*# Forward substitution */
	    
	    // The first element is the inverse of the diagonal
	    v[k] = cmplx(diag_g[k],zip);
	    
	    for(int i = k+1; i < N; ++i) {
	      zero_rep(v[i]);
	      
	      for(int j = k; j < i; ++j) {
		int elem_ij = i*(i-1)/2+j;	
		
		// subtract l_ij*d_j*x_{kj}
		v[i] -= inv_offd[elem_ij] *inv_d[j]*v[j];
		
	      }
	
	      // scale out by 1/d_i
	      v[i] *= diag_g[i];
	    }
	    
	    /*# Backward substitution */
	    // V[N-1] remains unchanged
	    // Start from V[N-2]
	    
	    for(int i = N-2; (int)i >= (int)k; --i) {
	      for(int j = i+1; j < N; ++j) {
		int elem_ji = j*(j-1)/2 + i;
		// Subtract terms of typ (l_ji)*x_kj
		v[i] -= adj(inv_offd[elem_ji]) * v[j];
	      }
	    }
	    
	    /*# Overwrite column k of invcl.offd */
	    inv_d[k] = real(v[k]);
	    for(int i = k+1; i < N; ++i) {
	      
	      int elem_ik = i*(i-1)/2+k;
	      inv_offd[elem_ik] = v[i];
	    }
	  }
	  
	  
	  // Overwrite original data
	  for(int i=0; i < N; i++) { 
	    tri_diag[site][block][i] = inv_d[i];
	  }
	  for(int i=0; i < 15; i++) { 
	    tri_off_diag[site][block][i] = inv_offd[i];
	  }
	}
	
	if( site_neg_logdet != 0 ) { 
	  // Report if site had odd number of negative terms. (-ve def)
	  std::cout << "WARNING: Odd number of negative terms in Clover DET (" 
		    << site_neg_logdet<< ") at site: " << site << std::endl;
	}
      } // End Site Loop
    } // End Function
  } // End Namespace

   /*! An LDL^\dag decomposition and inversion? */
  void SSEDCloverTerm::ldagdlinv(LatticeReal& tr_log_diag, int cb)
  {
    START_CODE();

    if ( 2*Nc < 3 )
      QDP_error_exit("Matrix is too small", Nc, Ns);

    // Zero trace log
    tr_log_diag = zero;

    SSEDCloverEnv::LDagDLInvArgs a = { tr_log_diag, tri_diag, tri_off_diag, cb };
    dispatch_to_threads(rb[cb].numSiteTable(), a, SSEDCloverEnv::lDagDLInvSiteLoop);

    
    // This comes from the days when we used to do Cholesky
    choles_done[cb] = true;
    END_CODE();
  }


  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  extern void ssed_clover_apply(REAL64* diag, REAL64* offd, REAL64* psiptr, REAL64* chiptr, int n_sites);

  namespace SSEDCloverEnv{
    struct ApplyArgs { 
      LatticeFermion& chi;
      const LatticeFermion& psi;
      PrimitiveClovOffDiag* tri_off;
      PrimitiveClovDiag* tri_diag;
      int cb;
    };

    // Dispatch function for threading
    inline
    void orderedApplySiteLoop(int lo, int hi, int myID, 
			      ApplyArgs* a)
    {
      int n_4vec=hi-lo;
      int start=rb[ a->cb ].start()+lo;
      
      REAL64* chiptr = (REAL64 *)&( a->chi.elem(start).elem(0).elem(0).real());
      REAL64* psiptr = (REAL64 *)&( a->psi.elem(start).elem(0).elem(0).real());
      REAL64* offd = (REAL64 *)&(a->tri_off[start][0][0].real());
      REAL64* diag = (REAL64 *)&(a->tri_diag[start][0][0].elem());
      ssed_clover_apply(diag, offd, psiptr, chiptr, n_4vec);
    }

    inline 
    void unorderedApplySiteLoop(int lo, int hi, int myId, ApplyArgs* a)
    {
      int cb = a->cb;
      const multi1d<int>& tab = rb[cb].siteTable();


      for(int ssite=lo; ssite < hi; ++ssite) {
	
	int site = tab[ssite];
	unsigned long n_sites=1;
		
	REAL64* chiptr = (REAL64 *)&( a->chi.elem(site).elem(0).elem(0).real());
	REAL64* psiptr = (REAL64 *)&(a->psi.elem(site).elem(0).elem(0).real());
	REAL64* offd = (REAL64 *)&(a->tri_off[site][0][0].real());
	REAL64* diag = (REAL64 *)&(a->tri_diag[site][0][0].elem());
	ssed_clover_apply(diag, offd, psiptr, chiptr, n_sites);
	
      }
    }
  }

  void SSEDCloverTerm::apply(LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if ( Ns != 4 )
      QDP_error_exit("code requires Ns == 4", Ns);

    int n = 2*Nc;

    SSEDCloverEnv::ApplyArgs a = {chi, psi, tri_off_diag, tri_diag,cb};
    
    if( rb[cb].hasOrderedRep() ) {

      dispatch_to_threads(rb[cb].siteTable().size(), a, 
			  SSEDCloverEnv::orderedApplySiteLoop);

      
    }
    else {
      dispatch_to_threads(rb[cb].siteTable().size(), a, 
			  SSEDCloverEnv::unorderedApplySiteLoop);
    }

    getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  void SSEDCloverTerm::applySite(LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign, int site) const
  {
    START_CODE();

    if ( Ns != 4 )
      QDP_error_exit("code requires Ns == 4", Ns);

    int n = 2*Nc;

    RComplex<REAL>* cchi = (RComplex<REAL>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REAL>* ppsi = (const RComplex<REAL>*)&(psi.elem(site).elem(0).elem(0));


    cchi[ 0] = tri_diag[site][0][ 0]  * ppsi[ 0]
      +   conj(tri_off_diag[site][0][ 0]) * ppsi[ 1]
      +   conj(tri_off_diag[site][0][ 1]) * ppsi[ 2]
      +   conj(tri_off_diag[site][0][ 3]) * ppsi[ 3]
      +   conj(tri_off_diag[site][0][ 6]) * ppsi[ 4]
      +   conj(tri_off_diag[site][0][10]) * ppsi[ 5];
    
    cchi[ 1] = tri_diag[site][0][ 1]  * ppsi[ 1]
      +        tri_off_diag[site][0][ 0]  * ppsi[ 0]
      +   conj(tri_off_diag[site][0][ 2]) * ppsi[ 2]
      +   conj(tri_off_diag[site][0][ 4]) * ppsi[ 3]
      +   conj(tri_off_diag[site][0][ 7]) * ppsi[ 4]
      +   conj(tri_off_diag[site][0][11]) * ppsi[ 5];
    
    cchi[ 2] = tri_diag[site][0][ 2]  * ppsi[ 2]
      +        tri_off_diag[site][0][ 1]  * ppsi[ 0]
      +        tri_off_diag[site][0][ 2]  * ppsi[ 1]
      +   conj(tri_off_diag[site][0][ 5]) * ppsi[ 3]
      +   conj(tri_off_diag[site][0][ 8]) * ppsi[ 4]
      +   conj(tri_off_diag[site][0][12]) * ppsi[ 5];
    
    cchi[ 3] = tri_diag[site][0][ 3]  * ppsi[ 3]
      +        tri_off_diag[site][0][ 3]  * ppsi[ 0]
      +        tri_off_diag[site][0][ 4]  * ppsi[ 1]
      +        tri_off_diag[site][0][ 5]  * ppsi[ 2]
      +   conj(tri_off_diag[site][0][ 9]) * ppsi[ 4]
      +   conj(tri_off_diag[site][0][13]) * ppsi[ 5];
    
    cchi[ 4] = tri_diag[site][0][ 4]  * ppsi[ 4]
      +        tri_off_diag[site][0][ 6]  * ppsi[ 0]
      +        tri_off_diag[site][0][ 7]  * ppsi[ 1]
      +        tri_off_diag[site][0][ 8]  * ppsi[ 2]
      +        tri_off_diag[site][0][ 9]  * ppsi[ 3]
      +   conj(tri_off_diag[site][0][14]) * ppsi[ 5];
    
    cchi[ 5] = tri_diag[site][0][ 5]  * ppsi[ 5]
      +        tri_off_diag[site][0][10]  * ppsi[ 0]
      +        tri_off_diag[site][0][11]  * ppsi[ 1]
      +        tri_off_diag[site][0][12]  * ppsi[ 2]
      +        tri_off_diag[site][0][13]  * ppsi[ 3]
      +        tri_off_diag[site][0][14]  * ppsi[ 4];
    
    cchi[ 6] = tri_diag[site][1][ 0]  * ppsi[ 6]
      +   conj(tri_off_diag[site][1][ 0]) * ppsi[ 7]
      +   conj(tri_off_diag[site][1][ 1]) * ppsi[ 8]
      +   conj(tri_off_diag[site][1][ 3]) * ppsi[ 9]
      +   conj(tri_off_diag[site][1][ 6]) * ppsi[10]
      +   conj(tri_off_diag[site][1][10]) * ppsi[11];
    
    cchi[ 7] = tri_diag[site][1][ 1]  * ppsi[ 7]
      +        tri_off_diag[site][1][ 0]  * ppsi[ 6]
      +   conj(tri_off_diag[site][1][ 2]) * ppsi[ 8]
      +   conj(tri_off_diag[site][1][ 4]) * ppsi[ 9]
      +   conj(tri_off_diag[site][1][ 7]) * ppsi[10]
      +   conj(tri_off_diag[site][1][11]) * ppsi[11];
    
    cchi[ 8] = tri_diag[site][1][ 2]  * ppsi[ 8]
      +        tri_off_diag[site][1][ 1]  * ppsi[ 6]
      +        tri_off_diag[site][1][ 2]  * ppsi[ 7]
      +   conj(tri_off_diag[site][1][ 5]) * ppsi[ 9]
      +   conj(tri_off_diag[site][1][ 8]) * ppsi[10]
      +   conj(tri_off_diag[site][1][12]) * ppsi[11];
    
    cchi[ 9] = tri_diag[site][1][ 3]  * ppsi[ 9]
      +        tri_off_diag[site][1][ 3]  * ppsi[ 6]
      +        tri_off_diag[site][1][ 4]  * ppsi[ 7]
      +        tri_off_diag[site][1][ 5]  * ppsi[ 8]
      +   conj(tri_off_diag[site][1][ 9]) * ppsi[10]
      +   conj(tri_off_diag[site][1][13]) * ppsi[11];
    
    cchi[10] = tri_diag[site][1][ 4]  * ppsi[10]
      +        tri_off_diag[site][1][ 6]  * ppsi[ 6]
      +        tri_off_diag[site][1][ 7]  * ppsi[ 7]
      +        tri_off_diag[site][1][ 8]  * ppsi[ 8]
      +        tri_off_diag[site][1][ 9]  * ppsi[ 9]
      +   conj(tri_off_diag[site][1][14]) * ppsi[11];
    
    cchi[11] = tri_diag[site][1][ 5]  * ppsi[11]
      +        tri_off_diag[site][1][10]  * ppsi[ 6]
      +        tri_off_diag[site][1][11]  * ppsi[ 7]
      +        tri_off_diag[site][1][12]  * ppsi[ 8]
      +        tri_off_diag[site][1][13]  * ppsi[ 9]
      +        tri_off_diag[site][1][14]  * ppsi[10];


    END_CODE();
  }




  //! TRIACNTR 
  /*! 
   * \ingroup linop
   *
   *  Calculates
   *     Tr_D ( Gamma_mat L )
   *
   * This routine is specific to Wilson fermions!
   * 
   *  the trace over the Dirac indices for one of the 16 Gamma matrices
   *  and a hermitian color x spin matrix A, stored as a block diagonal
   *  complex lower triangular matrix L and a real diagonal diag_L.

   *  Here 0 <= mat <= 15 and
   *  if mat = mat_1 + mat_2 * 2 + mat_3 * 4 + mat_4 * 8
   *
   *  Gamma(mat) = gamma(1)^(mat_1) * gamma(2)^(mat_2) * gamma(3)^(mat_3)
   *             * gamma(4)^(mat_4)
   *
   *  Further, in basis for the Gamma matrices used, A is of the form
   *
   *      | A_0 |  0  |
   *  A = | --------- |
   *      |  0  | A_1 |
   *
   *
   * Arguments:
   *
   *  \param B         the resulting SU(N) color matrix	  (Write) 
   *  \param clov      clover term                        (Read) 
   *  \param mat       label of the Gamma matrix          (Read)
   */

  namespace SSEDCloverEnv {
    struct TriaCntrArgs { 
      LatticeColorMatrix& B;
      PrimitiveClovDiag* tri_diag;
      PrimitiveClovOffDiag* tri_off;
      int mat;
      int cb;
    };
  

    inline 
    void triaCntrSiteLoop(int lo, int hi, int myId, TriaCntrArgs* a) 
    {
      PrimitiveClovDiag* tri_diag = a->tri_diag;
      PrimitiveClovOffDiag* tri_off_diag = a->tri_off;
      LatticeColorMatrix& B = a->B;
      int mat = a->mat;
      int cb = a->cb;


      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) {
      
	int site = rb[cb].siteTable()[ssite];
	
	switch( mat ) {

	case 0:
	  /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
	  /*#               0  1  0  0 */
	  /*#               0  0  1  0 */
	  /*#               0  0  0  1 */
	  /*# From diagonal part */
	  {
	    RComplex<REAL> lctmp0;
	    RScalar<REAL> lr_zero0;
	    RScalar<REAL> lrtmp0;
	    
	    lr_zero0 = 0;
	    
	    for(int i0 = 0; i0 < Nc; ++i0) {
	      lrtmp0 = tri_diag[site][0][i0];
	      lrtmp0 += tri_diag[site][0][i0+Nc];
	      lrtmp0 += tri_diag[site][1][i0];
	      lrtmp0 += tri_diag[site][1][i0+Nc];
	      B.elem(site).elem().elem(i0,i0) = cmplx(lrtmp0,lr_zero0);
	    }
	    
	    /*# From lower triangular portion */
	    int elem_ij0 = 0;
	    for(int i0 = 1; i0 < Nc; ++i0) {
	      
	      int elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc;
	      
	      for(int j0 = 0; j0 < i0; ++j0) {
		
		lctmp0 = tri_off_diag[site][0][elem_ij0];
		lctmp0 += tri_off_diag[site][0][elem_ijb0];
		lctmp0 += tri_off_diag[site][1][elem_ij0];
		lctmp0 += tri_off_diag[site][1][elem_ijb0];
		
		B.elem(site).elem().elem(j0,i0) = lctmp0;
		B.elem(site).elem().elem(i0,j0) = adj(lctmp0);
		
		
		elem_ij0++;
		elem_ijb0++;
	      }
	    }
	  }
	  break;
	  
	case 3:
	  /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */
	  /*#               0  i  0  0 */
	  /*#               0  0 -i  0 */
	  /*#               0  0  0  i */
	  /*# From diagonal part */
	  {
	    
	    RComplex<REAL> lctmp3;
	    RScalar<REAL> lr_zero3;
	    RScalar<REAL> lrtmp3;
	    
	    lr_zero3 = 0;
	    
	    for(int i3 = 0; i3 < Nc; ++i3) {
	      
	      lrtmp3 = tri_diag[site][0][i3+Nc];
	      lrtmp3 -= tri_diag[site][0][i3];
	      lrtmp3 -= tri_diag[site][1][i3];
	      lrtmp3 += tri_diag[site][1][i3+Nc];
	      B.elem(site).elem().elem(i3,i3) = cmplx(lr_zero3,lrtmp3);
	    }
	    
	    /*# From lower triangular portion */
	    int elem_ij3 = 0;
	    for(int i3 = 1; i3 < Nc; ++i3) {
	      
	      int elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc;
	      
	      for(int j3 = 0; j3 < i3; ++j3) {
		
		lctmp3 = tri_off_diag[site][0][elem_ijb3];
		lctmp3 -= tri_off_diag[site][0][elem_ij3];
		lctmp3 -= tri_off_diag[site][1][elem_ij3];
		lctmp3 += tri_off_diag[site][1][elem_ijb3];
		
		B.elem(site).elem().elem(j3,i3) = timesI(adj(lctmp3));
		B.elem(site).elem().elem(i3,j3) = timesI(lctmp3);
		
		elem_ij3++;
		elem_ijb3++;
	      }
	    }
	  }
	  break;
	  
	case 5:
	  /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */
	  /*#               1  0  0  0 */
	  /*#               0  0  0 -1 */
	  /*#               0  0  1  0 */
	  {
	    
	    
	    RComplex<REAL> lctmp5;
	    RScalar<REAL> lrtmp5;
	    
	    for(int i5 = 0; i5 < Nc; ++i5) {
	      
	      int elem_ij5 = (i5+Nc)*(i5+Nc-1)/2;
	      
	      for(int j5 = 0; j5 < Nc; ++j5) {
		
		int elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5;
		
		
		lctmp5 = adj(tri_off_diag[site][0][elem_ji5]);
		lctmp5 -= tri_off_diag[site][0][elem_ij5];
		lctmp5 += adj(tri_off_diag[site][1][elem_ji5]);
		lctmp5 -= tri_off_diag[site][1][elem_ij5];
		
		
		B.elem(site).elem().elem(i5,j5) = lctmp5;
		
		elem_ij5++;
	      }
	    }
	  }
	  break;
	  
	case 6:
	  /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */
	  /*#              -i  0  0  0 */
	  /*#               0  0  0 -i */
	  /*#               0  0 -i  0 */
	  {
	    
	    RComplex<REAL> lctmp6;
	    RScalar<REAL> lrtmp6;
	    
	    for(int i6 = 0; i6 < Nc; ++i6) {
	      
	      int elem_ij6 = (i6+Nc)*(i6+Nc-1)/2;
	      
	      for(int j6 = 0; j6 < Nc; ++j6) {
		
		int elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6;
		
		lctmp6 = adj(tri_off_diag[site][0][elem_ji6]);
		lctmp6 += tri_off_diag[site][0][elem_ij6];
		lctmp6 += adj(tri_off_diag[site][1][elem_ji6]);
		lctmp6 += tri_off_diag[site][1][elem_ij6];
		
		B.elem(site).elem().elem(i6,j6) = timesMinusI(lctmp6);
		
		elem_ij6++;
	      }
	    }
	  }
	  break;
	  
	case 9:
	  /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */
	  /*#               i  0  0  0 */
	  /*#               0  0  0 -i */
	  /*#               0  0 -i  0 */
	  {
	    RComplex<REAL> lctmp9;
	    RScalar<REAL> lrtmp9;
	    
	    for(int i9 = 0; i9 < Nc; ++i9) {
	      
	      int elem_ij9 = (i9+Nc)*(i9+Nc-1)/2;
	      
	      for(int j9 = 0; j9 < Nc; ++j9) {
		
		int elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9;
		
		lctmp9 = adj(tri_off_diag[site][0][elem_ji9]);
		lctmp9 += tri_off_diag[site][0][elem_ij9];
		lctmp9 -= adj(tri_off_diag[site][1][elem_ji9]);
		lctmp9 -= tri_off_diag[site][1][elem_ij9];
		
		B.elem(site).elem().elem(i9,j9) = timesI(lctmp9);
		
		elem_ij9++;
	      }
	    }
	  }
	  break;
	  
	case 10:
	/*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */
	/*#               1  0  0  0 */
	/*#               0  0  0  1 */
	/*#               0  0 -1  0 */
	  {
	    
	    RComplex<REAL> lctmp10;
	    RScalar<REAL> lrtmp10;
	    
	    for(int i10 = 0; i10 < Nc; ++i10) {
	      
	      int elem_ij10 = (i10+Nc)*(i10+Nc-1)/2; 
		
	      for(int j10 = 0; j10 < Nc; ++j10) {
		
		int elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10;
		
		lctmp10 = adj(tri_off_diag[site][0][elem_ji10]);
		lctmp10 -= tri_off_diag[site][0][elem_ij10];
		lctmp10 -= adj(tri_off_diag[site][1][elem_ji10]);
		lctmp10 += tri_off_diag[site][1][elem_ij10];
		
		B.elem(site).elem().elem(i10,j10) = lctmp10;
		
		elem_ij10++;
	      }
	    }
	  }
	  
	  break;
	  
	case 12:
	  /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */
	  /*#               0 -i  0  0 */
	  /*#               0  0 -i  0 */
	  /*#               0  0  0  i */
	  /*# From diagonal part */
	  {
	  
	  
	    RComplex<REAL> lctmp12;
	    RScalar<REAL> lr_zero12;
	    RScalar<REAL> lrtmp12;
	    
	    lr_zero12 = 0;
	    
	    for(int i12 = 0; i12 < Nc; ++i12) {
	      
	      lrtmp12 = tri_diag[site][0][i12];
	      lrtmp12 -= tri_diag[site][0][i12+Nc];
	      lrtmp12 -= tri_diag[site][1][i12];
	      lrtmp12 += tri_diag[site][1][i12+Nc];
	      B.elem(site).elem().elem(i12,i12) = cmplx(lr_zero12,lrtmp12);
	    }
	    
	    /*# From lower triangular portion */
	    int elem_ij12 = 0;
	    for(int i12 = 1; i12 < Nc; ++i12) {
	      
	      int elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc;
	      
	      for(int j12 = 0; j12 < i12; ++j12) {
		
		lctmp12 = tri_off_diag[site][0][elem_ij12];
		lctmp12 -= tri_off_diag[site][0][elem_ijb12];
		lctmp12 -= tri_off_diag[site][1][elem_ij12];
		lctmp12 += tri_off_diag[site][1][elem_ijb12];
		
		B.elem(site).elem().elem(i12,j12) = timesI(lctmp12);
		B.elem(site).elem().elem(j12,i12) = timesI(adj(lctmp12));
		
		elem_ij12++;
		elem_ijb12++;
	      }
	    }
	  }
	  break;
	
	default:
	  {
	    B = zero;
	    QDPIO::cout << "BAD DEFAULT CASE HIT" << std::endl;
	  }
	} // End Switch
	
      } // End Site Loop.
    } // End Function
  } // End Namespace 

  void SSEDCloverTerm::triacntr(LatticeColorMatrix& B, int mat, int cb) const
  {
    START_CODE();

    B = zero;

    if ( mat < 0  ||  mat > 15 )
    {
      QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << std::endl;
      QDP_abort(1);
    }
    
    SSEDCloverEnv::TriaCntrArgs a = { B, tri_diag, tri_off_diag, mat, cb};
    dispatch_to_threads(rb[cb].numSiteTable(), a, 
			SSEDCloverEnv::triaCntrSiteLoop);
  

    END_CODE();
  }

  //! Returns the appropriate clover coefficient for indices mu and nu
  Real SSEDCloverTerm::getCloverCoeff(int mu, int nu) const 
  { 
    START_CODE();

    if( param.anisoParam.anisoP ) 
    { 
      if (mu==param.anisoParam.t_dir || nu == param.anisoParam.t_dir) { 
	return param.clovCoeffT;
      }
      else { 
	// Otherwise return the spatial coeff
	return param.clovCoeffR;
      }
    }
    else { 
      // If there is no anisotropy just return the spatial one, it will
      // be the same as the temporal one
      return param.clovCoeffR; 
    } 
    
    END_CODE();
  }
}


 

