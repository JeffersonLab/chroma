// -*- C++ -*-
// $Id: clover_term_qdp_w.h,v 3.12 2009-10-06 20:35:48 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_jit_w_h__
#define __clover_term_jit_w_h__

#include "state.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/clover_term_base_w.h"
#include "meas/glue/mesfield.h"
namespace Chroma 
{ 

#if 0
  //! Special structure used for triangular objects
  template<typename R>
  struct PrimitiveClovTriang
  {
    RScalar<R>   diag[2][2*Nc];
    RComplex<R>  offd[2][2*Nc*Nc-Nc];
  };

  template<typename R>
  struct QUDAPackedClovSite {
    R diag1[6];
    R offDiag1[15][2];
    R diag2[6];
    R offDiag2[15][2];
  };
#endif

  // Reader/writers
  /*! \ingroup linop */
#if 0
  void read(XMLReader& xml, const string& path, PrimitiveClovTriang& param);

  /*! \ingroup linop */
  void write(XMLWriter& xml, const string& path, const PrimitiveClovTriang& param);
#endif

  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  template<typename T, typename U>
  class JitCloverTermT : public CloverTermBase<T, U>
  {
  public:
    // Typedefs to save typing
    typedef typename WordType<T>::Type_t REALT;

    typedef OLattice< PScalar< PScalar< RScalar< typename WordType<T>::Type_t> > > > LatticeREAL;
    typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;

    //! Constructor.
    JitCloverTermT();

    //! Real need for cleanup here
    ~JitCloverTermT() {
      if (Layout::primaryNode())
	QDP_info_primary("JIT Clover: Signing off triangular field tri_id=%u",(unsigned)tri_id);
      QDPCache::Instance().signoff( tri_id );
    }

    //! Creation routine
    void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
		const CloverFermActParams& param_);

    virtual void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
			const CloverFermActParams& param_,
			const JitCloverTermT<T,U>& from_);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     */
    void choles(int cb);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     * \return logarithm of the determinant  
     */
    Double cholesDet(int cb) const ;

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
     * \param cb      Checkerboard of OUTPUT vector               (Read) 
     */
    void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const;


    void applySite(T& chi, const T& psi, enum PlusMinus isign, int site) const;

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(U& B, int mat, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T, multi1d<U>, multi1d<U> >& getFermBC() const {return *fbc;}

    //! PACK UP the Clover term for QUDA library:
    void packForQUDA(multi1d<QUDAPackedClovSite<REALT> >& quda_pack, int cb) const; 


      
  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<U>& f, const RealT& diag_mass);

    //! Invert the clover term on cb
    void chlclovms(LatticeREAL& log_diag, int cb);
    void ldagdlinv(LatticeREAL& tr_log_diag, int cb);

    //! Get the u field
    const multi1d<U>& getU() const {return u;}

    //! Calculates Tr_D ( Gamma_mat L )
    Real getCloverCoeff(int mu, int nu) const;

  private:
    Handle< FermBC<T,multi1d<U>,multi1d<U> > >      fbc;
    multi1d<U>  u;
    CloverFermActParams          param;
    LatticeREAL                  tr_log_diag_; // Fill this out during create
                                                  // but save the global sum until needed.
    multi1d<bool> choles_done;   // Keep note of whether the decomposition has been done
                                 // on a particular checkerboard. 

    //multi1d<PrimitiveClovTriang<REALT> >  tri;

    int tri_id;
    
  };


  // Constructor.
  template<typename T, typename U>
  JitCloverTermT<T,U>::JitCloverTermT() 
  {
    const int nodeSites = QDP::Layout::sitesOnNode();
    size_t bytes_total = sizeof( PrimitiveClovTriang< typename JitCloverTermT::REALT > ) * nodeSites;
    tri_id = QDPCache::Instance().registrate( bytes_total , 1 );
    if (Layout::primaryNode())
      QDP_info_primary("JIT Clover: Registering triangular field for clover %u  tri_id=%u",(unsigned)bytes_total,(unsigned)tri_id);
  }


  // Now copy
  template<typename T, typename U>
  void JitCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_,
				   const JitCloverTermT<T,U>& from)
  {
    START_CODE();
    u.resize(Nd);

    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "QDPCloverTerm: error: fbc is null" << endl;
      QDP_abort(1);
    }
    
    {
      RealT ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= Real(0.5) * ff;
      param.clovCoeffT *= Real(0.5);
    }
    
    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
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

    const int nodeSites = QDP::Layout::sitesOnNode();
    size_t bytes_total = sizeof( PrimitiveClovTriang< typename JitCloverTermT::REALT > ) * nodeSites;

    CudaMemcpy( QDPCache::Instance().getDevicePtr( tri_id ),
	        QDPCache::Instance().getDevicePtr( from.tri_id ), 
	        bytes_total );

    END_CODE();  
  }


  //! Creation routine
  template<typename T, typename U>
  void JitCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_)
  {
    START_CODE();
   
    u.resize(Nd);
    
    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "QDPCloverTerm: error: fbc is null" << endl;
      QDP_abort(1);
    }

    {
      RealT ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= RealT(0.5) * ff;
      param.clovCoeffT *= RealT(0.5);
    }
    
    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }
    
    /* Calculate F(mu,nu) */
    multi1d<U> f;
    mesField(f, u);
    makeClov(f, diag_mass);
    
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++) {
      choles_done[i] = false;
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
   *  The clover mass term is suppose to act on a vector like
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

  /*Threaded this. It needs a QMT arg struct and I've extracted the site loop */
  namespace JITCloverEnv { 
    
    template<typename U>
    struct JITCloverMakeClovArg {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;
      const RealT& diag_mass;
      const U& f0;
      const U& f1;
      const U& f2;
      const U& f3;
      const U& f4;
      const U& f5;
      int tri_id;
    };


    template<typename U>
    inline 
    void makeClovJIT(int nodeSites, JITCloverMakeClovArg<U> *a)
    {
      typedef typename JITCloverMakeClovArg<U>::RealT RealT;
      typedef typename JITCloverMakeClovArg<U>::REALT REALT;
      
      const RealT& diag_mass = a->diag_mass;
      const U& f0=a->f0;
      const U& f1=a->f1;
      const U& f2=a->f2;
      const U& f3=a->f3;
      const U& f4=a->f4;
      const U& f5=a->f5;
      int tri_id=a->tri_id;

      while (1) {

	static QDPJit::SharedLibEntry sharedLibEntry;
	static MapVolumes*              mapVolumes;
	static string                   strId;
	static string                   prg;

	QDPJitArgs cudaArgs;
	string typeU,strREALT,codeDiagMass;

	int argDestPtr = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( tri_id ) );
	int argNum = cudaArgs.addInt( nodeSites );

	QDP_debug("makeClov dev!");

	if (!mapVolumes) {

	  string codeF0;
	  string codeF1;
	  string codeF2;
	  string codeF3;
	  string codeF4;
	  string codeF5;

	  getTypeString( typeU , f1 , cudaArgs );
	  getTypeString( strREALT , REALT(0) );

	  if (!getCodeString( codeDiagMass , diag_mass , "idx", cudaArgs )) {QDP_info(": could not cache diag_mass");break;}
	  if (!getCodeString( codeF0 , f0 , "idx", cudaArgs )) { QDP_info(": could not cache f0"); break;  }
	  if (!getCodeString( codeF1 , f1 , "idx", cudaArgs )) { QDP_info(": could not cache f1"); break;  }
	  if (!getCodeString( codeF2 , f2 , "idx", cudaArgs )) { QDP_info(": could not cache f2"); break;  }
	  if (!getCodeString( codeF3 , f3 , "idx", cudaArgs )) { QDP_info(": could not cache f3"); break;  }
	  if (!getCodeString( codeF4 , f4 , "idx", cudaArgs )) { QDP_info(": could not cache f4"); break;  }
	  if (!getCodeString( codeF5 , f5 , "idx", cudaArgs )) { QDP_info(": could not cache f5"); break;  }	  

	  ostringstream osId;
	  osId << "makeClov "  << typeU;
	  strId = osId.str();
	  xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	  cout << "strId = " << strId << endl;
#endif

	  std::ostringstream sTri;
	  sTri << "((PrimitiveClovTriang*)( " << cudaArgs.getCode(argDestPtr) << " ))[idx]";
	  string strTri = sTri.str();
      
	  std::ostringstream sprg;

	  sprg << " typedef " << strREALT << " REALT;" << endl;
	  sprg << " const int Nc = " << Nc << ";" << endl;

	  sprg << "	struct PrimitiveClovTriang" << endl;
	  sprg << "	{" << endl;
	  sprg << "   RScalar<REALT>   diag[2][2*Nc];" << endl;
	  sprg << "   RComplex<REALT>  offd[2][2*Nc*Nc-Nc];" << endl;
	  sprg << "	};" << endl;

	  sprg << "int idx = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	  sprg << "if (idx < " << cudaArgs.getCode(argNum) << " ) {" << endl;

	  sprg << "for(int jj = 0; jj < 2; jj++) {" << endl;
	  sprg << "  for(int ii = 0; ii < 2*Nc; ii++) {" << endl;
	  sprg << "    " << strTri << ".diag[jj][ii] = (" << codeDiagMass << ").elem().elem();" << endl;
	  sprg << "  }" << endl;
	  sprg << "}" << endl;

	  sprg << "	RComplex<REALT> E_minus;" << endl;
	  sprg << "	RComplex<REALT> B_minus;" << endl;
	  sprg << "	RComplex<REALT> ctmp_0;" << endl;
	  sprg << "	RComplex<REALT> ctmp_1;" << endl;
	  sprg << "	RScalar<REALT> rtmp_0;" << endl;
	  sprg << "	RScalar<REALT> rtmp_1;" << endl;

	  sprg << "	for(int i = 0; i < Nc; ++i) {" << endl;
	  sprg << "	  ctmp_0 = " << codeF5 << ".elem().elem(i,i);" << endl;
	  sprg << "	  ctmp_0 -= " << codeF0 << ".elem().elem(i,i);" << endl;
	  sprg << "	  rtmp_0 = imag(ctmp_0);" << endl;
	  sprg << "	  " << strTri << ".diag[0][i] += rtmp_0;" << endl;
	  sprg << "	  " << strTri << ".diag[0][i+Nc] -= rtmp_0;" << endl;
	  sprg << "	  ctmp_1 = " << codeF5 << ".elem().elem(i,i);" << endl;
	  sprg << "	  ctmp_1 += " << codeF0 << ".elem().elem(i,i);" << endl;
	  sprg << "	  rtmp_1 = imag(ctmp_1);" << endl;
	  sprg << "	  " << strTri << ".diag[1][i] -= rtmp_1;" << endl;
	  sprg << "	  " << strTri << ".diag[1][i+Nc] += rtmp_1;" << endl;
	  sprg << "	}" << endl;
	  sprg << "	for(int i = 1; i < Nc; ++i) {" << endl;
	  sprg << "	  for(int j = 0; j < i; ++j) {" << endl;
	  sprg << "	    int elem_ij  = i*(i-1)/2 + j;" << endl;
	  sprg << "	    int elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;" << endl;
	  sprg << "	    ctmp_0 = " << codeF0 << ".elem().elem(i,j);" << endl;
	  sprg << "	    ctmp_0 -= " << codeF5 << ".elem().elem(i,j);" << endl;
	  sprg << "	    " << strTri << ".offd[0][elem_ij] = timesI(ctmp_0);" << endl;
	  sprg << "	    " << strTri << ".offd[0][elem_tmp] = -" << strTri << ".offd[0][elem_ij];" << endl;
	  sprg << "	    ctmp_1 = " << codeF5 << ".elem().elem(i,j);" << endl;
	  sprg << "	    ctmp_1 += " << codeF0 << ".elem().elem(i,j);" << endl;
	  sprg << "	    " << strTri << ".offd[1][elem_ij] = timesI(ctmp_1);" << endl;
	  sprg << "	    " << strTri << ".offd[1][elem_tmp] = -" << strTri << ".offd[1][elem_ij];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	}" << endl;
	  sprg << "	for(int i = 0; i < Nc; ++i) {" << endl;
	  sprg << "	  for(int j = 0; j < Nc; ++j) {" << endl;
	  sprg << "	    int elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;" << endl;
	  sprg << "	    E_minus = timesI(" << codeF2 << ".elem().elem(i,j));" << endl;
	  sprg << "	    E_minus += " << codeF4 << ".elem().elem(i,j);" << endl;
	  sprg << "	    B_minus = timesI(" << codeF3 << ".elem().elem(i,j));" << endl;
	  sprg << "	    B_minus -= " << codeF1 << ".elem().elem(i,j);" << endl;
	  sprg << "	    " << strTri << ".offd[0][elem_ij] = B_minus - E_minus;" << endl;
	  sprg << "	    " << strTri << ".offd[1][elem_ij] = E_minus + B_minus;" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	}" << endl;
	  sprg << "}" << endl;
	  prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	  cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	} else {
	  if (!cacheLock( diag_mass  , cudaArgs )) { QDP_info("eval: could not cache diag_mass");  break;   }
	  if (!cacheLock( f0 , cudaArgs )) { QDP_info("eval: could not cache f0");  break;   }
	  if (!cacheLock( f1 , cudaArgs )) { QDP_info("eval: could not cache f1");  break;   }
	  if (!cacheLock( f2 , cudaArgs )) { QDP_info("eval: could not cache f2");  break;   }
	  if (!cacheLock( f3 , cudaArgs )) { QDP_info("eval: could not cache f3");  break;   }
	  if (!cacheLock( f4 , cudaArgs )) { QDP_info("eval: could not cache f4");  break;   }
	  if (!cacheLock( f5 , cudaArgs )) { QDP_info("eval: could not cache f5");  break;   }
	}

	if (!QDPJit::Instance()( strId , prg , cudaArgs.getDevPtr() , nodeSites , sharedLibEntry  , mapVolumes )) {
	  QDP_info("makeClov call to cuda jitter failed");
	  break;
	}

	return;
      }
      QDP_error("makeClov host not implemented");
    }
  }

  
  /* This now just sets up and dispatches... */
  template<typename T, typename U>
  void JitCloverTermT<T,U>::makeClov(const multi1d<U>& f, const RealT& diag_mass)
  {
    START_CODE();
    
    if ( Nd != 4 ){
      QDPIO::cerr << __func__ << ": expecting Nd==4" << endl;
      QDP_abort(1);
    }
    
    if ( Ns != 4 ){
      QDPIO::cerr << __func__ << ": expecting Ns==4" << endl;
      QDP_abort(1);
    }
  
    U f0 = f[0] * getCloverCoeff(0,1);
    U f1 = f[1] * getCloverCoeff(0,2);
    U f2 = f[2] * getCloverCoeff(0,3);
    U f3 = f[3] * getCloverCoeff(1,2);
    U f4 = f[4] * getCloverCoeff(1,3);
    U f5 = f[5] * getCloverCoeff(2,3);    

    const int nodeSites = QDP::Layout::sitesOnNode();

    //tri.resize(nodeSites);  // hold local lattice

    JITCloverEnv::JITCloverMakeClovArg<U> arg = {diag_mass, f0,f1,f2,f3,f4,f5,tri_id };
    makeClovJIT( nodeSites , &arg );

    END_CODE();
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  template<typename T, typename U>
  void JitCloverTermT<T,U>::choles(int cb)
  {
    START_CODE();

    // When you are doing the cholesky - also fill out the trace_log_diag piece)
    // chlclovms(tr_log_diag_, cb);
    // Switch to LDL^\dag inversion
    ldagdlinv(tr_log_diag_,cb);

    END_CODE();
  }


  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   *
   * \return logarithm of the determinant  
   */
  template<typename T, typename U>
  Double JitCloverTermT<T,U>::cholesDet(int cb) const
  {
    START_CODE();

    if( choles_done[cb] == false ) 
    {
      QDPIO::cout << __func__ << ": Error: you have not done the Cholesky.on this operator on this subset" << endl;
      QDPIO::cout << "You sure you should not be asking invclov?" << endl;
      QDP_abort(1);
    }

    LatticeREAL ff=tr_log_diag_;

    if( param.sub_zero_usedP ) { 
 	QDPIO::cout << "Subtracting "<< param.sub_zero<<endl;
	LatticeREAL tmp;
	tmp[rb[cb]] = param.sub_zero;
	ff[rb[cb]] -= tmp;
    }
    END_CODE();

    // Need to thread generic sums in QDP++?
    // Need to thread generic norm2() in QDP++?
    return sum(ff, rb[cb]);
  }

  namespace JITCloverEnv { 
    template<typename U>
    struct LDagDLInvArgs { 
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;
      typedef OLattice< PScalar< PScalar< RScalar<REALT> > > > LatticeRealT;
      LatticeRealT& tr_log_diag;
      int tri_id;
      int cb;
    };


    template<typename U>
    inline 
    void LDagDLInvJIT(LDagDLInvArgs<U>* a) 
    {
      typedef typename LDagDLInvArgs<U>::REALT REALT;
      typedef typename LDagDLInvArgs<U>::RealT RealT;
      typedef typename LDagDLInvArgs<U>::LatticeRealT LatticeRealT;

      LatticeRealT& tr_log_diag = a->tr_log_diag;
      int tri_id = a->tri_id;
      int cb = a->cb;

      int num_sites = rb[cb].siteTable().size();

      while (1) {

	static QDPJit::SharedLibEntry sharedLibEntry;
	static MapVolumes*              mapVolumes;
	static string                   strId;
	static string                   prg;

	const int nodeSites = QDP::Layout::sitesOnNode();

	QDPJitArgs cudaArgs;

	int argNum = cudaArgs.addInt( rb[cb].numSiteTable() );
	int argOrd = cudaArgs.addBool( rb[cb].hasOrderedRep() );
	int argStart = cudaArgs.addInt( rb[cb].start() );
	int argSubset = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( rb[cb].getId() ) );
	int argDestPtr = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( tri_id ) );

	if (!mapVolumes) {
	  string strREALT;
	  string codeTr;

	  getTypeString( strREALT , REALT(0) );
	  if (!getCodeString( codeTr , tr_log_diag , "site", cudaArgs )) { QDP_info("LDagDLInvJIT: could not cache tr_log_diag"); break;  }      

	  ostringstream osId;
	  osId << "LDagDLInv "  << strREALT;
	  strId = osId.str();
	  xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	  cout << "strId = " << strId << endl;
#endif

	  std::ostringstream sprg;

	  sprg << " typedef " << strREALT << " REALT;" << endl;

	  sprg << " const int Nc = " << Nc << ";" << endl;
	  sprg << " RScalar<REALT> zip=0;" << endl;
	  sprg << " int N = 2*Nc;" << endl;

	  sprg << "	struct PrimitiveClovTriang" << endl;
	  sprg << "	{" << endl;
	  sprg << "   RScalar<REALT>   diag[2][2*Nc];" << endl;
	  sprg << "   RComplex<REALT>  offd[2][2*Nc*Nc-Nc];" << endl;
	  sprg << "	};" << endl;

	  sprg << "#define tri ((PrimitiveClovTriang*)( " << cudaArgs.getCode( argDestPtr ) << " ))" << endl;

	  sprg << "int site;" << endl;

	  sprg << "  if (" << cudaArgs.getCode(argOrd) << ") {" << endl;
	  sprg << "    site = " << cudaArgs.getCode(argStart) << ";" << endl;
	  sprg << "    site += blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	  sprg << "    if (site >= " << cudaArgs.getCode(argNum) << "+" <<  cudaArgs.getCode(argStart) << ") return;" << endl;
	  sprg << "  } else {" << endl;
	  sprg << "    int idx0 = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x; \n";
	  sprg << "    if (idx0 >= " << cudaArgs.getCode(argNum) << ") return;" << endl;
	  sprg << "    site = ((int*)" << cudaArgs.getCode(argSubset) << ")[idx0];" << endl;
	  sprg << "  }" << endl;

	  sprg << "	int site_neg_logdet=0;" << endl;

	  //
	  // Use a rolled version of the block loop ?
	  //
#if 0
	  sprg << "	for(int block=0; block < 2; block++) { " << endl;
	  sprg << "	  RScalar<REALT> inv_d[6];" << endl;
	  sprg << "	  RComplex<REALT> inv_offd[15];" << endl;
	  sprg << "	  RComplex<REALT> v[6];" << endl;
	  sprg << "	  RScalar<REALT>  diag_g[6];" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    inv_d[i] = tri[site].diag[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    inv_offd[i]  =tri[site].offd[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int j=0; j < N; ++j) { " << endl;
	  sprg << "	    for(int i=0; i < j; i++) { " << endl;
	  sprg << "	      int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "	      RComplex<REALT> A_ii = cmplx( inv_d[i], zip );" << endl;
	  sprg << "	      v[i] = A_ii*adj(inv_offd[elem_ji]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[j] = cmplx(inv_d[j],zip);" << endl;
	  sprg << "	    for(int k=0; k < j; k++) { " << endl;
	  sprg << "	      int elem_jk = j*(j-1)/2 + k;" << endl;
	  sprg << "	      v[j] -= inv_offd[elem_jk]*v[k];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[j] = real( v[j] );" << endl;
	  sprg << "	    for(int k=j+1; k < N; k++) { " << endl;
	  sprg << "	      int elem_kj = k*(k-1)/2 + j;" << endl;
	  sprg << "	      for(int l=0; l < j; l++) { " << endl;
	  sprg << "		int elem_kl = k*(k-1)/2 + l;" << endl;
	  sprg << "		inv_offd[elem_kj] -= inv_offd[elem_kl] * v[l];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      inv_offd[elem_kj] /= v[j];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  RScalar<REALT> one;" << endl;
	  sprg << "	  one.elem() = (REALT)1;" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    diag_g[i] = one/inv_d[i];" << endl;
	  sprg << "	    " << codeTr << ".elem().elem().elem() += log(fabs(inv_d[i].elem()));" << endl;
	  sprg << "	    if( inv_d[i].elem() < 0 ) { " << endl;
	  sprg << "	      site_neg_logdet++;" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  RComplex<REALT> sum;" << endl;
	  sprg << "	  for(int k = 0; k < N; ++k) {" << endl;
	  sprg << "	    for(int i = 0; i < k; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[k] = cmplx(diag_g[k],zip);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	      for(int j = k; j < i; ++j) {" << endl;
	  sprg << "		int elem_ij = i*(i-1)/2+j;	" << endl;
	  sprg << "		v[i] -= inv_offd[elem_ij] *inv_d[j]*v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      v[i] *= diag_g[i];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    for(int i = N-2; (int)i >= (int)k; --i) {" << endl;
	  sprg << "	      for(int j = i+1; j < N; ++j) {" << endl;
	  sprg << "		int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "		v[i] -= adj(inv_offd[elem_ji]) * v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[k] = real(v[k]);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      int elem_ik = i*(i-1)/2+k;" << endl;
	  sprg << "	      inv_offd[elem_ik] = v[i];" << endl;        // error
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    tri[site].diag[block][i] = inv_d[i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    tri[site].offd[block][i] = inv_offd[i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	}" << endl;
#else
	  sprg << "	int block=0;" << endl;
	  sprg << "	  RScalar<REALT> inv_d[6];" << endl;
	  sprg << "	  RComplex<REALT> inv_offd[15];" << endl;
	  sprg << "	  RComplex<REALT> v[6];" << endl;
	  sprg << "	  RScalar<REALT>  diag_g[6];" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    inv_d[i] = tri[site].diag[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    inv_offd[i]  =tri[site].offd[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int j=0; j < N; ++j) { " << endl;
	  sprg << "	    for(int i=0; i < j; i++) { " << endl;
	  sprg << "	      int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "	      RComplex<REALT> A_ii = cmplx( inv_d[i], zip );" << endl;
	  sprg << "	      v[i] = A_ii*adj(inv_offd[elem_ji]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[j] = cmplx(inv_d[j],zip);" << endl;
	  sprg << "	    for(int k=0; k < j; k++) { " << endl;
	  sprg << "	      int elem_jk = j*(j-1)/2 + k;" << endl;
	  sprg << "	      v[j] -= inv_offd[elem_jk]*v[k];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[j] = real( v[j] );" << endl;
	  sprg << "	    for(int k=j+1; k < N; k++) { " << endl;
	  sprg << "	      int elem_kj = k*(k-1)/2 + j;" << endl;
	  sprg << "	      for(int l=0; l < j; l++) { " << endl;
	  sprg << "		int elem_kl = k*(k-1)/2 + l;" << endl;
	  sprg << "		inv_offd[elem_kj] -= inv_offd[elem_kl] * v[l];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      inv_offd[elem_kj] /= v[j];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  RScalar<REALT> one;" << endl;
	  sprg << "	  one.elem() = (REALT)1;" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    diag_g[i] = one/inv_d[i];" << endl;
	  sprg << "	    " << codeTr << ".elem().elem().elem() += log(fabs(inv_d[i].elem()));" << endl;
	  sprg << "	    if( inv_d[i].elem() < 0 ) { " << endl;
	  sprg << "	      site_neg_logdet++;" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  RComplex<REALT> sum;" << endl;
	  sprg << "	  for(int k = 0; k < N; ++k) {" << endl;
	  sprg << "	    for(int i = 0; i < k; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[k] = cmplx(diag_g[k],zip);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	      for(int j = k; j < i; ++j) {" << endl;
	  sprg << "		int elem_ij = i*(i-1)/2+j;	" << endl;
	  sprg << "		v[i] -= inv_offd[elem_ij] *inv_d[j]*v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      v[i] *= diag_g[i];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    for(int i = N-2; (int)i >= (int)k; --i) {" << endl;
	  sprg << "	      for(int j = i+1; j < N; ++j) {" << endl;
	  sprg << "		int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "		v[i] -= adj(inv_offd[elem_ji]) * v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[k] = real(v[k]);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      int elem_ik = i*(i-1)/2+k;" << endl;
	  sprg << "	      inv_offd[elem_ik] = v[i];" << endl;        // error
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    tri[site].diag[block][i] = inv_d[i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    tri[site].offd[block][i] = inv_offd[i];" << endl;
	  sprg << "	  }" << endl;

	  sprg << "	block=1; " << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    inv_d[i] = tri[site].diag[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    inv_offd[i]  =tri[site].offd[block][i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int j=0; j < N; ++j) { " << endl;
	  sprg << "	    for(int i=0; i < j; i++) { " << endl;
	  sprg << "	      int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "	      RComplex<REALT> A_ii = cmplx( inv_d[i], zip );" << endl;
	  sprg << "	      v[i] = A_ii*adj(inv_offd[elem_ji]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[j] = cmplx(inv_d[j],zip);" << endl;
	  sprg << "	    for(int k=0; k < j; k++) { " << endl;
	  sprg << "	      int elem_jk = j*(j-1)/2 + k;" << endl;
	  sprg << "	      v[j] -= inv_offd[elem_jk]*v[k];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[j] = real( v[j] );" << endl;
	  sprg << "	    for(int k=j+1; k < N; k++) { " << endl;
	  sprg << "	      int elem_kj = k*(k-1)/2 + j;" << endl;
	  sprg << "	      for(int l=0; l < j; l++) { " << endl;
	  sprg << "		int elem_kl = k*(k-1)/2 + l;" << endl;
	  sprg << "		inv_offd[elem_kj] -= inv_offd[elem_kl] * v[l];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      inv_offd[elem_kj] /= v[j];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  one.elem() = (REALT)1;" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    diag_g[i] = one/inv_d[i];" << endl;
	  sprg << "	    " << codeTr << ".elem().elem().elem() += log(fabs(inv_d[i].elem()));" << endl;
	  sprg << "	    if( inv_d[i].elem() < 0 ) { " << endl;
	  sprg << "	      site_neg_logdet++;" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int k = 0; k < N; ++k) {" << endl;
	  sprg << "	    for(int i = 0; i < k; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    v[k] = cmplx(diag_g[k],zip);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      zero_rep(v[i]);" << endl;
	  sprg << "	      for(int j = k; j < i; ++j) {" << endl;
	  sprg << "		int elem_ij = i*(i-1)/2+j;	" << endl;
	  sprg << "		v[i] -= inv_offd[elem_ij] *inv_d[j]*v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	      v[i] *= diag_g[i];" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    for(int i = N-2; (int)i >= (int)k; --i) {" << endl;
	  sprg << "	      for(int j = i+1; j < N; ++j) {" << endl;
	  sprg << "		int elem_ji = j*(j-1)/2 + i;" << endl;
	  sprg << "		v[i] -= adj(inv_offd[elem_ji]) * v[j];" << endl;
	  sprg << "	      }" << endl;
	  sprg << "	    }" << endl;
	  sprg << "	    inv_d[k] = real(v[k]);" << endl;
	  sprg << "	    for(int i = k+1; i < N; ++i) {" << endl;
	  sprg << "	      int elem_ik = i*(i-1)/2+k;" << endl;
	  sprg << "	      inv_offd[elem_ik] = v[i];" << endl;        // error
	  sprg << "	    }" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < N; i++) { " << endl;
	  sprg << "	    tri[site].diag[block][i] = inv_d[i];" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	  for(int i=0; i < 15; i++) { " << endl;
	  sprg << "	    tri[site].offd[block][i] = inv_offd[i];" << endl;
	  sprg << "	  }" << endl;
#endif

	  sprg << "	if( site_neg_logdet != 0 ) { " << endl;
	  sprg << "	}" << endl;

	  prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	  cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	} else {
	  if (!cacheLock( tr_log_diag , cudaArgs )) { QDP_info("eval: could not cache tr_log_diag");  break;   }
	}

	if (!QDPJit::Instance()( strId , prg , cudaArgs.getDevPtr() , rb[cb].numSiteTable() , sharedLibEntry  , mapVolumes )) {
	  QDP_info("packClov call to cuda jitter failed");
	  break;
	}

	return;
      }
      QDP_error("packClov host not implemented");
    }



  } /* End Namespace */


  /*! An LDL^\dag decomposition and inversion? */
  template<typename T, typename U>
  void JitCloverTermT<T,U>::ldagdlinv(LatticeREAL& tr_log_diag, int cb)
  {
    START_CODE();

    if ( 2*Nc < 3 )
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << endl;
      QDP_abort(1);
    }

    // Zero trace log
    tr_log_diag = zero;

    JITCloverEnv::LDagDLInvArgs<U> a = { tr_log_diag, tri_id, cb };
    LDagDLInvJIT(&a);

    
    // This comes from the days when we used to do Cholesky
    choles_done[cb] = true;
    END_CODE();
  }
 
  /*! CHLCLOVMS - Cholesky decompose the clover mass term and uses it to
   *              compute  lower(A^-1) = lower((L.L^dag)^-1)
   *              Adapted from Golub and Van Loan, Matrix Computations, 2nd, Sec 4.2.4
   *
   * Arguments:
   *
   * \param DetP         flag whether to compute determinant (Read)
   * \param logdet       logarithm of the determinant        (Write)
   * \param cb           checkerboard of work                (Read)
   */

  namespace JITCloverEnv { 

    template<typename U>
    inline 
    void cholesSiteLoop(int lo, int hi, int myId, LDagDLInvArgs<U>* a)
    {
      QDP_error_exit("Clover: cholesSiteLoop n.i.");
#if 0
      typedef typename LDagDLInvArgs<U>::REALT REALT;
      typedef typename LDagDLInvArgs<U>::RealT RealT;
      typedef typename LDagDLInvArgs<U>::LatticeRealT LatticeRealT;

      LatticeRealT& tr_log_diag = a->tr_log_diag;
      multi1d<PrimitiveClovTriang <REALT> >& tri = a->tri;
      int cb = a->cb;
      
      int n = 2*Nc;
      /*# Cholesky decompose  A = L.L^dag */
      /*# NOTE!!: I can store this matrix in  invclov, but will need a */
      /*#   temporary  diag */
      for(int ssite=lo; ssite < hi; ++ssite)  {
	int site = rb[cb].siteTable()[ssite];

	PrimitiveClovTriang<REALT>  invcl;

	multi1d< RScalar<REALT> > diag_g(n);
	multi1d< RComplex<REALT> > v1(n);
	RComplex<REALT> sum;
	RScalar<REALT> one;
	RScalar<REALT> zero;
	RScalar<REALT> lrtmp;
	
	one = 1;
	zero = 0;
  
	for(int s = 0; s < 2; ++s) {
	  
	  int elem_jk = 0;
	  int elem_ij;
	  
	  for(int j = 0; j <  n; ++j) {

	    /*# Multiply clover mass term against basis vector.  */
	    /*# Actually, I need a column of the lower triang matrix clov. */
	    v1[j] = cmplx(tri[site].diag[s][j],zero);
    
	    elem_ij = elem_jk + 2*j;
	    
	    for(int i = j+1; i < n; ++i) {
	      v1[i] = tri[site].offd[s][elem_ij];
	      elem_ij += i;
	    }
      
	    /*# Back to cholesky */
	    /*# forward substitute */
	    for(int k = 0; k < j; ++k) {

	      int elem_ik = elem_jk;
	
	      for(int i = j; i < n; ++i) {
		
		v1[i] -= adj(invcl.offd[s][elem_jk]) * invcl.offd[s][elem_ik];
		elem_ik += i;
	      }
	      elem_jk++;
	    }

	    /*# The diagonal is (should be!!) real and positive */
	    diag_g[j] = real(v1[j]);
	  
	    /*#+ */
	    /*# Squeeze in computation of the trace log of the diagonal term */
	    /*#- */
	    if ( diag_g[j].elem() > 0 ) {
	      
	      lrtmp = log(diag_g[j]);
	    }
	    else {

	      // Make sure any node can print this message
	      cerr << "Clover term has negative diagonal element: "
		   << "diag_g[" << j << "]= " << diag_g[j] 
		   << " at site: " << site << endl;
	      QDP_abort(1);
	    }

	    tr_log_diag.elem(site).elem().elem() += lrtmp;
	  
	    diag_g[j] = sqrt(diag_g[j]);
	    diag_g[j] = one / diag_g[j];
      
	    /*# backward substitute */
	    elem_ij = elem_jk + j;
	    for(int i = j+1; i < n; ++i) {

	      invcl.offd[s][elem_ij] = v1[i] * diag_g[j];
	      elem_ij += i;
	    }
	  }

	  /*# Use forward and back substitution to construct  invcl.offd = lower(A^-1) */
	  for(int k = 0; k < n; ++k) {

	    for(int i = 0; i < k; ++i)
	      zero_rep(v1[i]);
	  
	    /*# Forward substitution */
	    v1[k] = cmplx(diag_g[k],zero);
      
	    for(int i = k+1; i < n; ++i) {
	  
	      zero_rep(sum);
	      elem_ij = i*(i-1)/2+k;	
	      for(int j = k; j < i; ++j) {

		sum -= invcl.offd[s][elem_ij] * v1[j];
		elem_ij++;
	      }
	
	      v1[i] = sum * diag_g[i];
	    }
      
	    /*# Backward substitution */
	    v1[n-1] = v1[n-1] * diag_g[n-1];
     
	    for(int i = n-2; (int)i >= (int)k; --i) {

	      sum = v1[i];
	    
	      int elem_ji = ((i+1)*i)/2+i;
	      for(int j = i+1; j < n; ++j) {

		sum -= adj(invcl.offd[s][elem_ji]) * v1[j];
		elem_ji += j;
	      }
	      v1[i] = sum * diag_g[i];
	    }

	    /*# Overwrite column k of invcl.offd */
	    invcl.diag[s][k] = real(v1[k]);

	    int elem_ik = ((k+1)*k)/2+k;
      
	    for(int i = k+1; i < n; ++i) {

	      invcl.offd[s][elem_ik] = v1[i];
	      elem_ik += i;
	    }
	  }
	}

	// Overwrite original element
	tri[site] = invcl;
      }
#endif
    } // End function

  } // End Namespace


  template<typename T, typename U>
  void JitCloverTermT<T,U>::chlclovms(LatticeREAL& tr_log_diag, int cb)
  {
    QDP_error_exit("Clover: chlclovms n.i.");
#if 0
    START_CODE();

    if ( 2*Nc < 3 )
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << endl;
      QDP_abort(1);
    }
  
    tr_log_diag = zero;
    JITCloverEnv::LDagDLInvArgs<U> a = { tr_log_diag, tri, cb};
    dispatch_to_threads(rb[cb].numSiteTable(), a, JITCloverEnv::cholesSiteLoop<U>);
    
    choles_done[cb] = true;
    END_CODE();
#endif
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
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  template<typename T, typename U>
  void JitCloverTermT<T,U>::applySite(T& chi, const T& psi, 
			    enum PlusMinus isign, int site) const
  {
    QDP_error_exit("Clover: applySite n.i.");
#if 0
    START_CODE();

    if ( Ns != 4 )
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << endl;
      QDP_abort(1);
    }

    int n = 2*Nc;

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));


    cchi[ 0] = tri[site].diag[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 0]) * ppsi[ 1]
      +   conj(tri[site].offd[0][ 1]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 3]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 6]) * ppsi[ 4]
      +   conj(tri[site].offd[0][10]) * ppsi[ 5];
    
    cchi[ 1] = tri[site].diag[0][ 1]  * ppsi[ 1]
      +        tri[site].offd[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 2]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 4]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 7]) * ppsi[ 4]
      +   conj(tri[site].offd[0][11]) * ppsi[ 5];
    
    cchi[ 2] = tri[site].diag[0][ 2]  * ppsi[ 2]
      +        tri[site].offd[0][ 1]  * ppsi[ 0]
      +        tri[site].offd[0][ 2]  * ppsi[ 1]
      +   conj(tri[site].offd[0][ 5]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 8]) * ppsi[ 4]
      +   conj(tri[site].offd[0][12]) * ppsi[ 5];
    
    cchi[ 3] = tri[site].diag[0][ 3]  * ppsi[ 3]
      +        tri[site].offd[0][ 3]  * ppsi[ 0]
      +        tri[site].offd[0][ 4]  * ppsi[ 1]
      +        tri[site].offd[0][ 5]  * ppsi[ 2]
      +   conj(tri[site].offd[0][ 9]) * ppsi[ 4]
      +   conj(tri[site].offd[0][13]) * ppsi[ 5];
    
    cchi[ 4] = tri[site].diag[0][ 4]  * ppsi[ 4]
      +        tri[site].offd[0][ 6]  * ppsi[ 0]
      +        tri[site].offd[0][ 7]  * ppsi[ 1]
      +        tri[site].offd[0][ 8]  * ppsi[ 2]
      +        tri[site].offd[0][ 9]  * ppsi[ 3]
      +   conj(tri[site].offd[0][14]) * ppsi[ 5];
    
    cchi[ 5] = tri[site].diag[0][ 5]  * ppsi[ 5]
      +        tri[site].offd[0][10]  * ppsi[ 0]
      +        tri[site].offd[0][11]  * ppsi[ 1]
      +        tri[site].offd[0][12]  * ppsi[ 2]
      +        tri[site].offd[0][13]  * ppsi[ 3]
      +        tri[site].offd[0][14]  * ppsi[ 4];
    
    cchi[ 6] = tri[site].diag[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 0]) * ppsi[ 7]
      +   conj(tri[site].offd[1][ 1]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 3]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 6]) * ppsi[10]
      +   conj(tri[site].offd[1][10]) * ppsi[11];
    
    cchi[ 7] = tri[site].diag[1][ 1]  * ppsi[ 7]
      +        tri[site].offd[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 2]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 4]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 7]) * ppsi[10]
      +   conj(tri[site].offd[1][11]) * ppsi[11];
    
    cchi[ 8] = tri[site].diag[1][ 2]  * ppsi[ 8]
      +        tri[site].offd[1][ 1]  * ppsi[ 6]
      +        tri[site].offd[1][ 2]  * ppsi[ 7]
      +   conj(tri[site].offd[1][ 5]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 8]) * ppsi[10]
      +   conj(tri[site].offd[1][12]) * ppsi[11];
    
    cchi[ 9] = tri[site].diag[1][ 3]  * ppsi[ 9]
      +        tri[site].offd[1][ 3]  * ppsi[ 6]
      +        tri[site].offd[1][ 4]  * ppsi[ 7]
      +        tri[site].offd[1][ 5]  * ppsi[ 8]
      +   conj(tri[site].offd[1][ 9]) * ppsi[10]
      +   conj(tri[site].offd[1][13]) * ppsi[11];
    
    cchi[10] = tri[site].diag[1][ 4]  * ppsi[10]
      +        tri[site].offd[1][ 6]  * ppsi[ 6]
      +        tri[site].offd[1][ 7]  * ppsi[ 7]
      +        tri[site].offd[1][ 8]  * ppsi[ 8]
      +        tri[site].offd[1][ 9]  * ppsi[ 9]
      +   conj(tri[site].offd[1][14]) * ppsi[11];
    
    cchi[11] = tri[site].diag[1][ 5]  * ppsi[11]
      +        tri[site].offd[1][10]  * ppsi[ 6]
      +        tri[site].offd[1][11]  * ppsi[ 7]
      +        tri[site].offd[1][12]  * ppsi[ 8]
      +        tri[site].offd[1][13]  * ppsi[ 9]
      +        tri[site].offd[1][14]  * ppsi[10];


    END_CODE();
#endif
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

  namespace JITCloverEnv { 
    template<typename U>
    struct TriaCntrArgs {
      typedef typename WordType<U>::Type_t REALT;
      
      U& B;
      int tri_id;
      int mat;
      int cb;
    };
  


    template<typename U>
    inline 
    void triaCntrJIT(TriaCntrArgs<U>* a)
    {
      typedef typename WordType<U>::Type_t REALT;
      U& B = a->B;
      int tri_id = a->tri_id;
      int mat = a->mat;
      int cb  = a->cb;

      switch( mat ){
      case 0:
      case 3:
      case 5:
      case 6:
      case 9:
      case 10:
      case 12:
	break;
      default:
	QDPIO::cout << __func__ << ": invalid Gamma matrix int" << endl;
	QDP_abort(1);
      }
      
      while (1) {

	static QDPJit::SharedLibEntry sharedLibEntry;
	static MapVolumes*              mapVolumes;
	static string                   strId;
	static string                   prg;

	QDPJitArgs cudaArgs;

	const int nodeSites = QDP::Layout::sitesOnNode();

	int argNum = cudaArgs.addInt( rb[cb].numSiteTable() );
	int argOrd = cudaArgs.addBool( rb[cb].hasOrderedRep() );
	int argStart = cudaArgs.addInt( rb[cb].start() );
	int argSubset = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( rb[cb].getId() ) );
	int argDestPtr = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( tri_id ) );
	int argMat = cudaArgs.addInt( mat );

	if (!mapVolumes) {
	  string strREALT,typeT;
	  string codeB;

	  getTypeString( strREALT , REALT(0) );
	  getTypeString( typeT , B, cudaArgs );

	  if (!getCodeString( codeB , B , "site", cudaArgs )) { QDP_info("Clover: triaCntrJIT: could not cache B"); break;  }      

	  ostringstream osId;
	  osId << "triaCntrJIT "  << typeT;
	  strId = osId.str();
	  xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	  cout << "strId = " << strId << endl;
#endif

	  std::ostringstream sprg;

	  sprg << " typedef " << strREALT << " REALT;" << endl;
	  sprg << " const int Nc = " << Nc << ";" << endl;

	  sprg << "	struct PrimitiveClovTriang" << endl;
	  sprg << "	{" << endl;
	  sprg << "   RScalar<REALT>   diag[2][2*Nc];" << endl;
	  sprg << "   RComplex<REALT>  offd[2][2*Nc*Nc-Nc];" << endl;
	  sprg << "	};" << endl;

	  sprg << "#define tri ((PrimitiveClovTriang*)( " << cudaArgs.getCode( argDestPtr ) << " ))" << endl;

	  sprg << "int site;" << endl;

	  sprg << "  if (" << cudaArgs.getCode(argOrd) << ") {" << endl;
	  sprg << "    site = " << cudaArgs.getCode(argStart) << ";" << endl;
	  sprg << "    site += blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	  sprg << "    if (site >= " << cudaArgs.getCode(argNum) << "+" <<  cudaArgs.getCode(argStart) << ") return;" << endl;
	  sprg << "  } else {" << endl;
	  sprg << "    int idx0 = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x; \n";
	  sprg << "    if (idx0 >= " << cudaArgs.getCode(argNum) << ") return;" << endl;
	  sprg << "    site = ((int*)" << cudaArgs.getCode(argSubset) << ")[idx0];" << endl;
	  sprg << "  }" << endl;

	  sprg << " int mat = " << cudaArgs.getCode( argMat ) << ";" << endl;

	  sprg << "     switch( mat ){ " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 0: " << endl;
	  sprg << "       /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */ " << endl;
	  sprg << "       /*#               0  1  0  0 */ " << endl;
	  sprg << "       /*#               0  0  1  0 */ " << endl;
	  sprg << "       /*#               0  0  0  1 */ " << endl;
	  sprg << "       /*# From diagonal part */ " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	RComplex<REALT> lctmp0; " << endl;
	  sprg << " 	RScalar<REALT> lr_zero0; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp0; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	lr_zero0 = 0; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i0 = 0; i0 < Nc; ++i0) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  lrtmp0 = tri[site].diag[0][i0]; " << endl;
	  sprg << " 	  lrtmp0 += tri[site].diag[0][i0+Nc]; " << endl;
	  sprg << " 	  lrtmp0 += tri[site].diag[1][i0]; " << endl;
	  sprg << " 	  lrtmp0 += tri[site].diag[1][i0+Nc]; " << endl;
	  sprg << " 	  " << codeB << ".elem().elem(i0,i0) = cmplx(lrtmp0,lr_zero0); " << endl;
	  sprg << " 	} " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	/*# From lower triangular portion */ " << endl;
	  sprg << " 	int elem_ij0 = 0; " << endl;
	  sprg << " 	for(int i0 = 1; i0 < Nc; ++i0) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j0 = 0; j0 < i0; ++j0) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp0 = tri[site].offd[0][elem_ij0]; " << endl;
	  sprg << " 	    lctmp0 += tri[site].offd[0][elem_ijb0]; " << endl;
	  sprg << " 	    lctmp0 += tri[site].offd[1][elem_ij0]; " << endl;
	  sprg << " 	    lctmp0 += tri[site].offd[1][elem_ijb0]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(j0,i0) = lctmp0; " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i0,j0) = adj(lctmp0); " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij0++; " << endl;
	  sprg << " 	    elem_ijb0++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << " 	   " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 3: " << endl;
	  sprg << "       /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */ " << endl;
	  sprg << "       /*#               0  i  0  0 */ " << endl;
	  sprg << "       /*#               0  0 -i  0 */ " << endl;
	  sprg << "       /*#               0  0  0  i */ " << endl;
	  sprg << "       /*# From diagonal part */ " << endl;
	  sprg << " 	   " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	RComplex<REALT> lctmp3; " << endl;
	  sprg << " 	RScalar<REALT> lr_zero3; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp3; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	lr_zero3 = 0; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i3 = 0; i3 < Nc; ++i3) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  lrtmp3 = tri[site].diag[0][i3+Nc]; " << endl;
	  sprg << " 	  lrtmp3 -= tri[site].diag[0][i3]; " << endl;
	  sprg << " 	  lrtmp3 -= tri[site].diag[1][i3]; " << endl;
	  sprg << " 	  lrtmp3 += tri[site].diag[1][i3+Nc]; " << endl;
	  sprg << " 	  " << codeB << ".elem().elem(i3,i3) = cmplx(lr_zero3,lrtmp3); " << endl;
	  sprg << " 	} " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	/*# From lower triangular portion */ " << endl;
	  sprg << " 	int elem_ij3 = 0; " << endl;
	  sprg << " 	for(int i3 = 1; i3 < Nc; ++i3) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j3 = 0; j3 < i3; ++j3) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp3 = tri[site].offd[0][elem_ijb3]; " << endl;
	  sprg << " 	    lctmp3 -= tri[site].offd[0][elem_ij3]; " << endl;
	  sprg << " 	    lctmp3 -= tri[site].offd[1][elem_ij3]; " << endl;
	  sprg << " 	    lctmp3 += tri[site].offd[1][elem_ijb3]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(j3,i3) = timesI(adj(lctmp3)); " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i3,j3) = timesI(lctmp3); " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij3++; " << endl;
	  sprg << " 	    elem_ijb3++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 5: " << endl;
	  sprg << "       /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */ " << endl;
	  sprg << "       /*#               1  0  0  0 */ " << endl;
	  sprg << "       /*#               0  0  0 -1 */ " << endl;
	  sprg << "       /*#               0  0  1  0 */ " << endl;
	  sprg << " 	   " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	RComplex<REALT> lctmp5; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp5; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i5 = 0; i5 < Nc; ++i5) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ij5 = (i5+Nc)*(i5+Nc-1)/2; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j5 = 0; j5 < Nc; ++j5) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    int elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp5 = adj(tri[site].offd[0][elem_ji5]); " << endl;
	  sprg << " 	    lctmp5 -= tri[site].offd[0][elem_ij5]; " << endl;
	  sprg << " 	    lctmp5 += adj(tri[site].offd[1][elem_ji5]); " << endl;
	  sprg << " 	    lctmp5 -= tri[site].offd[1][elem_ij5]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i5,j5) = lctmp5; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij5++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 6: " << endl;
	  sprg << "       /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */ " << endl;
	  sprg << "       /*#              -i  0  0  0 */ " << endl;
	  sprg << "       /*#               0  0  0 -i */ " << endl;
	  sprg << "       /*#               0  0 -i  0 */ " << endl;
	  sprg << " 	   " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	RComplex<REALT> lctmp6; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp6; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i6 = 0; i6 < Nc; ++i6) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ij6 = (i6+Nc)*(i6+Nc-1)/2; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j6 = 0; j6 < Nc; ++j6) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    int elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp6 = adj(tri[site].offd[0][elem_ji6]); " << endl;
	  sprg << " 	    lctmp6 += tri[site].offd[0][elem_ij6]; " << endl;
	  sprg << " 	    lctmp6 += adj(tri[site].offd[1][elem_ji6]); " << endl;
	  sprg << " 	    lctmp6 += tri[site].offd[1][elem_ij6]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i6,j6) = timesMinusI(lctmp6); " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij6++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 9: " << endl;
	  sprg << "       /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */ " << endl;
	  sprg << "       /*#               i  0  0  0 */ " << endl;
	  sprg << "       /*#               0  0  0 -i */ " << endl;
	  sprg << "       /*#               0  0 -i  0 */ " << endl;
	  sprg << " 	   " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	RComplex<REALT> lctmp9; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp9; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i9 = 0; i9 < Nc; ++i9) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ij9 = (i9+Nc)*(i9+Nc-1)/2; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j9 = 0; j9 < Nc; ++j9) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    int elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp9 = adj(tri[site].offd[0][elem_ji9]); " << endl;
	  sprg << " 	    lctmp9 += tri[site].offd[0][elem_ij9]; " << endl;
	  sprg << " 	    lctmp9 -= adj(tri[site].offd[1][elem_ji9]); " << endl;
	  sprg << " 	    lctmp9 -= tri[site].offd[1][elem_ij9]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i9,j9) = timesI(lctmp9); " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij9++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	   " << endl;
	  sprg << "     case 10: " << endl;
	  sprg << "       /*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */ " << endl;
	  sprg << "       /*#               1  0  0  0 */ " << endl;
	  sprg << "       /*#               0  0  0  1 */ " << endl;
	  sprg << "       /*#               0  0 -1  0 */ " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	RComplex<REALT> lctmp10; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp10; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i10 = 0; i10 < Nc; ++i10) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ij10 = (i10+Nc)*(i10+Nc-1)/2; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j10 = 0; j10 < Nc; ++j10) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    int elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp10 = adj(tri[site].offd[0][elem_ji10]); " << endl;
	  sprg << " 	    lctmp10 -= tri[site].offd[0][elem_ij10]; " << endl;
	  sprg << " 	    lctmp10 -= adj(tri[site].offd[1][elem_ji10]); " << endl;
	  sprg << " 	    lctmp10 += tri[site].offd[1][elem_ij10]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i10,j10) = lctmp10; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij10++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     case 12: " << endl;
	  sprg << "       /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */ " << endl;
	  sprg << "       /*#               0 -i  0  0 */ " << endl;
	  sprg << "       /*#               0  0 -i  0 */ " << endl;
	  sprg << "       /*#               0  0  0  i */ " << endl;
	  sprg << "       /*# From diagonal part */ " << endl;
	  sprg << "       { " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	RComplex<REALT> lctmp12; " << endl;
	  sprg << " 	RScalar<REALT> lr_zero12; " << endl;
	  sprg << " 	RScalar<REALT> lrtmp12; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	lr_zero12 = 0; " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	for(int i12 = 0; i12 < Nc; ++i12) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  lrtmp12 = tri[site].diag[0][i12]; " << endl;
	  sprg << " 	  lrtmp12 -= tri[site].diag[0][i12+Nc]; " << endl;
	  sprg << " 	  lrtmp12 -= tri[site].diag[1][i12]; " << endl;
	  sprg << " 	  lrtmp12 += tri[site].diag[1][i12+Nc]; " << endl;
	  sprg << " 	  " << codeB << ".elem().elem(i12,i12) = cmplx(lr_zero12,lrtmp12); " << endl;
	  sprg << " 	} " << endl;
	  sprg << " 	     " << endl;
	  sprg << " 	/*# From lower triangular portion */ " << endl;
	  sprg << " 	int elem_ij12 = 0; " << endl;
	  sprg << " 	for(int i12 = 1; i12 < Nc; ++i12) { " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  int elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc; " << endl;
	  sprg << " 	       " << endl;
	  sprg << " 	  for(int j12 = 0; j12 < i12; ++j12) { " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    lctmp12 = tri[site].offd[0][elem_ij12]; " << endl;
	  sprg << " 	    lctmp12 -= tri[site].offd[0][elem_ijb12]; " << endl;
	  sprg << " 	    lctmp12 -= tri[site].offd[1][elem_ij12]; " << endl;
	  sprg << " 	    lctmp12 += tri[site].offd[1][elem_ijb12]; " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(i12,j12) = timesI(lctmp12); " << endl;
	  sprg << " 	    " << codeB << ".elem().elem(j12,i12) = timesI(adj(lctmp12)); " << endl;
	  sprg << " 		 " << endl;
	  sprg << " 	    elem_ij12++; " << endl;
	  sprg << " 	    elem_ijb12++; " << endl;
	  sprg << " 	  } " << endl;
	  sprg << " 	} " << endl;
	  sprg << "       } " << endl;
	  sprg << "       break; " << endl;
	  sprg << "     } " << endl;

	  prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	  cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	} else {
	  if (!cacheLock(  B , cudaArgs )) { QDP_info("eval: could not cache B");  break;   }
	}

	if (!QDPJit::Instance()( strId , prg , cudaArgs.getDevPtr() , rb[cb].numSiteTable() , sharedLibEntry  , mapVolumes )) {
	  QDP_info("packClov call to cuda jitter failed");
	  break;
	}

	return;
      }
      QDP_error("Clover: triaCntrJIT: host not implemented");
    }


  }



   
  template<typename T, typename U>
  void JitCloverTermT<T,U>::triacntr(U& B, int mat, int cb) const
  {
    START_CODE();

    B = zero;

    if ( mat < 0  ||  mat > 15 )
    {
      QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << endl;
      QDP_abort(1);
    }

    JITCloverEnv::TriaCntrArgs<U> a = { B, tri_id , mat, cb };
    triaCntrJIT( &a );

    END_CODE();
  }

  //! Returns the appropriate clover coefficient for indices mu and nu
  template<typename T, typename U>
  Real
  JitCloverTermT<T,U>::getCloverCoeff(int mu, int nu) const 
  { 
    START_CODE();

    if( param.anisoParam.anisoP )  {
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


  namespace JITCloverEnv { 

    template<typename T>
    struct ApplyArgs {
      typedef typename WordType<T>::Type_t REALT;
      T& chi;
      const T& psi;
      int tri_id;
      int cb;
    };



    template<typename T>
    void applyJIT( ApplyArgs<T>* arg)
    {
      typedef typename WordType<T>::Type_t REALT;

      T& chi=arg->chi;
      const T& psi=arg->psi;
      int tri_id = arg->tri_id;
      int cb = arg->cb;

      while (1) {

	static QDPJit::SharedLibEntry sharedLibEntry;
	static MapVolumes*              mapVolumes;
	static string                   strId;
	static string                   prg;

	const int nodeSites = QDP::Layout::sitesOnNode();

	QDPJitArgs cudaArgs;

	int argNum = cudaArgs.addInt( rb[cb].numSiteTable() );
	int argOrd = cudaArgs.addBool( rb[cb].hasOrderedRep() );
	int argStart = cudaArgs.addInt( rb[cb].start() );
	int argSubset = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( rb[cb].getId() ) );
	int argDestPtr = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( tri_id ) );

	if (!mapVolumes) {
	  string strREALT,typeT;
	  string codeChi,codePsi;

	  getTypeString( strREALT , REALT(0) );
	  getTypeString( typeT , chi, cudaArgs );

	  if (!getCodeString( codePsi , psi , "site", cudaArgs )) { QDP_info("Clover: applyJIT: could not cache psi");break;}
	  if (!getCodeString( codeChi , chi , "site", cudaArgs )) { QDP_info("Clover: applyJIT: could not cache chi");break;}

	  ostringstream osId;
	  osId << "applyJIT "  << typeT;
	  strId = osId.str();
	  xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	  cout << "strId = " << strId << endl;
#endif

	  std::ostringstream sprg;

	  sprg << " typedef " << strREALT << " REALT;" << endl;
	  sprg << " const int Nc = " << Nc << ";" << endl;

	  sprg << "	struct PrimitiveClovTriang" << endl;
	  sprg << "	{" << endl;
	  sprg << "   RScalar<REALT>   diag[2][2*Nc];" << endl;
	  sprg << "   RComplex<REALT>  offd[2][2*Nc*Nc-Nc];" << endl;
	  sprg << "	};" << endl;

	  sprg << "#define tri ((PrimitiveClovTriang*)( " << cudaArgs.getCode( argDestPtr ) << " ))" << endl;

	  sprg << "int site;" << endl;

	  sprg << "  if (" << cudaArgs.getCode(argOrd) << ") {" << endl;
	  sprg << "    site = " << cudaArgs.getCode(argStart) << ";" << endl;
	  sprg << "    site += blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	  sprg << "    if (site >= " << cudaArgs.getCode(argNum) << "+" <<  cudaArgs.getCode(argStart) << ") return;" << endl;
	  sprg << "  } else {" << endl;
	  sprg << "    int idx0 = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x; \n";
	  sprg << "    if (idx0 >= " << cudaArgs.getCode(argNum) << ") return;" << endl;
	  sprg << "    site = ((int*)" << cudaArgs.getCode(argSubset) << ")[idx0];" << endl;
	  sprg << "  }" << endl;

	  sprg << "RComplex<REALT>* cchi = (RComplex<REALT>*)&(" << codeChi << ".elem(0).elem(0));" << endl;
	  sprg << "const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(" << codePsi << ".elem(0).elem(0));" << endl;

	  sprg << "     // Unrolled version 3.  " << endl;
	  sprg << "  " << endl;
	  sprg << "     cchi[ 0].real()  = tri[site].diag[0][0].elem()  * ppsi[0].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][0].real()  * ppsi[1].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][0].imag()  * ppsi[1].imag(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][1].real()  * ppsi[2].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][1].imag()  * ppsi[2].imag(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][3].real()  * ppsi[3].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][3].imag()  * ppsi[3].imag(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][6].real()  * ppsi[4].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][6].imag()  * ppsi[4].imag(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][10].real() * ppsi[5].real(); " << endl;
	  sprg << "     cchi[ 0].real() += tri[site].offd[0][10].imag() * ppsi[5].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 0].imag()  = tri[site].diag[0][0].elem()  * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() += tri[site].offd[0][0].real()  * ppsi[1].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() -= tri[site].offd[0][0].imag()  * ppsi[1].real(); " << endl;
	  sprg << "     cchi[ 0].imag() += tri[site].offd[0][3].real()  * ppsi[3].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() -= tri[site].offd[0][3].imag()  * ppsi[3].real(); " << endl;
	  sprg << "     cchi[ 0].imag() += tri[site].offd[0][1].real()  * ppsi[2].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() -= tri[site].offd[0][1].imag()  * ppsi[2].real(); " << endl;
	  sprg << "     cchi[ 0].imag() += tri[site].offd[0][6].real()  * ppsi[4].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() -= tri[site].offd[0][6].imag()  * ppsi[4].real(); " << endl;
	  sprg << "     cchi[ 0].imag() += tri[site].offd[0][10].real() * ppsi[5].imag(); " << endl;
	  sprg << "     cchi[ 0].imag() -= tri[site].offd[0][10].imag() * ppsi[5].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 1].real()  = tri[site].diag[0][ 1].elem() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 0].real() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 1].real() -= tri[site].offd[0][ 0].imag() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 2].real() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 2].imag() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 4].real() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 4].imag() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 7].real() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][ 7].imag() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][11].real() * ppsi[ 5].real(); " << endl;
	  sprg << "     cchi[ 1].real() += tri[site].offd[0][11].imag() * ppsi[ 5].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 1].imag()  = tri[site].diag[0][ 1].elem() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][ 0].real() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][ 0].imag() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][ 2].real() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() -= tri[site].offd[0][ 2].imag() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][ 4].real() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() -= tri[site].offd[0][ 4].imag() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][ 7].real() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() -= tri[site].offd[0][ 7].imag() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 1].imag() += tri[site].offd[0][11].real() * ppsi[ 5].imag(); " << endl;
	  sprg << "     cchi[ 1].imag() -= tri[site].offd[0][11].imag() * ppsi[ 5].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 2].real() = tri[site].diag[0][ 2].elem()  * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][ 1].real() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 2].real() -= tri[site].offd[0][ 1].imag() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][ 2].real() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 2].real() -= tri[site].offd[0][ 2].imag() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][5].real()  * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][5].imag()  * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][8].real()  * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][8].imag()  * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][12].real() * ppsi[ 5].real(); " << endl;
	  sprg << "     cchi[ 2].real() += tri[site].offd[0][12].imag() * ppsi[ 5].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 2].imag() = tri[site].diag[0][ 2].elem()  * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][ 1].real() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][ 1].imag() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][ 2].real() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][ 2].imag() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][5].real()  * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() -= tri[site].offd[0][5].imag()  * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][8].real()  * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() -= tri[site].offd[0][8].imag()  * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 2].imag() += tri[site].offd[0][12].real() * ppsi[ 5].imag(); " << endl;
	  sprg << "     cchi[ 2].imag() -= tri[site].offd[0][12].imag() * ppsi[ 5].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 3].real()  = tri[site].diag[0][ 3].elem() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][ 3].real() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 3].real() -= tri[site].offd[0][ 3].imag() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][ 4].real() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 3].real() -= tri[site].offd[0][ 4].imag() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][ 5].real() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 3].real() -= tri[site].offd[0][ 5].imag() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][ 9].real() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][ 9].imag() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][13].real() * ppsi[ 5].real(); " << endl;
	  sprg << "     cchi[ 3].real() += tri[site].offd[0][13].imag() * ppsi[ 5].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 3].imag()  = tri[site].diag[0][ 3].elem() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 3].real() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 3].imag() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 4].real() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 4].imag() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 5].real() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 5].imag() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][ 9].real() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() -= tri[site].offd[0][ 9].imag() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 3].imag() += tri[site].offd[0][13].real() * ppsi[ 5].imag(); " << endl;
	  sprg << "     cchi[ 3].imag() -= tri[site].offd[0][13].imag() * ppsi[ 5].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 4].real()  = tri[site].diag[0][ 4].elem() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][ 6].real() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 4].real() -= tri[site].offd[0][ 6].imag() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][ 7].real() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 4].real() -= tri[site].offd[0][ 7].imag() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][ 8].real() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 4].real() -= tri[site].offd[0][ 8].imag() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][ 9].real() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 4].real() -= tri[site].offd[0][ 9].imag() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][14].real() * ppsi[ 5].real(); " << endl;
	  sprg << "     cchi[ 4].real() += tri[site].offd[0][14].imag() * ppsi[ 5].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 4].imag()  = tri[site].diag[0][ 4].elem() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 6].real() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 6].imag() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 7].real() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 7].imag() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 8].real() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 8].imag() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 9].real() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][ 9].imag() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 4].imag() += tri[site].offd[0][14].real() * ppsi[ 5].imag(); " << endl;
	  sprg << "     cchi[ 4].imag() -= tri[site].offd[0][14].imag() * ppsi[ 5].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 5].real()  = tri[site].diag[0][ 5].elem() * ppsi[ 5].real(); " << endl;
	  sprg << "     cchi[ 5].real() += tri[site].offd[0][10].real() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 5].real() -= tri[site].offd[0][10].imag() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 5].real() += tri[site].offd[0][11].real() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 5].real() -= tri[site].offd[0][11].imag() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 5].real() += tri[site].offd[0][12].real() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 5].real() -= tri[site].offd[0][12].imag() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 5].real() += tri[site].offd[0][13].real() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 5].real() -= tri[site].offd[0][13].imag() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 5].real() += tri[site].offd[0][14].real() * ppsi[ 4].real(); " << endl;
	  sprg << "     cchi[ 5].real() -= tri[site].offd[0][14].imag() * ppsi[ 4].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 5].imag()  = tri[site].diag[0][ 5].elem() * ppsi[ 5].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][10].real() * ppsi[ 0].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][10].imag() * ppsi[ 0].real(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][11].real() * ppsi[ 1].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][11].imag() * ppsi[ 1].real(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][12].real() * ppsi[ 2].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][12].imag() * ppsi[ 2].real(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][13].real() * ppsi[ 3].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][13].imag() * ppsi[ 3].real(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][14].real() * ppsi[ 4].imag(); " << endl;
	  sprg << "     cchi[ 5].imag() += tri[site].offd[0][14].imag() * ppsi[ 4].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 6].real()  = tri[site].diag[1][0].elem()  * ppsi[6].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][0].real()  * ppsi[7].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][0].imag()  * ppsi[7].imag(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][1].real()  * ppsi[8].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][1].imag()  * ppsi[8].imag(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][3].real()  * ppsi[9].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][3].imag()  * ppsi[9].imag(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][6].real()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][6].imag()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][10].real() * ppsi[11].real(); " << endl;
	  sprg << "     cchi[ 6].real() += tri[site].offd[1][10].imag() * ppsi[11].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 6].imag()  = tri[site].diag[1][0].elem()  * ppsi[6].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() += tri[site].offd[1][0].real()  * ppsi[7].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() -= tri[site].offd[1][0].imag()  * ppsi[7].real(); " << endl;
	  sprg << "     cchi[ 6].imag() += tri[site].offd[1][1].real()  * ppsi[8].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() -= tri[site].offd[1][1].imag()  * ppsi[8].real(); " << endl;
	  sprg << "     cchi[ 6].imag() += tri[site].offd[1][3].real()  * ppsi[9].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() -= tri[site].offd[1][3].imag()  * ppsi[9].real(); " << endl;
	  sprg << "     cchi[ 6].imag() += tri[site].offd[1][6].real()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() -= tri[site].offd[1][6].imag()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 6].imag() += tri[site].offd[1][10].real() * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[ 6].imag() -= tri[site].offd[1][10].imag() * ppsi[11].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 7].real()  = tri[site].diag[1][ 1].elem()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 0].real()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 7].real() -= tri[site].offd[1][ 0].imag()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 2].real()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 2].imag()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 4].real()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 4].imag()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 7].real()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][ 7].imag()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][11].real()  * ppsi[11].real(); " << endl;
	  sprg << "     cchi[ 7].real() += tri[site].offd[1][11].imag()  * ppsi[11].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 7].imag()  = tri[site].diag[1][ 1].elem()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][ 0].real()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][ 0].imag()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][ 2].real()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() -= tri[site].offd[1][ 2].imag()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][ 4].real()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() -= tri[site].offd[1][ 4].imag()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][ 7].real()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() -= tri[site].offd[1][ 7].imag()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 7].imag() += tri[site].offd[1][11].real()  * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[ 7].imag() -= tri[site].offd[1][11].imag()  * ppsi[11].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 8].real()  = tri[site].diag[1][ 2].elem()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][ 1].real()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 8].real() -= tri[site].offd[1][ 1].imag()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][ 2].real()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[ 8].real() -= tri[site].offd[1][ 2].imag()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][5].real()   * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][5].imag()   * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][8].real()   * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][8].imag()   * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][12].real()  * ppsi[11].real(); " << endl;
	  sprg << "     cchi[ 8].real() += tri[site].offd[1][12].imag()  * ppsi[11].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 8].imag() = tri[site].diag[1][ 2].elem()   * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][ 1].real()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][ 1].imag()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][ 2].real()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][ 2].imag()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][5].real()   * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() -= tri[site].offd[1][5].imag()   * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][8].real()   * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() -= tri[site].offd[1][8].imag()   * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 8].imag() += tri[site].offd[1][12].real()  * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[ 8].imag() -= tri[site].offd[1][12].imag()  * ppsi[11].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 9].real()  = tri[site].diag[1][ 3].elem()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][ 3].real()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 9].real() -= tri[site].offd[1][ 3].imag()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][ 4].real()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[ 9].real() -= tri[site].offd[1][ 4].imag()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][ 5].real()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[ 9].real() -= tri[site].offd[1][ 5].imag()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][ 9].real()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][ 9].imag()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][13].real()  * ppsi[11].real(); " << endl;
	  sprg << "     cchi[ 9].real() += tri[site].offd[1][13].imag()  * ppsi[11].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[ 9].imag()  = tri[site].diag[1][ 3].elem()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 3].real()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 3].imag()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 4].real()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 4].imag()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 5].real()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 5].imag()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][ 9].real()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() -= tri[site].offd[1][ 9].imag()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[ 9].imag() += tri[site].offd[1][13].real()  * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[ 9].imag() -= tri[site].offd[1][13].imag()  * ppsi[11].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[10].real()  = tri[site].diag[1][ 4].elem()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][ 6].real()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[10].real() -= tri[site].offd[1][ 6].imag()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][ 7].real()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[10].real() -= tri[site].offd[1][ 7].imag()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][ 8].real()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[10].real() -= tri[site].offd[1][ 8].imag()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][ 9].real()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[10].real() -= tri[site].offd[1][ 9].imag()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][14].real()  * ppsi[11].real(); " << endl;
	  sprg << "     cchi[10].real() += tri[site].offd[1][14].imag()  * ppsi[11].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[10].imag()  = tri[site].diag[1][ 4].elem()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 6].real()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 6].imag()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 7].real()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 7].imag()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 8].real()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 8].imag()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 9].real()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][ 9].imag()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[10].imag() += tri[site].offd[1][14].real()  * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[10].imag() -= tri[site].offd[1][14].imag()  * ppsi[11].real(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[11].real()  = tri[site].diag[1][ 5].elem()  * ppsi[11].real(); " << endl;
	  sprg << "     cchi[11].real() += tri[site].offd[1][10].real()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[11].real() -= tri[site].offd[1][10].imag()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[11].real() += tri[site].offd[1][11].real()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[11].real() -= tri[site].offd[1][11].imag()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[11].real() += tri[site].offd[1][12].real()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[11].real() -= tri[site].offd[1][12].imag()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[11].real() += tri[site].offd[1][13].real()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[11].real() -= tri[site].offd[1][13].imag()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[11].real() += tri[site].offd[1][14].real()  * ppsi[10].real(); " << endl;
	  sprg << "     cchi[11].real() -= tri[site].offd[1][14].imag()  * ppsi[10].imag(); " << endl;
	  sprg << " 	 " << endl;
	  sprg << "     cchi[11].imag()  = tri[site].diag[1][ 5].elem()  * ppsi[11].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][10].real()  * ppsi[ 6].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][10].imag()  * ppsi[ 6].real(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][11].real()  * ppsi[ 7].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][11].imag()  * ppsi[ 7].real(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][12].real()  * ppsi[ 8].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][12].imag()  * ppsi[ 8].real(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][13].real()  * ppsi[ 9].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][13].imag()  * ppsi[ 9].real(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][14].real()  * ppsi[10].imag(); " << endl;
	  sprg << "     cchi[11].imag() += tri[site].offd[1][14].imag()  * ppsi[10].real(); " << endl;

	  prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	  cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	} else {
	  if (!cacheLock( psi , cudaArgs )) { QDP_info("eval: could not cache psi");  break;   }
	  if (!cacheLock( chi , cudaArgs )) { QDP_info("eval: could not cache chi");  break;   }
	}

	if (!QDPJit::Instance()( strId , prg , cudaArgs.getDevPtr() , rb[cb].numSiteTable() , sharedLibEntry  , mapVolumes )) {
	  QDP_info("packClov call to cuda jitter failed");
	  break;
	}

	return;
      }
      QDP_error("Clover: applyJIT host not implemented");
    }






  } // Namespace 

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
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  template<typename T, typename U>
  void JitCloverTermT<T,U>::apply(T& chi, const T& psi, 
			    enum PlusMinus isign, int cb) const
  {
    START_CODE();
    
    if ( Ns != 4 ) {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << endl;
      QDP_abort(1);
    }

    JITCloverEnv::ApplyArgs<T> arg = { chi,psi,tri_id,cb };
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    applyJIT( &arg );

    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }


  namespace JITCloverEnv {
    template<typename R> 
    struct QUDAPackArgs { 
      int cb;
      multi1d<QUDAPackedClovSite<R> >& quda_array;
      int tri_id;
    };
    


    template<typename U>
    inline 
    void qudaPackJIT( QUDAPackArgs<U>* a)
    {
      typedef U REALT;

      int cb = a->cb;
      int num_sites = rb[cb].siteTable().size();
      int Ns2 = Ns/2;

      multi1d<QUDAPackedClovSite<U> >& quda_array = a->quda_array;
      int tri_id = a->tri_id;

      while (1) {

	static QDPJit::SharedLibEntry sharedLibEntry;
	static MapVolumes*              mapVolumes;
	static string                   strId;
	static string                   prg;

	const int nodeSites = QDP::Layout::sitesOnNode();

	void* quda_array_dev;

	if (!QDPCache::Instance().allocate_device_static( &quda_array_dev , nodeSites * sizeof(QUDAPackedClovSite<U>) )) {
	  QDP_debug("Pack Clover JIT: no GPU memory left for quda_array field");
	  break;
	}

	CudaMemcpy( quda_array_dev , (void*)quda_array.slice() , nodeSites * sizeof(QUDAPackedClovSite<U>) );

	QDPJitArgs cudaArgs;

	int argNum = cudaArgs.addInt( rb[cb].numSiteTable() );
	int argOrd = cudaArgs.addBool( rb[cb].hasOrderedRep() );
	int argStart = cudaArgs.addInt( rb[cb].start() );
	int argSubset = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( rb[cb].getId() ) );
	int argDestPtr = cudaArgs.addPtr( quda_array_dev );
	int argMiscPtr = cudaArgs.addPtr( QDPCache::Instance().getDevicePtr( tri_id ) );

	if (!mapVolumes) {
	  string strREALT;
	  getTypeString( strREALT , REALT(0) );

	  ostringstream osId;
	  osId << "packClov "  << strREALT;
	  strId = osId.str();
	  xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	  cout << "strId = " << strId << endl;
#endif

	  std::ostringstream sprg;

	  sprg << " typedef " << strREALT << " REALT;" << endl;
	  sprg << " const int Nc = " << Nc << ";" << endl;
	  sprg << " const int Ns2 = " << Ns2 << ";" << endl;
	  //      sprg << " const int idtab[15]={0,1,3,6,10,2,4,7,11,5,8,12,9,13,14};" << endl;

	  sprg << "	struct PrimitiveClovTriang" << endl;
	  sprg << "	{" << endl;
	  sprg << "   RScalar<REALT>   diag[2][2*Nc];" << endl;
	  sprg << "   RComplex<REALT>  offd[2][2*Nc*Nc-Nc];" << endl;
	  sprg << "	};" << endl;

	  sprg << "  struct QUDAPackedClovSite {" << endl;
	  sprg << "    REALT diag1[6];" << endl;
	  sprg << "    REALT offDiag1[15][2];" << endl;
	  sprg << "    REALT diag2[6];" << endl;
	  sprg << "    REALT offDiag2[15][2];" << endl;
	  sprg << "  };" << endl;

	  sprg << "#define tri ((PrimitiveClovTriang*)( " << cudaArgs.getCode(argMiscPtr) << " ))" << endl;
	  sprg << "#define quda_array ((QUDAPackedClovSite*)( " << cudaArgs.getCode(argDestPtr) << " ))" << endl;

	  sprg << "int site;" << endl;

	  sprg << "  if (" << cudaArgs.getCode(argOrd) << ") {" << endl;
	  sprg << "    site = " << cudaArgs.getCode(argStart) << ";" << endl;
	  sprg << "    site += blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	  sprg << "    if (site >= " << cudaArgs.getCode(argNum) << "+" <<  cudaArgs.getCode(argStart) << ") return;" << endl;
	  sprg << "  } else {" << endl;
	  sprg << "    int idx0 = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x; \n";
	  sprg << "    if (idx0 >= " << cudaArgs.getCode(argNum) << ") return;" << endl;
	  sprg << "    site = ((int*)" << cudaArgs.getCode(argSubset) << ")[idx0];" << endl;
	  sprg << "  }" << endl;

	  sprg << "	for(int i=0; i < 6; i++) { " << endl;
	  sprg << "	  quda_array[site].diag1[i] = tri[site].diag[0][i].elem();" << endl;
	  sprg << "	}" << endl;
	  sprg << "	int target_index=0;" << endl;
	  sprg << "	for(int col=0; col < Nc*Ns2-1; col++) { " << endl;
	  sprg << "	  for(int row=col+1; row < Nc*Ns2; row++) {" << endl;
	  sprg << "	    int source_index = row*(row-1)/2 + col;" << endl;
	  sprg << "	    quda_array[site].offDiag1[target_index][0] = tri[site].offd[0][source_index].real();" << endl;
	  sprg << "	    quda_array[site].offDiag1[target_index][1] = tri[site].offd[0][source_index].imag();" << endl;
	  sprg << "	    target_index++;" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	}" << endl;
	  sprg << "	for(int i=0; i < 6; i++) { " << endl;
	  sprg << "	  quda_array[site].diag2[i] = tri[site].diag[1][i].elem();" << endl;
	  sprg << "	}" << endl;
	  sprg << "	target_index=0;" << endl;
	  sprg << "	for(int col=0; col < Nc*Ns2-1; col++) { " << endl;
	  sprg << "	  for(int row=col+1; row < Nc*Ns2; row++) {" << endl;
	  sprg << "	    int source_index = row*(row-1)/2 + col;" << endl;
	  sprg << "	    quda_array[site].offDiag2[target_index][0] = tri[site].offd[1][source_index].real();" << endl;
	  sprg << "	    quda_array[site].offDiag2[target_index][1] = tri[site].offd[1][source_index].imag();" << endl;
	  sprg << "	    target_index++;" << endl;
	  sprg << "	  }" << endl;
	  sprg << "	}" << endl;

	  prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	  cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	}

	if (!QDPJit::Instance()( strId , prg , cudaArgs.getDevPtr() , rb[cb].numSiteTable() , sharedLibEntry  , mapVolumes )) {
	  QDP_info("packClov call to cuda jitter failed");
	  break;
	}

	CudaMemcpy( (void*)quda_array.slice() , quda_array_dev , nodeSites * sizeof(QUDAPackedClovSite<U>) );

	QDPCache::Instance().free_device_static( quda_array_dev );
	return;
      }
      QDP_error("packClov host not implemented");
    }


  }


  template<typename T, typename U>
  void JitCloverTermT<T,U>::packForQUDA(multi1d<QUDAPackedClovSite<typename WordType<T>::Type_t> >& quda_array, int cb) const
    {
      typedef typename WordType<T>::Type_t REALT;
      int num_sites = rb[cb].siteTable().size();

      JITCloverEnv::QUDAPackArgs<REALT> args = { cb, quda_array,tri_id };
      qudaPackJIT( &args );
      
    }  



  typedef JitCloverTermT<LatticeFermion, LatticeColorMatrix> JitCloverTerm;
  typedef JitCloverTermT<LatticeFermionF, LatticeColorMatrixF> JitCloverTermF;
  typedef JitCloverTermT<LatticeFermionD, LatticeColorMatrixD> JitCloverTermD;
} // End Namespace Chroma


#endif
