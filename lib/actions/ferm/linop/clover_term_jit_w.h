// -*- C++ -*-
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_jit_w_h__
#define __clover_term_jit_w_h__

//#warning "Using QDP-JIT clover term"

#include "state.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/clover_term_base_w.h"
#include "meas/glue/mesfield.h"



namespace QDP
{
  class PackForQUDATimer {
    double acc_time;
    PackForQUDATimer(): acc_time(0.0) {}
  public:
    static PackForQUDATimer& Instance() {
      static PackForQUDATimer singleton;
      return singleton;
    }

    double&        get()       { return acc_time; }
    const double&  get() const { return acc_time; }
  };


  template<typename T>
  struct PComp
  {
    typedef T Sub_t;
    enum { ThisSize = 2 };
    T comp[2];
  };

  template<class T>  struct PCompREG;

  template<typename T>
  struct PCompJIT: public BaseJIT<T,2>
  {
    template<class T1>
    PCompJIT& operator=( const PCompREG<T1>& rhs) {
      //std::cout << __PRETTY_FUNCTION__ << "\n";
      elem(0) = rhs.elem(0);
      elem(1) = rhs.elem(1);
      return *this;
    }


    inline       T elem(int i)       { return this->arrayF(i); }
  };

  template<class T>
  struct PCompREG 
  {
    T F[2];
    void setup(PCompJIT< typename JITType<T>::Type_t > rhs ) {
      F[0].setup( rhs.elem(0) );
      F[1].setup( rhs.elem(1) );
    }
    inline       T& elem(int i)       { return F[i]; }
    inline const T& elem(int i) const { return F[i]; }
  };

  template<class T> 
  struct ScalarType<PComp<T> >
  {
    typedef PComp<typename ScalarType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct ScalarType<PCompJIT<T> >
  {
    typedef PCompJIT<typename ScalarType<T>::Type_t>  Type_t;
  };
  
  template<class T> 
  struct JITType<PComp<T> >
  {
    typedef PCompJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct JITType<PCompREG<T> >
  {
    typedef PCompJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct REGType<PCompJIT<T> >
  {
    typedef PCompREG<typename REGType<T>::Type_t>  Type_t;
  };

  template<class T>
  struct WordType<PComp<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };

  template<class T>
  struct WordType<PCompJIT<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };






  template<typename T>
  struct PTriDia
  {
    typedef T Sub_t;
    enum { ThisSize = 2*Nc };
    T diag[2*Nc];
  };



  template<class T> struct PTriDiaREG;

  template<typename T>
  struct PTriDiaJIT: public BaseJIT<T,2*Nc>
  {
    template<class T1>
    PTriDiaJIT& operator=( const PTriDiaREG<T1>& rhs) {
      //std::cout << __PRETTY_FUNCTION__ << "\n";
      for ( int i = 0 ; i < 2 * Nc ; i++ )
	elem(i) = rhs.elem(i);
      return *this;
    }

    inline       T elem(int i)       { return this->arrayF(i); }
  };

  template<class T>
  struct PTriDiaREG 
  {
    T F[2*Nc];
    void setup( PTriDiaJIT< typename JITType<T>::Type_t > rhs ) {
      for (int i=0;i<2*Nc;++i)
	F[i].setup( rhs.elem(i) );
    }
    inline       T& elem(int i)       { return F[i]; }
    inline const T& elem(int i) const { return F[i]; }
  };


  template<class T> 
  struct ScalarType<PTriDia<T> >
  {
    typedef PTriDia<typename ScalarType<T>::Type_t>  Type_t;
  };
  
  template<class T> 
  struct ScalarType<PTriDiaJIT<T> >
  {
    typedef PTriDiaJIT<typename ScalarType<T>::Type_t>  Type_t;
  };

  
  template<class T> 
  struct JITType<PTriDia<T> >
  {
    typedef PTriDiaJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct JITType<PTriDiaREG<T> >
  {
    typedef PTriDiaJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct REGType<PTriDiaJIT<T> >
  {
    typedef PTriDiaREG<typename REGType<T>::Type_t>  Type_t;
  };

  template<class T>
  struct WordType<PTriDia<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };

  template<class T>
  struct WordType<PTriDiaJIT<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };






  template<typename T>
  struct PTriOff
  {
    typedef T Sub_t;
    enum { ThisSize = 2*Nc*Nc-Nc };
    T offd[2*Nc*Nc-Nc];
  };

  template<class T> struct PTriOffREG;

  template<typename T>
  struct PTriOffJIT: public BaseJIT<T,2*Nc*Nc-Nc>
  {
    template<class T1>
    PTriOffJIT& operator=( const PTriOffREG<T1>& rhs) {
      //std::cout << __PRETTY_FUNCTION__ << "\n";
      for ( int i = 0 ; i < 2*Nc*Nc-Nc ; i++ )
	elem(i) = rhs.elem(i);
      return *this;
    }

    inline       T elem(int i)       { return this->arrayF(i); }
  };

  template<class T>
  struct PTriOffREG 
  {
    T F[2*Nc*Nc-Nc];
    void setup( PTriOffJIT< typename JITType<T>::Type_t > rhs ) {
      for (int i=0;i<2*Nc*Nc-Nc;++i)
	F[i].setup( rhs.elem(i) );
    }
    inline       T& elem(int i)       { return F[i]; }
    inline const T& elem(int i) const { return F[i]; }
  };


  template<class T> 
  struct ScalarType<PTriOff<T> >
  {
    typedef PTriOff<typename ScalarType<T>::Type_t>  Type_t;
  };
  
  template<class T> 
  struct ScalarType<PTriOffJIT<T> >
  {
    typedef PTriOffJIT<typename ScalarType<T>::Type_t>  Type_t;
  };


  template<class T> 
  struct JITType<PTriOff<T> >
  {
    typedef PTriOffJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct JITType<PTriOffREG<T> >
  {
    typedef PTriOffJIT<typename JITType<T>::Type_t>  Type_t;
  };

  template<class T> 
  struct REGType<PTriOffJIT<T> >
  {
    typedef PTriOffREG<typename REGType<T>::Type_t>  Type_t;
  };

  template<class T>
  struct WordType<PTriOff<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };

  template<class T>
  struct WordType<PTriOffJIT<T> > 
  {
    typedef typename WordType<T>::Type_t  Type_t;
  };



  template<class T>
  struct LeafFunctor<PComp<T>, PrintTag>
  {
    typedef int Type_t;
    static int apply(const PrintTag &f)
    { 
      f.os_m << "PComp<";
      LeafFunctor<T,PrintTag>::apply(f);
      f.os_m << ">"; 
      return 0;
    }
  };

  template<class T>
  struct LeafFunctor<PTriDia<T>, PrintTag>
  {
    typedef int Type_t;
    static int apply(const PrintTag &f)
    { 
      f.os_m << "PTriDia<";
      LeafFunctor<T,PrintTag>::apply(f);
      f.os_m << ">"; 
      return 0;
    }
  };

  template<class T>
  struct LeafFunctor<PTriOff<T>, PrintTag>
  {
    typedef int Type_t;
    static int apply(const PrintTag &f)
    { 
      f.os_m << "PTriOff<";
      LeafFunctor<T,PrintTag>::apply(f);
      f.os_m << ">"; 
      return 0;
    }
  };
} // QDP


#if defined (QDP_BACKEND_AVX)
#define WORD WordVec
#else
#define WORD Word
#endif


namespace Chroma 
{ 
  template<typename R>
  struct QUDAPackedClovSite {
    R diag1[6];
    R offDiag1[15][2];
    R diag2[6];
    R offDiag2[15][2];
  };


  template<typename T, typename U>
  class JITCloverTermT : public CloverTermBase<T, U>
  {
  public:
    // Typedefs to save typing
    typedef typename WordType<T>::Type_t REALT;

    typedef OLattice< PScalar< PScalar< RScalar< WORD< REALT> > > > > LatticeREAL;
    typedef OScalar<  PScalar< PScalar< RScalar< Word< REALT> > > > > RealT;

    //! Empty constructor. Must use create later
    JITCloverTermT();

    //! No real need for cleanup here
    ~JITCloverTermT() {}

    //! Creation routine
    void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
		const CloverFermActParams& param_);

    virtual void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
			const CloverFermActParams& param_,
			const JITCloverTermT<T,U>& from_);

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
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     */
    void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const;


    void applySite(T& chi, const T& psi, enum PlusMinus isign, int site) const;

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(U& B, int mat, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T, multi1d<U>, multi1d<U> >& getFermBC() const {return *fbc;}

    //! PACK UP the Clover term for QUDA library:
    void packForQUDA(multi1d<QUDAPackedClovSite<REALT> >& quda_pack, int cb) const; 

    int getDiaId() const { return tri_dia.getId(); }
    int getOffId() const { return tri_off.getId(); }

      
  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<U>& f, const RealT& diag_mass);

    //! Invert the clover term on cb
    //void chlclovms(LatticeREAL& log_diag, int cb);
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

    OLattice<PComp<PTriDia<RScalar <WORD<REALT> > > > >  tri_dia;
    OLattice<PComp<PTriOff<RComplex<WORD<REALT> > > > >  tri_off;    
  };

#undef WORD

  

   // Empty constructor. Must use create later
  template<typename T, typename U>
  JITCloverTermT<T,U>::JITCloverTermT() {}

  // Now copy
  template<typename T, typename U>
  void JITCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_,
				   const JITCloverTermT<T,U>& from)
  {
    START_CODE();

    //std::cout << "PTX Clover create from other "  << (void*)this << "\n";

    u.resize(Nd);

    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "JITCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }
    
    {
      RealT ff = param.anisoParam.anisoP ? Real(1) / param.anisoParam.xi_0 : Real(1);
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
      RealT ff = param.anisoParam.anisoP ? param.anisoParam.nu / param.anisoParam.xi_0 : Real(1);
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
    
    tri_dia = from.tri_dia;
    tri_off = from.tri_off;

    END_CODE();  
  }


  //! Creation routine
  template<typename T, typename U>
  void JITCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_)
  {
    START_CODE();

    //std::cout << "PTX Clover create "  << (void*)this << "\n";
   
    u.resize(Nd);
    
    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "JITCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }

    {
      RealT ff = param.anisoParam.anisoP ? Real(1) / param.anisoParam.xi_0 : Real(1);
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
      RealT ff = param.anisoParam.anisoP ? param.anisoParam.nu / param.anisoParam.xi_0 : Real(1);
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

  template<typename RealT,typename U,typename X,typename Y>
  void function_make_clov_exec(JitFunction& function, 
			       const RealT& diag_mass, 
			       const U& f0,
			       const U& f1,
			       const U& f2,
			       const U& f3,
			       const U& f4,
			       const U& f5,
			       X& tri_dia,
			       Y& tri_off)
  {
#ifdef QDP_DEEP_LOG
    function.type_W = typeid(REAL).name();
    //function.set_dest_id( tri_dia.getId() );
    function.set_dest_id( tri_off.getId() );
    function.set_is_lat(true);
#endif
    
    AddressLeaf addr_leaf(all);

    forEach(diag_mass, addr_leaf, NullCombine());
    forEach(f0, addr_leaf, NullCombine());
    forEach(f1, addr_leaf, NullCombine());
    forEach(f2, addr_leaf, NullCombine());
    forEach(f3, addr_leaf, NullCombine());
    forEach(f4, addr_leaf, NullCombine());
    forEach(f5, addr_leaf, NullCombine());
    forEach(tri_dia, addr_leaf, NullCombine());
    forEach(tri_off, addr_leaf, NullCombine());

    int th_count = Layout::sitesOnNode();

    WorkgroupGuardExec workgroupGuardExec(th_count);

    std::vector<QDPCache::ArgKey> ids;
    workgroupGuardExec.check(ids);
    //ids.push_back( s.getIdSiteTable() );
    for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
      ids.push_back( addr_leaf.ids[i] );
    jit_launch(function,th_count,ids);
  }



  template<typename RealT,typename U,typename X,typename Y>
  void function_make_clov_build(JitFunction& function,
				const RealT& diag_mass, 
				const U& f0,
				const U& f1,
				const U& f2,
				const U& f3,
				const U& f4,
				const U& f5,
				const X& tri_dia,
				const Y& tri_off)
  {
    //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

    typedef typename WordType<RealT>::Type_t REALT;

    llvm_start_new_function("make_clov",__PRETTY_FUNCTION__ );

    WorkgroupGuard workgroupGuard;

    ParamLeafScalar param_leaf;
    
    typedef typename LeafFunctor<RealT, ParamLeafScalar>::Type_t  RealTJIT;
    RealTJIT diag_mass_jit(forEach(diag_mass, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<U, ParamLeafScalar>::Type_t  UJIT;
    UJIT f0_jit(forEach(f0, param_leaf, TreeCombine()));
    UJIT f1_jit(forEach(f1, param_leaf, TreeCombine()));
    UJIT f2_jit(forEach(f2, param_leaf, TreeCombine()));
    UJIT f3_jit(forEach(f3, param_leaf, TreeCombine()));
    UJIT f4_jit(forEach(f4, param_leaf, TreeCombine()));
    UJIT f5_jit(forEach(f5, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<X, ParamLeafScalar>::Type_t  XJIT;
    XJIT tri_dia_jit(forEach(tri_dia, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<Y, ParamLeafScalar>::Type_t  YJIT;
    YJIT tri_off_jit(forEach(tri_off, param_leaf, TreeCombine()));

    llvm::Value* r_idx = llvm_thread_idx();

    workgroupGuard.check(r_idx);

    auto f0_j = f0_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto f1_j = f1_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto f2_j = f2_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto f3_j = f3_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto f4_j = f4_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto f5_j = f5_jit.elem(JitDeviceLayout::Coalesced , r_idx );

    auto tri_dia_j = tri_dia_jit.elem(JitDeviceLayout::Coalesced , r_idx );
    auto tri_off_j = tri_off_jit.elem(JitDeviceLayout::Coalesced , r_idx );

    typename REGType< typename RealTJIT::Subtype_t >::Type_t diag_mass_reg;
    //diag_mass_reg.setup( diag_mass_jit.elem() );
    diag_mass_reg.setup_value( diag_mass_jit.elem() );


    for(int jj = 0; jj < 2; jj++) {
      for(int ii = 0; ii < 2*Nc; ii++) {
	tri_dia_j.elem(jj).elem(ii) = diag_mass_reg.elem().elem();
	//tri[site].diag[jj][ii] = diag_mass.elem().elem().elem();
      }
    }


    RComplexREG<WordREG<REALT> > E_minus;
    RComplexREG<WordREG<REALT> > B_minus;
    RComplexREG<WordREG<REALT> > ctmp_0;
    RComplexREG<WordREG<REALT> > ctmp_1;
    RScalarREG<WordREG<REALT> > rtmp_0;
    RScalarREG<WordREG<REALT> > rtmp_1;


    for(int i = 0; i < Nc; ++i) {
      ctmp_0 = f5_j.elem().elem(i,i);
      ctmp_0 -= f0_j.elem().elem(i,i);
      rtmp_0 = imag(ctmp_0);
      tri_dia_j.elem(0).elem(i) += rtmp_0;
	  
      tri_dia_j.elem(0).elem(i+Nc) -= rtmp_0;
	  
      ctmp_1 = f5_j.elem().elem(i,i);
      ctmp_1 += f0_j.elem().elem(i,i);
      rtmp_1 = imag(ctmp_1);
      tri_dia_j.elem(1).elem(i) -= rtmp_1;
	  
      tri_dia_j.elem(1).elem(i+Nc) += rtmp_1;
    }

    for(int i = 1; i < Nc; ++i) {
      for(int j = 0; j < i; ++j) {
	    
	int elem_ij  = i*(i-1)/2 + j;
	int elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;
	    
	ctmp_0 = f0_j.elem().elem(i,j);
	ctmp_0 -= f5_j.elem().elem(i,j);
	tri_off_j.elem(0).elem(elem_ij) = timesI(ctmp_0);
	    
	zero_rep( tri_off_j.elem(0).elem(elem_tmp) );
	tri_off_j.elem(0).elem(elem_tmp) -= tri_off_j.elem(0).elem(elem_ij);// * -1.0;
	    
	ctmp_1 = f5_j.elem().elem(i,j);
	ctmp_1 += f0_j.elem().elem(i,j);
	tri_off_j.elem(1).elem(elem_ij) = timesI(ctmp_1);
	    
	zero_rep( tri_off_j.elem(1).elem(elem_tmp) );
	tri_off_j.elem(1).elem(elem_tmp) -= tri_off_j.elem(1).elem(elem_ij);
      }
    }

    for(int i = 0; i < Nc; ++i) {
      for(int j = 0; j < Nc; ++j) {
	    
	int elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;
	    
	//E_minus = timesI(f2_j.elem().elem(i,j));
	E_minus = f2_j.elem().elem(i,j);
	E_minus = timesI( E_minus );

	E_minus += f4_j.elem().elem(i,j);
	    
	//B_minus = timesI(f3_j.elem().elem(i,j));
	B_minus = f3_j.elem().elem(i,j);
	B_minus = timesI( B_minus );

	B_minus -= f1_j.elem().elem(i,j);
	    
	tri_off_j.elem(0).elem(elem_ij) = B_minus - E_minus;
	    
	tri_off_j.elem(1).elem(elem_ij) = E_minus + B_minus;
      }
    }

    //    std::cout << __PRETTY_FUNCTION__ << ": leaving\n";

    jit_get_function(function);
  }



  
  /* This now just sets up and dispatches... */
  template<typename T, typename U>
  void JITCloverTermT<T,U>::makeClov(const multi1d<U>& f, const RealT& diag_mass)
  {
    START_CODE();
    
    if ( Nd != 4 ){
      QDPIO::cerr << __func__ << ": expecting Nd==4" << std::endl;
      QDP_abort(1);
    }
    
    if ( Ns != 4 ){
      QDPIO::cerr << __func__ << ": expecting Ns==4" << std::endl;
      QDP_abort(1);
    }
  
    U f0 = f[0] * getCloverCoeff(0,1);
    U f1 = f[1] * getCloverCoeff(0,2);
    U f2 = f[2] * getCloverCoeff(0,3);
    U f3 = f[3] * getCloverCoeff(1,2);
    U f4 = f[4] * getCloverCoeff(1,3);
    U f5 = f[5] * getCloverCoeff(2,3);    

    
    //QDPIO::cout << "PTX Clover make "  << (void*)this << "\n";
    //std::cout << "PTX Clover make "  << (void*)this << "\n";
    static JitFunction function;

    if (function.empty())
      function_make_clov_build(function, diag_mass, f0,f1,f2,f3,f4,f5, tri_dia , tri_off );

    // Execute the function
    function_make_clov_exec(function, diag_mass, f0,f1,f2,f3,f4,f5,tri_dia, tri_off);

    END_CODE();
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  template<typename T, typename U>
  void JITCloverTermT<T,U>::choles(int cb)
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
  Double JITCloverTermT<T,U>::cholesDet(int cb) const
  {
    START_CODE();

    if( choles_done[cb] == false ) 
      {
	QDPIO::cout << __func__ << ": Error: you have not done the Cholesky.on this operator on this subset" << std::endl;
	QDPIO::cout << "You sure you should not be asking invclov?" << std::endl;
	QDP_abort(1);
      }

    LatticeREAL ff=tr_log_diag_;


    END_CODE();

    // Need to thread generic sums in QDP++?
    // Need to thread generic norm2() in QDP++?
    return sum(tr_log_diag_, rb[cb]);
  }


  template<typename T,typename X,typename Y>
  void function_ldagdlinv_exec( JitFunction& function,
				T& tr_log_diag,
				X& tri_dia,
				Y& tri_off,
				const Subset& s)
  {
#ifdef QDP_DEEP_LOG
    function.type_W = typeid(REAL).name();
    function.set_dest_id( tr_log_diag.getId() );
    function.set_is_lat(true);
#endif
    
    AddressLeaf addr_leaf(s);

    forEach(tr_log_diag, addr_leaf, NullCombine());
    forEach(tri_dia, addr_leaf, NullCombine());
    forEach(tri_off, addr_leaf, NullCombine());

    int th_count = s.numSiteTable();

    WorkgroupGuardExec workgroupGuardExec(th_count);

    std::vector<QDPCache::ArgKey> ids;
    workgroupGuardExec.check(ids);
    ids.push_back( s.getIdSiteTable() );
    for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
      ids.push_back( addr_leaf.ids[i] );
    jit_launch(function,th_count,ids);
  }





  template<typename U,typename T,typename X,typename Y>
  void function_ldagdlinv_build(JitFunction& function,
				const T& tr_log_diag,
				const X& tri_dia,
				const Y& tri_off,
				const Subset& s)
  {
    typedef typename WordType<U>::Type_t REALT;

    //std::cout << __PRETTY_FUNCTION__ << " entering\n";

    llvm_start_new_function("ldagdlinv",__PRETTY_FUNCTION__);

    WorkgroupGuard workgroupGuard;
    ParamRef p_site_table = llvm_add_param<int*>();

    ParamLeafScalar param_leaf;

    typedef typename LeafFunctor<T, ParamLeafScalar>::Type_t  TJIT;
    TJIT tr_log_diag_jit(forEach(tr_log_diag, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<X, ParamLeafScalar>::Type_t  XJIT;
    XJIT tri_dia_jit(forEach(tri_dia, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<Y, ParamLeafScalar>::Type_t  YJIT;
    YJIT tri_off_jit(forEach(tri_off, param_leaf, TreeCombine()));

    llvm::Value* r_idx_thread = llvm_thread_idx();

    workgroupGuard.check(r_idx_thread);

    llvm::Value* r_idx = llvm_array_type_indirection<int>( p_site_table , r_idx_thread );

    auto tr_log_diag_j = tr_log_diag_jit.elem(JitDeviceLayout::Coalesced,r_idx);
    auto tri_dia_j     = tri_dia_jit.elem(JitDeviceLayout::Coalesced,r_idx);
    auto tri_off_j     = tri_off_jit.elem(JitDeviceLayout::Coalesced,r_idx);

    typename REGType< typename XJIT::Subtype_t >::Type_t tri_dia_r;
    typename REGType< typename YJIT::Subtype_t >::Type_t tri_off_r;

    tri_dia_r.setup( tri_dia_j );
    tri_off_r.setup( tri_off_j );


    RScalarREG<WordREG<REALT> > zip;
    zero_rep(zip);
    int N = 2*Nc;
      
    //int site_neg_logdet=0;
  
    for(int block=0; block < 2; block++) {
	  
      RScalarREG<WordREG<REALT> >   inv_d[6] ;
      RComplexREG<WordREG<REALT> >  inv_offd[15] ;
      RComplexREG<WordREG<REALT> >  v[6] ;
      RScalarREG<WordREG<REALT> >   diag_g[6] ;

      for(int i=0; i < N; i++) {
	inv_d[i] = tri_dia_r.elem(block).elem(i);
      }

      for(int i=0; i < 15; i++) { 
	inv_offd[i] = tri_off_r.elem(block).elem(i);
      }


      for(int j=0; j < N; ++j) { 
	    
	for(int i=0; i < j; i++) { 
	  int elem_ji = j*(j-1)/2 + i;
	      
	  RComplexREG<WordREG<REALT> > A_ii = cmplx( inv_d[i], zip );
	  v[i] = A_ii*adj(inv_offd[elem_ji]);
	}
	    

	v[j] = cmplx(inv_d[j],zip);
	    
	for(int k=0; k < j; k++) { 
	  int elem_jk = j*(j-1)/2 + k;
	  v[j] -= inv_offd[elem_jk]*v[k];
	}
	    
	inv_d[j] = real( v[j] );
	    
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
      RScalarREG<WordREG<REALT> > one(1.0);
      //one.elem() = (REALT)1;
	  
      for(int i=0; i < N; i++) { 
	diag_g[i] = one/inv_d[i];
      
	// Compute the trace log
	// NB we are always doing trace log | A | 
	// (because we are always working with actually A^\dagger A
	//  even in one flavour case where we square root)
	tr_log_diag_j.elem().elem() += log(fabs(inv_d[i]));
	// However, it is worth counting just the no of negative logdets
	// on site
#if 0
	if( inv_d[i].elem() < 0 ) { 
	  site_neg_logdet++;
	}
#endif
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
      RComplexREG<WordREG<REALT> > sum;
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
	tri_dia_j.elem(block).elem(i) = inv_d[i];
      }
      for(int i=0; i < 15; i++) { 
	tri_off_j.elem(block).elem(i) = inv_offd[i];
      }
    }

    //    std::cout << __PRETTY_FUNCTION__ << " leaving\n";

    jit_get_function(function);
  }





  /*! An LDL^\dag decomposition and inversion? */
  template<typename T, typename U>
  void JITCloverTermT<T,U>::ldagdlinv(LatticeREAL& tr_log_diag, int cb)
  {
    START_CODE();

    if ( 2*Nc < 3 )
      {
	QDPIO::cerr << __func__ << ": Matrix is too small" << std::endl;
	QDP_abort(1);
      }

    // Zero trace log
    tr_log_diag[rb[cb]] = zero;

    //QDPIO::cout << "PTX Clover ldagdlinv " << (void*)this << "\n";
    //std::cout << "PTX Clover ldagdlinv " << (void*)this << "\n";
    static JitFunction function;

    if (function.empty())
      function_ldagdlinv_build<U>(function, tr_log_diag, tri_dia, tri_off, rb[cb] );

    // Execute the function
    function_ldagdlinv_exec(function, tr_log_diag, tri_dia, tri_off, rb[cb] );

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


  template<typename U,typename X,typename Y>
  void function_triacntr_exec( JitFunction& function,
			       U& B,
			       const X& tri_dia,
			       const Y& tri_off,
			       int mat,
			       const Subset& s)
  {
#ifdef QDP_DEEP_LOG
    function.type_W = typeid(REAL).name();
    function.set_dest_id( B.getId() );
    function.set_is_lat(true);
#endif
    
    AddressLeaf addr_leaf(s);

    forEach(B, addr_leaf, NullCombine());
    forEach(tri_dia, addr_leaf, NullCombine());
    forEach(tri_off, addr_leaf, NullCombine());

    int th_count = s.numSiteTable();

    WorkgroupGuardExec workgroupGuardExec(th_count);

    JitParam jit_mat( QDP_get_global_cache().addJitParamInt( mat ) );

    std::vector<QDPCache::ArgKey> ids;
    workgroupGuardExec.check(ids);
    ids.push_back( s.getIdSiteTable() );
    ids.push_back( jit_mat.get_id() );
    for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
      ids.push_back( addr_leaf.ids[i] );
    jit_launch(function,th_count,ids);
  }




  template<typename U,typename X,typename Y>
  void function_triacntr_build( JitFunction& function,
				const U& B,
				const X& tri_dia,
				const Y& tri_off,
				int mat,
				const Subset& s)
  {
    //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

    typedef typename WordType<U>::Type_t REALT;

    llvm_start_new_function( "triacntr" , __PRETTY_FUNCTION__ );

    WorkgroupGuard workgroupGuard;
    ParamRef p_site_table = llvm_add_param<int*>();

    ParamRef p_mat    = llvm_add_param<int>();

    ParamLeafScalar param_leaf;

    typedef typename LeafFunctor<U, ParamLeafScalar>::Type_t  UJIT;
    UJIT B_jit(forEach(B, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<X, ParamLeafScalar>::Type_t  XJIT;
    XJIT tri_dia_jit(forEach(tri_dia, param_leaf, TreeCombine()));

    typedef typename LeafFunctor<Y, ParamLeafScalar>::Type_t  YJIT;
    YJIT tri_off_jit(forEach(tri_off, param_leaf, TreeCombine()));

    llvm::Value* r_idx_thread = llvm_thread_idx();

    workgroupGuard.check(r_idx_thread);

    llvm::Value* r_idx = llvm_array_type_indirection<int>( p_site_table , r_idx_thread );

    llvm::Value * r_mat    = llvm_derefParam( p_mat );

    auto B_j = B_jit.elem(JitDeviceLayout::Coalesced,r_idx);
    auto tri_dia_j = tri_dia_jit.elem(JitDeviceLayout::Coalesced,r_idx);
    auto tri_off_j = tri_off_jit.elem(JitDeviceLayout::Coalesced,r_idx);

    typename REGType< typename XJIT::Subtype_t >::Type_t tri_dia_r;
    typename REGType< typename YJIT::Subtype_t >::Type_t tri_off_r;

    tri_dia_r.setup( tri_dia_j );
    tri_off_r.setup( tri_off_j );

    JitSwitch sw(r_mat);
    {
      /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
      /*#               0  1  0  0 */
      /*#               0  0  1  0 */
      /*#               0  0  0  1 */
      /*# From diagonal part */
      sw.case_begin(0);
      {
	RComplexREG<WordREG<REALT> > lctmp0;
	RScalarREG< WordREG<REALT> > lr_zero0;
	RScalarREG< WordREG<REALT> > lrtmp0;
    
	zero_rep(lr_zero0);
	    
	for(int i0 = 0; i0 < Nc; ++i0)
	  {
	    lrtmp0 = tri_dia_r.elem(0).elem(i0);
	    lrtmp0 += tri_dia_r.elem(0).elem(i0+Nc);
	    lrtmp0 += tri_dia_r.elem(1).elem(i0);
	    lrtmp0 += tri_dia_r.elem(1).elem(i0+Nc);
	    B_j.elem().elem(i0,i0) = cmplx(lrtmp0,lr_zero0);
	  }
	    
	/*# From lower triangular portion */
	int elem_ij0 = 0;
	for(int i0 = 1; i0 < Nc; ++i0) {
	      
	  int elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc;
	      
	  for(int j0 = 0; j0 < i0; ++j0) {
		
	    lctmp0 = tri_off_r.elem(0).elem(elem_ij0);
	    lctmp0 += tri_off_r.elem(0).elem(elem_ijb0);
	    lctmp0 += tri_off_r.elem(1).elem(elem_ij0);
	    lctmp0 += tri_off_r.elem(1).elem(elem_ijb0);
		
	    B_j.elem().elem(j0,i0) = lctmp0;
	    B_j.elem().elem(i0,j0) = adj(lctmp0);
		
	    elem_ij0++;
	    elem_ijb0++;
	  }
	}
      }
      sw.case_end();


      /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */
      /*#               0  i  0  0 */
      /*#               0  0 -i  0 */
      /*#               0  0  0  i */
      /*# From diagonal part */
      sw.case_begin( 3 );
      {
	RComplexREG<WordREG<REALT> > lctmp3;
	RScalarREG<WordREG<REALT> > lr_zero3;
	RScalarREG<WordREG<REALT> > lrtmp3;
	    
	lr_zero3 = 0;
	    
	for(int i3 = 0; i3 < Nc; ++i3) {
	      
	  lrtmp3 = tri_dia_r.elem(0).elem(i3+Nc);
	  lrtmp3 -= tri_dia_r.elem(0).elem(i3);
	  lrtmp3 -= tri_dia_r.elem(1).elem(i3);
	  lrtmp3 += tri_dia_r.elem(1).elem(i3+Nc);
	  B_j.elem().elem(i3,i3) = cmplx(lr_zero3,lrtmp3);
	}
	    
	/*# From lower triangular portion */
	int elem_ij3 = 0;
	for(int i3 = 1; i3 < Nc; ++i3) {
	      
	  int elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc;
	      
	  for(int j3 = 0; j3 < i3; ++j3) {
		
	    lctmp3 = tri_off_r.elem(0).elem(elem_ijb3);
	    lctmp3 -= tri_off_r.elem(0).elem(elem_ij3);
	    lctmp3 -= tri_off_r.elem(1).elem(elem_ij3);
	    lctmp3 += tri_off_r.elem(1).elem(elem_ijb3);
		
	    B_j.elem().elem(j3,i3) = timesI(adj(lctmp3));
	    B_j.elem().elem(i3,j3) = timesI(lctmp3);
		
	    elem_ij3++;
	    elem_ijb3++;
	  }
	}
      }
      sw.case_end();

      /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */
      /*#               1  0  0  0 */
      /*#               0  0  0 -1 */
      /*#               0  0  1  0 */
      sw.case_begin( 5 );
      {
	RComplexREG<WordREG<REALT> > lctmp5;
	RScalarREG<WordREG<REALT> > lrtmp5;
	    
	for(int i5 = 0; i5 < Nc; ++i5) {
	      
	  int elem_ij5 = (i5+Nc)*(i5+Nc-1)/2;
	      
	  for(int j5 = 0; j5 < Nc; ++j5) {
		
	    int elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5;
		
	    lctmp5 = adj(tri_off_r.elem(0).elem(elem_ji5));
	    lctmp5 -= tri_off_r.elem(0).elem(elem_ij5);
	    lctmp5 += adj(tri_off_r.elem(1).elem(elem_ji5));
	    lctmp5 -= tri_off_r.elem(1).elem(elem_ij5);
		
	    B_j.elem().elem(i5,j5) = lctmp5;
		
	    elem_ij5++;
	  }
	}
      }
      sw.case_end();

      /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */
      /*#              -i  0  0  0 */
      /*#               0  0  0 -i */
      /*#               0  0 -i  0 */
      sw.case_begin( 6 );
      {
	RComplexREG<WordREG<REALT> > lctmp6;
	RScalarREG<WordREG<REALT> > lrtmp6;
	    
	for(int i6 = 0; i6 < Nc; ++i6) {
	      
	  int elem_ij6 = (i6+Nc)*(i6+Nc-1)/2;
	      
	  for(int j6 = 0; j6 < Nc; ++j6) {
		
	    int elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6;
	
	    lctmp6 = adj(tri_off_r.elem(0).elem(elem_ji6));
	    lctmp6 += tri_off_r.elem(0).elem(elem_ij6);
	    lctmp6 += adj(tri_off_r.elem(1).elem(elem_ji6));
	    lctmp6 += tri_off_r.elem(1).elem(elem_ij6);
		
	    B_j.elem().elem(i6,j6) = timesMinusI(lctmp6);
		
	    elem_ij6++;
	  }
	}
      }
      sw.case_end();

      /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */
      /*#               i  0  0  0 */
      /*#               0  0  0 -i */
      /*#               0  0 -i  0 */
      sw.case_begin( 9 );
      {
	RComplexREG<WordREG<REALT> > lctmp9;
	RScalarREG<WordREG<REALT> > lrtmp9;
	    
	for(int i9 = 0; i9 < Nc; ++i9) {
	      
	  int elem_ij9 = (i9+Nc)*(i9+Nc-1)/2;
	      
	  for(int j9 = 0; j9 < Nc; ++j9) {
		
	    int elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9;
	
	    lctmp9 = adj(tri_off_r.elem(0).elem(elem_ji9));
	    lctmp9 += tri_off_r.elem(0).elem(elem_ij9);
	    lctmp9 -= adj(tri_off_r.elem(1).elem(elem_ji9));
	    lctmp9 -= tri_off_r.elem(1).elem(elem_ij9);
		
	    B_j.elem().elem(i9,j9) = timesI(lctmp9);
		
	    elem_ij9++;
	  }
	}
      }
      sw.case_end();


      /*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */
      /*#               1  0  0  0 */
      /*#               0  0  0  1 */
      /*#               0  0 -1  0 */
      sw.case_begin( 10 );
      {
	RComplexREG<WordREG<REALT> > lctmp10;
	RScalarREG<WordREG<REALT> > lrtmp10;
	    
	for(int i10 = 0; i10 < Nc; ++i10) {
	      
	  int elem_ij10 = (i10+Nc)*(i10+Nc-1)/2;
	      
	  for(int j10 = 0; j10 < Nc; ++j10) {
		
	    int elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10;
	
	    lctmp10 = adj(tri_off_r.elem(0).elem(elem_ji10));
	    lctmp10 -= tri_off_r.elem(0).elem(elem_ij10);
	    lctmp10 -= adj(tri_off_r.elem(1).elem(elem_ji10));
	    lctmp10 += tri_off_r.elem(1).elem(elem_ij10);
		
	    B_j.elem().elem(i10,j10) = lctmp10;
		
	    elem_ij10++;
	  }
	}
      }
      sw.case_end();


      /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */
      /*#               0 -i  0  0 */
      /*#               0  0 -i  0 */
      /*#               0  0  0  i */
      /*# From diagonal part */
      sw.case_begin( 12 );
      {
	RComplexREG<WordREG<REALT> > lctmp12;
	RScalarREG<WordREG<REALT> > lr_zero12;
	RScalarREG<WordREG<REALT> > lrtmp12;
	    
	lr_zero12 = 0;
	    
	for(int i12 = 0; i12 < Nc; ++i12) {
      
	  lrtmp12 = tri_dia_r.elem(0).elem(i12);
	  lrtmp12 -= tri_dia_r.elem(0).elem(i12+Nc);
	  lrtmp12 -= tri_dia_r.elem(1).elem(i12);
	  lrtmp12 += tri_dia_r.elem(1).elem(i12+Nc);
	  B_j.elem().elem(i12,i12) = cmplx(lr_zero12,lrtmp12);
	}
	    
	/*# From lower triangular portion */
	int elem_ij12 = 0;
	for(int i12 = 1; i12 < Nc; ++i12) {
	      
	  int elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc;
	      
	  for(int j12 = 0; j12 < i12; ++j12) {
	
	    lctmp12 = tri_off_r.elem(0).elem(elem_ij12);
	    lctmp12 -= tri_off_r.elem(0).elem(elem_ijb12);
	    lctmp12 -= tri_off_r.elem(1).elem(elem_ij12);
	    lctmp12 += tri_off_r.elem(1).elem(elem_ijb12);
		
	    B_j.elem().elem(i12,j12) = timesI(lctmp12);
	    B_j.elem().elem(j12,i12) = timesI(adj(lctmp12));
		
	    elem_ij12++;
	    elem_ijb12++;
	  }
	}
      }
      sw.case_end();

      sw.case_default();
      {
      }
      sw.case_end();
    }

    jit_get_function(function);
  }



   
  template<typename T, typename U>
  void JITCloverTermT<T,U>::triacntr(U& B, int mat, int cb) const
  {
    START_CODE();

    B = zero;

    if ( mat < 0  ||  mat > 15 )
      {
	QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << std::endl;
	QDP_abort(1);
      }

    //QDPIO::cout << "PTX Clover triacntr " << (void*)this << "\n";
    //std::cout << "PTX Clover triacntr " << (void*)this << "\n";
    static JitFunction function;

    if (function.empty())
      function_triacntr_build<U>( function, B, tri_dia, tri_off, mat, rb[cb] );

    // Execute the function
    function_triacntr_exec(function, B, tri_dia, tri_off, mat, rb[cb] );

    END_CODE();
  }

  //! Returns the appropriate clover coefficient for indices mu and nu
  template<typename T, typename U>
  Real
  JITCloverTermT<T,U>::getCloverCoeff(int mu, int nu) const 
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



  template<typename T,typename X,typename Y>
  void function_apply_clov_exec(JitFunction& function,
				T& chi,
				const T& psi,
				const X& tri_dia,
				const Y& tri_off,
				const Subset& s)
  {
#ifdef QDP_DEEP_LOG
    function.type_W = typeid(REAL).name();
    function.set_dest_id( chi.getId() );
    function.set_is_lat(true);
#endif
    
    AddressLeaf addr_leaf(s);

    forEach(chi, addr_leaf, NullCombine());
    forEach(psi, addr_leaf, NullCombine());
    forEach(tri_dia, addr_leaf, NullCombine());
    forEach(tri_off, addr_leaf, NullCombine());

    int th_count = s.numSiteTable();
    WorkgroupGuardExec workgroupGuardExec(th_count);

    std::vector<QDPCache::ArgKey> ids;
    workgroupGuardExec.check(ids);
    ids.push_back( s.getIdSiteTable() );
    for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
      ids.push_back( addr_leaf.ids[i] );
    jit_launch(function,th_count,ids);
  }




  template<typename T,typename X,typename Y>
  void function_apply_clov_build( JitFunction& function,
				  const T& chi,
				  const T& psi,
				  const X& tri_dia,
				  const Y& tri_off,
				  const Subset& s)
  {
    llvm_start_new_function("apply_clov",__PRETTY_FUNCTION__);

    WorkgroupGuard workgroupGuard;
    ParamRef p_site_table = llvm_add_param<int*>();

    ParamLeafScalar param_leaf;

    typedef typename LeafFunctor<T, ParamLeafScalar>::Type_t  TJIT;
    TJIT chi_jit(forEach(chi, param_leaf, TreeCombine()));
    TJIT psi_jit(forEach(psi, param_leaf, TreeCombine()));
    typename REGType< typename ScalarType<typename TJIT::Subtype_t>::Type_t >::Type_t psi_r;
    typename REGType< typename ScalarType<typename TJIT::Subtype_t>::Type_t >::Type_t chi_r;

    typedef typename LeafFunctor<X, ParamLeafScalar>::Type_t  XJIT;
    XJIT tri_dia_jit(forEach(tri_dia, param_leaf, TreeCombine()));
    typename REGType< typename XJIT::Subtype_t >::Type_t tri_dia_r;

    typedef typename LeafFunctor<Y, ParamLeafScalar>::Type_t  YJIT;
    YJIT tri_off_jit(forEach(tri_off, param_leaf, TreeCombine()));
    typename REGType< typename YJIT::Subtype_t >::Type_t tri_off_r;

    llvm::Value* r_idx_thread = llvm_thread_idx();

    workgroupGuard.check(r_idx_thread);

    llvm::Value* r_idx = llvm_array_type_indirection<int>( p_site_table , r_idx_thread );

    auto chi_j = chi_jit.elem(JitDeviceLayout::Coalesced,r_idx);
    psi_r.setup( psi_jit.elem(JitDeviceLayout::Coalesced,r_idx) );
    tri_dia_r.setup( tri_dia_jit.elem(JitDeviceLayout::Coalesced,r_idx) );
    tri_off_r.setup( tri_off_jit.elem(JitDeviceLayout::Coalesced,r_idx) );

    // RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    // const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));

    int n = 2*Nc;

    for(int i = 0; i < n; ++i)
      {
	chi_r.elem((0*n+i)/3).elem((0*n+i)%3) = tri_dia_r.elem(0).elem(i) * psi_r.elem((0*n+i)/3).elem((0*n+i)%3);
	// cchi[0*n+i] = tri[site].diag[0][i] * ppsi[0*n+i];

	chi_r.elem((1*n+i)/3).elem((1*n+i)%3) = tri_dia_r.elem(1).elem(i) * psi_r.elem((1*n+i)/3).elem((1*n+i)%3);
	// cchi[1*n+i] = tri[site].diag[1][i] * ppsi[1*n+i];
      }

    int kij = 0;  
    for(int i = 0; i < n; ++i)
      {
	for(int j = 0; j < i; j++)
	  {
	    chi_r.elem((0*n+i)/3).elem((0*n+i)%3) += tri_off_r.elem(0).elem(kij) * psi_r.elem((0*n+j)/3).elem((0*n+j)%3);
	    // cchi[0*n+i] += tri[site].offd[0][kij] * ppsi[0*n+j];

	    chi_r.elem((0*n+j)/3).elem((0*n+j)%3) += conj(tri_off_r.elem(0).elem(kij)) * psi_r.elem((0*n+i)/3).elem((0*n+i)%3);
	    // cchi[0*n+j] += conj(tri[site].offd[0][kij]) * ppsi[0*n+i];

	    chi_r.elem((1*n+i)/3).elem((1*n+i)%3) += tri_off_r.elem(1).elem(kij) * psi_r.elem((1*n+j)/3).elem((1*n+j)%3);
	    // cchi[1*n+i] += tri[site].offd[1][kij] * ppsi[1*n+j];

	    chi_r.elem((1*n+j)/3).elem((1*n+j)%3) += conj(tri_off_r.elem(1).elem(kij)) * psi_r.elem((1*n+i)/3).elem((1*n+i)%3);
	    // cchi[1*n+j] += conj(tri[site].offd[1][kij]) * ppsi[1*n+i];

	    kij++;
	  }
      }

    chi_j = chi_r;

    jit_get_function(function);
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
  template<typename T, typename U>
  void JITCloverTermT<T,U>::apply(T& chi, const T& psi, 
				  enum PlusMinus isign, int cb) const
  {
    START_CODE();
    
    if ( Ns != 4 ) {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    //QDPIO::cout << "PTX Clover apply"  << (void*)this << "\n";
    //std::cout << "PTX Clover apply"  << (void*)this << "\n";
    static JitFunction function;

    if (function.empty())
      function_apply_clov_build( function, chi, psi, tri_dia, tri_off, rb[cb] );

    // Execute the function
    function_apply_clov_exec(function, chi, psi, tri_dia, tri_off, rb[cb] );

    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }



#ifndef  BUILD_QUDA_DEVIFACE_CLOVER
  namespace QDPCloverEnv {
    template<typename R,typename TD,typename TO> 
    struct QUDAPackArgs { 
      int cb;
      multi1d<QUDAPackedClovSite<R> >& quda_array;
      const TD&  tri_dia;
      const TO&  tri_off;
    };
    
    template<typename R,typename TD,typename TO>
    void qudaPackSiteLoop(int lo, int hi, int myId, QUDAPackArgs<R,TD,TO>* a) {
      int cb = a->cb;
      int Ns2 = Ns/2;

      multi1d<QUDAPackedClovSite<R> >& quda_array = a->quda_array;

      const TD& tri_dia = a->tri_dia;
      const TO& tri_off = a->tri_off;

      const int idtab[15]={0,1,3,6,10,2,4,7,11,5,8,12,9,13,14};

      for(int ssite=lo; ssite < hi; ++ssite) {
	int site = rb[cb].siteTable()[ssite];
	// First Chiral Block
	for(int i=0; i < 6; i++) { 
	  quda_array[site].diag1[i] = tri_dia.elem(site).comp[0].diag[i].elem().elem();
	}

	int target_index=0;
	
	for(int col=0; col < Nc*Ns2-1; col++) { 
	  for(int row=col+1; row < Nc*Ns2; row++) {

	    int source_index = row*(row-1)/2 + col;

	    quda_array[site].offDiag1[target_index][0] = tri_off.elem(site).comp[0].offd[source_index].real().elem();
	    quda_array[site].offDiag1[target_index][1] = tri_off.elem(site).comp[0].offd[source_index].imag().elem();
	    target_index++;
	  }
	}
	// Second Chiral Block
	for(int i=0; i < 6; i++) { 
	  quda_array[site].diag2[i] = tri_dia.elem(site).comp[1].diag[i].elem().elem();
	}

	target_index=0;
	for(int col=0; col < Nc*Ns2-1; col++) { 
	  for(int row=col+1; row < Nc*Ns2; row++) {

	    int source_index = row*(row-1)/2 + col;

	    quda_array[site].offDiag2[target_index][0] = tri_off.elem(site).comp[1].offd[source_index].real().elem();
	    quda_array[site].offDiag2[target_index][1] = tri_off.elem(site).comp[1].offd[source_index].imag().elem();
	    target_index++;
	  }
	}
      }
      QDPIO::cout << "\n";
    }
  }

  template<typename T, typename U>
  void JITCloverTermT<T,U>::packForQUDA(multi1d<QUDAPackedClovSite<typename WordType<T>::Type_t> >& quda_array, int cb) const
    {
      typedef typename WordType<T>::Type_t REALT;
      int num_sites = rb[cb].siteTable().size();

      typedef OLattice<PComp<PTriDia<RScalar <Word<REALT> > > > > TD;
      typedef OLattice<PComp<PTriOff<RComplex<Word<REALT> > > > > TO;

      StopWatch watch;
      watch.start();

      QDPCloverEnv::QUDAPackArgs<REALT,TD,TO> args = { cb, quda_array , tri_dia , tri_off };
      dispatch_to_threads(num_sites, args, QDPCloverEnv::qudaPackSiteLoop<REALT,TD,TO>);

      watch.stop();
      PackForQUDATimer::Instance().get() += watch.getTimeInMicroseconds();
    }

#endif

  template<typename T, typename U>
  void JITCloverTermT<T,U>::applySite(T& chi, const T& psi, 
				      enum PlusMinus isign, int site) const
  {
    QDP_error_exit("JITCloverTermT<T,U>::applySite(T& chi, const T& psi,..) not implemented ");
  }

  typedef JITCloverTermT<LatticeFermion, LatticeColorMatrix> JITCloverTerm;
  typedef JITCloverTermT<LatticeFermionF, LatticeColorMatrixF> JITCloverTermF;
  typedef JITCloverTermT<LatticeFermionD, LatticeColorMatrixD> JITCloverTermD;
} // End Namespace Chroma




#endif
