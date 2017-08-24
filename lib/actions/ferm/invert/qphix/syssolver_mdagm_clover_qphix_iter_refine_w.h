// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_clover_qphix_iter_refine_h__
#define __syssolver_mdagm_clover_qphix_iter_refine_h__

#include "chroma_config.h"

#ifdef BUILD_QPHIX
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "update/molecdyn/predictor/null_predictor.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>
#include "util/gauge/reunit.h"
#include "io/aniso_io.h"

// Header files from the dslash package
#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/clover.h"
#include "qphix/invbicgstab.h"
#include "qphix/invcg.h"
#include "qphix/inv_richardson_multiprec.h"

#include "actions/ferm/invert/qphix/qphix_vec_traits.h"
#include "qphix_singleton.h"

#include <memory>

using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverQPhiXCloverIterRefineEnv
 {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Clover Fermion System using the QPhiX inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
  using namespace QPhiXVecTraits;

  template<typename T, typename U>  
  class MdagMSysSolverQPhiXCloverIterRefine : public MdagMSystemSolver<T>
  {
  public:
    using Q = multi1d<U>;
    using REALT =  typename WordType<T>::Type_t;

    using InnerReal =  CHROMA_QPHIX_INNER_TYPE;

    static const int OuterVec = MixedVecTraits<REALT,InnerReal>::Vec;
    static const int OuterSoa = MixedVecTraits<REALT,InnerReal>::Soa;
    static const int InnerVec = MixedVecTraits<REALT,InnerReal>::VecInner;
    static const int InnerSoa = MixedVecTraits<REALT,InnerReal>::SoaInner;
    static const bool comp12 = MixedVecTraits<REALT,InnerReal>::compress12;

    using OuterGeom = QPhiX::Geometry<REALT,OuterVec,OuterSoa,comp12>;
    using InnerGeom = QPhiX::Geometry<InnerReal,InnerVec,InnerSoa,comp12>;

    using QPhiX_Spinor =  typename OuterGeom::FourSpinorBlock;
    using QPhiX_Gauge  =  typename OuterGeom::SU3MatrixBlock;
    using QPhiX_Clover =  typename OuterGeom::CloverBlock;

    using QPhiX_InnerSpinor = typename InnerGeom::FourSpinorBlock;
    using QPhiX_InnerGauge  = typename InnerGeom::SU3MatrixBlock;
    using QPhiX_InnerClover = typename InnerGeom::CloverBlock;

    
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverQPhiXCloverIterRefine(Handle< LinearOperator<T> > A_,
			      Handle< FermState<T,Q,Q> > state_,
			      const SysSolverQPhiXCloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new CloverTermT<T, U>()), invclov(new CloverTermT<T, U>())
    {

      
      QDPIO::cout << "MdagMSysSolverQPhiXCloverIterRefine:" << std::endl;

      if ( toBool( invParam.Delta < 0 ) ) { 
	QDPIO::cout << "Error: Delta parameter for solve should be set" << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << std::endl;


      QDPIO::cout << "Veclen is " << OuterVec << std::endl;
      QDPIO::cout << "Soalen is " << OuterSoa << std::endl;
      QDPIO::cout << "Inner Veclen is " << InnerVec << std::endl;
      QDPIO::cout << "Inner Soalen is " << InnerSoa << std::endl;
      
      if ( OuterSoa > OuterVec ) {
	QDPIO::cerr << "PROBLEM: Soalen > Veclen. Please set soalen appropriately (<=VECLEN) at compile time" << std::endl;
	QDP_abort(1);
      }

      if ( InnerSoa > InnerVec ) {
	QDPIO::cerr << "PROBLEM: Inner Soalen > Inner Veclen. Please set soalen appropriately at compile time" << std::endl;
	QDP_abort(1);
      }

      Q u(Nd);
      for(int mu=0; mu < Nd; mu++) {
    	  u[mu] = state_->getLinks()[mu];
      }

      // Set up aniso coefficients
      multi1d<Real> aniso_coeffs(Nd);
      for(int mu=0; mu < Nd; mu++) aniso_coeffs[mu] = Real(1);

      bool anisotropy = invParam.CloverParams.anisoParam.anisoP;
      if( anisotropy ) { 
    	  aniso_coeffs = makeFermCoeffs( invParam.CloverParams.anisoParam );
      }
      
      double t_boundary=(double)(1);
      // NB: In this case, the state will have boundaries applied.
      // So we only need to apply our own boundaries if Compression is enabled
      //
      if (invParam.AntiPeriodicT) {
    	  t_boundary=(double)(-1);
    	  // Flip off the boundaries -- Dslash expects them off..
    	  u[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
    			  Real(t_boundary), Real(1));
      }

      cbsize_in_blocks = rb[0].numSiteTable()/OuterSoa;

      const QPhiX::QPhiXCLIArgs& QPhiXParams = TheQPhiXParams::Instance();
#ifdef QDP_IS_QDPJIT
      int pad_xy = 0;
      int pad_xyz = 0;
#else
      int pad_xy = QPhiXParams.getPxy();
      int pad_xyz = QPhiXParams.getPxyz();
#endif
      n_blas_simt = QPhiXParams.getSy()*QPhiXParams.getSz();

      // Grab a dslash from which we will get geometry.
      geom_outer = new OuterGeom(	Layout::subgridLattSize().slice(),
									QPhiXParams.getBy(),
									QPhiXParams.getBz(),
									QPhiXParams.getNCores(),
									QPhiXParams.getSy(),
									QPhiXParams.getSz(),
									pad_xy,
									pad_xyz,
									QPhiXParams.getMinCt());

      
      geom_inner = new InnerGeom( 	Layout::subgridLattSize().slice(),
    		  	  	  	  	  	  	QPhiXParams.getBy(),
									QPhiXParams.getBz(),
									QPhiXParams.getNCores(),
									QPhiXParams.getSy(),
									QPhiXParams.getSz(),
									pad_xy,
									pad_xyz,
									QPhiXParams.getMinCt());
						     

#ifndef QDP_IS_QDPJIT
      psi_qphix=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      chi_qphix=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      tmp_qphix=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
#else
       psi_qphix=nullptr;
       chi_qphix=nullptr;
       tmp_qphix=nullptr;

#endif

      // Pack the gauge field
      QDPIO::cout << "Packing gauge field..." ;
      // Alloc outer gauge field
      QPhiX_Gauge* packed_gauge_cb0=(QPhiX_Gauge *)geom_outer->allocCBGauge();
      QPhiX_Gauge* packed_gauge_cb1=(QPhiX_Gauge *)geom_outer->allocCBGauge();

      // Alloc inner gauge field
      QPhiX_InnerGauge* packed_gauge_cb0_inner=(QPhiX_InnerGauge *)geom_inner->allocCBGauge();
      QPhiX_InnerGauge* packed_gauge_cb1_inner=(QPhiX_InnerGauge *)geom_inner->allocCBGauge();

      // Pack outer gauge field
      QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, *geom_outer);
      u_packed[0] = packed_gauge_cb0;
      u_packed[1] = packed_gauge_cb1;
      
      // Pack inner gauge field
      QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0_inner,packed_gauge_cb1_inner, *geom_inner);
      u_packed_i[0] = packed_gauge_cb0_inner;
      u_packed_i[1] = packed_gauge_cb1_inner;

      
      QDPIO::cout << "Creating Clover Term" << std::endl;
      CloverTerm clov_qdp;
      clov->create(state_, invParam.CloverParams);
      QDPIO::cout << "Inverting Clover Term" << std::endl;
      invclov->create(state_, invParam.CloverParams, (*clov));
      for(int cb=0; cb < 2; cb++) { 
	invclov->choles(cb);
      }
      QDPIO::cout << "Done" << std::endl;
      // Create buffer for outer clover
      QPhiX_Clover* A_cb0=(QPhiX_Clover *)geom_outer->allocCBClov();
      QPhiX_Clover* A_cb1=(QPhiX_Clover *)geom_outer->allocCBClov();
      clov_packed[0] = A_cb0;
      clov_packed[1] = A_cb1;

      QPhiX_Clover* A_inv_cb0=(QPhiX_Clover *)geom_outer->allocCBClov();
      QPhiX_Clover* A_inv_cb1=(QPhiX_Clover *)geom_outer->allocCBClov();
      invclov_packed[0] = A_inv_cb0;
      invclov_packed[1] = A_inv_cb1;

      // Create buffer for inner clover
      QPhiX_InnerClover* A_cb0_inner=(QPhiX_InnerClover *)geom_inner->allocCBClov();
      QPhiX_InnerClover* A_cb1_inner=(QPhiX_InnerClover *)geom_inner->allocCBClov();
      clov_packed_i[0] = A_cb0_inner;
      clov_packed_i[1] = A_cb1_inner;

      QPhiX_InnerClover* A_inv_cb0_inner=(QPhiX_InnerClover *)geom_inner->allocCBClov();
      QPhiX_InnerClover* A_inv_cb1_inner=(QPhiX_InnerClover *)geom_inner->allocCBClov();
      invclov_packed_i[0] = A_inv_cb0_inner;
      invclov_packed_i[1] = A_inv_cb1_inner;

      // Do the actual packing
      QDPIO::cout << "Packing Clover term..." << std::endl;
      
      // Pack outer clover inverse
      for(int cb=0; cb < 2; cb++) { 
    	  QPhiX::qdp_pack_clover<>((*invclov), invclov_packed[cb], *geom_outer, cb);
      }
    
      // Pack outer clover
      for(int cb=0; cb < 2; cb++) { 
    	  QPhiX::qdp_pack_clover<>((*clov), clov_packed[cb], *geom_outer, cb);
      }

      // Pack inner clover inverse
      for(int cb=0; cb < 2; cb++) { 
    	  QPhiX::qdp_pack_clover<>((*invclov), invclov_packed_i[cb], *geom_inner, cb);
      }
    
      // Pack inner clover
      for(int cb=0; cb < 2; cb++) { 
    	  QPhiX::qdp_pack_clover<>((*clov), clov_packed_i[cb], *geom_inner, cb);
      }
      QDPIO::cout << "Done" << std::endl;


      M_outer=new QPhiX::EvenOddCloverOperator<REALT,OuterVec,OuterSoa,comp12>(u_packed,
						     clov_packed[1], 
						     invclov_packed[0],  
						     geom_outer,
						     t_boundary,
						     toDouble(aniso_coeffs[0]),
						     toDouble(aniso_coeffs[3]));

      M_inner = new QPhiX::EvenOddCloverOperator<InnerReal, InnerVec, InnerSoa, comp12>(
          u_packed_i,
          clov_packed_i[1],
          invclov_packed_i[0],
          geom_inner,
          t_boundary,
          toDouble(aniso_coeffs[0]),
          toDouble(aniso_coeffs[3]));

      switch (invParam.SolverType) {
      case CG: {
        QDPIO::cout << "Creating the CG Solver" << std::endl;
        inner_solver = new QPhiX::InvCG<InnerReal, InnerVec, InnerSoa, comp12>(
            *M_inner, invParam.MaxIter);

        mixed_solver.reset(new QPhiX::InvRichardsonMultiPrec<REALT,
                                                             OuterVec,
                                                             OuterSoa,
                                                             comp12,
                                                             InnerReal,
                                                             InnerVec,
                                                             InnerSoa,
                                                             comp12,
                                                             true>(
            *M_outer, *inner_solver, toDouble(invParam.Delta), invParam.MaxIter));

        break;
      }
      case BICGSTAB: {
        QDPIO::cout << "Creating the BiCGStab Solver" << std::endl;
        inner_solver = new QPhiX::InvBiCGStab<InnerReal, InnerVec, InnerSoa, comp12>(
            *M_inner, invParam.MaxIter);

        mixed_solver.reset(new QPhiX::InvRichardsonMultiPrec<REALT,
                                                             OuterVec,
                                                             OuterSoa,
                                                             comp12,
                                                             InnerReal,
                                                             InnerVec,
                                                             InnerSoa,
                                                             comp12>(
            *M_outer, *inner_solver, toDouble(invParam.Delta), invParam.MaxIter));
        break;
      }
      default: {
        QDPIO::cerr << "UNKNOWN Solver Type" << std::endl;
        QDP_abort(1);
        break;
      }
      }
    }

  //! Destructor
  ~MdagMSysSolverQPhiXCloverIterRefine()
  {

    // Need to unalloc all the memory...
    QDPIO::cout << "Destructing" << std::endl;

#ifndef QDP_IS_QDPJIT
      geom_outer->free(psi_qphix);
      geom_outer->free(chi_qphix);
      geom_outer->free(tmp_qphix);
#else
      psi_qphix = nullptr;
      chi_qphix = nullptr;
      tmp_qphix = nullptr;
#endif
      geom_outer->free(invclov_packed[0]);
      geom_outer->free(invclov_packed[1]);
      invclov_packed[0] = nullptr;
      invclov_packed[1] = nullptr;


      geom_outer->free(clov_packed[0]);
      geom_outer->free(clov_packed[1]);
      clov_packed[0] = nullptr;
      clov_packed[1] = nullptr;


      geom_outer->free(u_packed[0]);
      geom_outer->free(u_packed[1]);
      u_packed[0] = nullptr;
      u_packed[1] = nullptr;


      geom_inner->free(invclov_packed_i[0]);
      geom_inner->free(invclov_packed_i[1]);
      invclov_packed_i[0] = nullptr;
      invclov_packed_i[1] = nullptr;


      geom_inner->free(clov_packed_i[0]);
      geom_inner->free(clov_packed_i[1]);
      clov_packed_i[0] = nullptr;
      clov_packed_i[1] = nullptr;


      geom_inner->free(u_packed_i[0]);
      geom_inner->free(u_packed_i[1]);
      u_packed_i[0] = nullptr;
      u_packed_i[1] = nullptr;

      delete geom_inner;      
      delete geom_outer;


    }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    SystemSolverResults_t operator()(
        T &psi, const T &chi, Chroma::AbsChronologicalPredictor4D<T> &predictor) const
    {

      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      SystemSolverResults_t res;
      switch (invParam.SolverType) {
      case CG:
        res = cgSolve(psi, chi, predictor);
        break;
      case BICGSTAB:
        res = biCGStabSolve(psi, chi, predictor);
        break;
      default:
        QDPIO::cout << "Unknown Solver " << std::endl;
        break;
      }

      swatch.stop();
      QDPIO::cout << "QPHIX_MDAGM_ITER_REFINE_SOLVER: total time: "
                  << swatch.getTimeInSeconds() << " (sec)" << endl;
      END_CODE();
      return res;
    }

    SystemSolverResults_t operator()(T& psi, const T& chi) const
    {
    	START_CODE();
    	SystemSolverResults_t res;
    	Null4DChronoPredictor not_predicting;
    	res=(*this)(psi,chi,not_predicting);
    	END_CODE();
    	return res;

    }

    SystemSolverResults_t
    cgSolve(T &psi, const T &chi, AbsChronologicalPredictor4D<T> &predictor) const
    {
      SystemSolverResults_t res;
      Handle<LinearOperator<T>> MdagM(new MdagMLinOp<T>(A));
      predictor(psi, (*MdagM), chi);

#ifndef QDP_IS_QDPJIT
      // Pack Spinors psi and chi
      QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom_outer, 1);
      QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom_outer, 1);
#else
      psi_qphix = (QPhiX_Spinor *)(psi.getFjit()) + cbsize_in_blocks;
      chi_qphix = (QPhiX_Spinor *)(chi.getFjit()) + cbsize_in_blocks;
#endif
      double rsd_final;
      unsigned long site_flops = 0;
      unsigned long mv_apps = 0;

      double start = omp_get_wtime();
      int my_isign = 1;
      (*mixed_solver)(psi_qphix,
                      chi_qphix,
                      toDouble(invParam.RsdTarget),
                      res.n_count,
                      rsd_final,
                      site_flops,
                      mv_apps,
                      my_isign,
                      invParam.VerboseP);
      double end = omp_get_wtime();

#ifndef QDP_IS_QDPJIT
      QPhiX::qdp_unpack_cb_spinor<>(psi_qphix, psi, *geom_outer, 1);
#endif

      predictor.newVector(psi);

      // Chi Should now hold the result spinor
      // Check it against chroma.
      {
        T r = chi;
        T tmp, tmp2;
        (*A)(tmp, psi, PLUS);
        (*A)(tmp2, tmp, MINUS);
        r[A->subset()] -= tmp2;

        Double r2 = norm2(r, A->subset());
        Double b2 = norm2(chi, A->subset());
        Double rel_resid = sqrt(r2 / b2);
        res.resid = rel_resid;
        QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: " << res.n_count
                    << " iters,  rsd_sq_final=" << rel_resid << std::endl;

        QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: || r || / || b || = " << rel_resid
                    << std::endl;

#if 0
	if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	  QDPIO::cout << "SOLVE FAILED" << std::endl;
	  QDP_abort(1);
	}
#endif
      }

      int num_cb_sites = Layout::vol() / 2;
      unsigned long total_flops =
          (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;
      double gflops = (double)(total_flops) / (1.0e9);

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: Solver Time=" << total_time
                  << " (sec)  Performance=" << gflops / total_time << " GFLOPS"
                  << std::endl;

      END_CODE();
      return res;
    }

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t
    biCGStabSolve(T &psi, const T &chi, AbsChronologicalPredictor4D<T> &predictor) const
    {
    	SystemSolverResults_t res;
    	Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );

    	T Y;
    	Y[ A->subset() ] = psi; // Y is initial guess

    	try  {
    		// Try to cast the predictor to a two step predictor
    		AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor =
    				dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);


    		// Predict Y and X separately
    		two_step_predictor.predictY(Y,*A,chi);
    		two_step_predictor.predictX(psi,*MdagM, chi);
    	}
    	catch( std::bad_cast) {

    		// Not a 2 step predictor. Predict X
    		// Then MX = Y is a good guess.
    		predictor(psi,*MdagM, chi);
    		(*A)(Y,psi,PLUS);

    	}


#ifndef QDP_IS_QDPJIT
    	QDPIO::cout << "Packing" << std::endl << std::flush ;
    	QPhiX::qdp_pack_cb_spinor<>(Y, tmp_qphix, *geom_outer,1); // Initial Guess for Y
    	QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom_outer,1); // Initial Guess for X
    	QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom_outer,1); // RHS
#else
    	tmp_qphix = (QPhiX_Spinor *)(Y.getFjit()) + cbsize_in_blocks;
    	psi_qphix = (QPhiX_Spinor *)(psi.getFjit()) + cbsize_in_blocks;
    	chi_qphix = (QPhiX_Spinor *)(chi.getFjit()) + cbsize_in_blocks;
#endif
    	QDPIO::cout << "Done" << std::endl << std::flush;
    	double rsd_final;
    	int num_cb_sites = Layout::vol()/2;

    	unsigned long site_flops1=0;
    	unsigned long mv_apps1=0;
    	int n_count1=0;

    	unsigned long site_flops2=0;
    	unsigned long mv_apps2=0;
    	int n_count2=0;

    	QDPIO::cout << "Starting Y solve" << std::endl << std::flush ;
    	double start = omp_get_wtime();
    	(*mixed_solver)(tmp_qphix,chi_qphix, toDouble(invParam.RsdTarget), n_count1, rsd_final, site_flops1, mv_apps1, -1, invParam.VerboseP);
    	double end = omp_get_wtime();


    	unsigned long total_flops = (site_flops1 + (1320+504+1320+504+48)*mv_apps1)*num_cb_sites;
    	double gflops = (double)(total_flops)/(1.0e9);
    	double total_time = end - start;

    	Double r_final = sqrt(toDouble(rsd_final));
    	QDPIO::cout << "QPHIX_CLOVER_ITER_REFINE_BICGSTAB_SOLVER: " << n_count1 << " iters,  rsd_sq_final=" << rsd_final << " ||r||/||b|| (acc) = " << r_final <<std::endl;
    	QDPIO::cout << "QPHIX_CLOVER_ITER_REFINE_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

    	QDPIO::cout << "Starting X solve" << std::endl << std::flush ;
    	start = omp_get_wtime();
    	(*mixed_solver)(psi_qphix,tmp_qphix, toDouble(invParam.RsdTarget), n_count2, rsd_final, site_flops2, mv_apps2, +1, invParam.VerboseP);
    	end = omp_get_wtime();
    	total_flops = (site_flops2 + (1320+504+1320+504+48)*mv_apps2)*num_cb_sites;
    	gflops = (double)(total_flops)/(1.0e9);
    	total_time = end - start;

#ifndef QDP_IS_QDPJIT
    	// Want Norm of Y
    	QPhiX::qdp_unpack_cb_spinor<>(tmp_qphix, Y, *geom_outer,1);
#endif
    	double norm2Y;
    	QPhiX::norm2Spinor(norm2Y, tmp_qphix, *geom_outer, n_blas_simt);
    	r_final = sqrt(toDouble(rsd_final));

    	QDPIO::cout << "QPHIX_CLOVER_ITER_REFINE_BICGSTAB_SOLVER: " << n_count2 << " iters,  rsd_sq_final=" << rsd_final << " ||r||/||b|| (acc) = " << r_final << std::endl;
    	QDPIO::cout << "QPHIX_CLOVER_ITER_REFINE_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

#ifndef QDP_IS_QDPJIT
    	QPhiX::qdp_unpack_cb_spinor<>(psi_qphix, psi, *geom_outer,1);
#endif

    	try  {
    		// Try to cast the predictor to a two step predictor
    		AbsTwoStepChronologicalPredictor4D<T>& two_step_predictor =
    				dynamic_cast<AbsTwoStepChronologicalPredictor4D<T>& >(predictor);
    		two_step_predictor.newYVector(Y);
    		two_step_predictor.newXVector(psi);

    	}
    	catch( std::bad_cast) {

    		// Not a 2 step predictor. Predict X
    		// Then MX = Y is a good guess.
    		predictor.newVector(psi);
    	}

    	// Chi Should now hold the result spinor
    	// Check it against chroma. -- reuse Y as the residuum
    	Y[ A->subset() ] = chi;
    	{
    		T tmp,tmp2;
    		(*A)(tmp, psi, PLUS);
    		(*A)(tmp2, tmp, MINUS);

    		Y[ A->subset() ] -= tmp2;
    	}

    	Double r2 = norm2(Y,A->subset());
    	Double b2 = norm2(chi, A->subset());
    	Double rel_resid = sqrt(r2/b2);
    	res.resid=rel_resid;
    	res.n_count = n_count1 + n_count2;
    	QDPIO::cout << "QPHIX_CLOVER_ITER_REFINE_BICGSTAB_SOLVER: total_iters="<<res.n_count<<" || r || / || b || = " << res.resid << std::endl;


#if 0
    	if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
    		QDPIO::cout << "SOLVE FAILED" << std::endl;
    		QDP_abort(1);
    	}
#endif
    	return res;
    }






  private:
    // Hide default constructor
    MdagMSysSolverQPhiXCloverIterRefine() {}
    

    // The LinOp
    Handle< LinearOperator<T> > A;

    // Params
    const SysSolverQPhiXCloverParams invParam;

    // Clover Terms
    Handle< CloverTermT<T, U> > clov;
    Handle< CloverTermT<T, U> > invclov;

    // Outer Geom
    QPhiX::Geometry<REALT, 
		    MixedVecTraits<REALT,InnerReal>::Vec, 
		    MixedVecTraits<REALT,InnerReal>::Soa, 
		    MixedVecTraits<REALT,InnerReal>::compress12>* geom_outer;

    // Inner Geom
    QPhiX::Geometry<InnerReal, 
		    MixedVecTraits<REALT,InnerReal>::VecInner, 
		    MixedVecTraits<REALT,InnerReal>::SoaInner, 
		    MixedVecTraits<REALT,InnerReal>::compress12>* geom_inner;

    // Outer M
    Handle< QPhiX::EvenOddCloverOperator<REALT, 
					 MixedVecTraits<REALT,InnerReal>::Vec, 
					 MixedVecTraits<REALT,InnerReal>::Soa, 
					 MixedVecTraits<REALT,InnerReal>::compress12> > M_outer;


    // Inner M
    Handle< QPhiX::EvenOddCloverOperator<InnerReal,
				     MixedVecTraits<REALT,InnerReal>::VecInner, 
				     MixedVecTraits<REALT,InnerReal>::SoaInner, 
				     MixedVecTraits<REALT,InnerReal>::compress12> > M_inner;
    
    // Inner solver, can be CG or BiCGStab
    Handle<QPhiX::AbstractSolver<InnerReal,
                                 MixedVecTraits<REALT, InnerReal>::VecInner,
                                 MixedVecTraits<REALT, InnerReal>::SoaInner,
                                 MixedVecTraits<REALT, InnerReal>::compress12>>
        inner_solver;

    // Outer solver, will be Richardson solver.
    std::unique_ptr<QPhiX::AbstractSolver<REALT, OuterVec, OuterSoa, comp12>>
        mixed_solver;

    QPhiX_Clover* invclov_packed[2];
    QPhiX_Clover* clov_packed[2];
    QPhiX_Gauge* u_packed[2];
    
    QPhiX_InnerClover* invclov_packed_i[2];
    QPhiX_InnerClover* clov_packed_i[2];
    QPhiX_InnerGauge* u_packed_i[2];
    

    mutable QPhiX_Spinor* psi_qphix;
    mutable QPhiX_Spinor* chi_qphix;
    mutable QPhiX_Spinor* tmp_qphix;
    size_t cbsize_in_blocks;
    int n_blas_simt;
    
  };


} // End namespace

#endif // BUILD_QPHIX
#endif 

