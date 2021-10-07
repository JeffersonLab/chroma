// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_mdagm_clover_qphix_h__
#define __syssolver_mdagm_clover_qphix_h__

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

#include "chroma_config.h"
#include "util/gauge/reunit.h"
#include "io/aniso_io.h"

// Header files from the dslash package
#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"

#include "qphix_singleton.h"
#include "actions/ferm/invert/qphix/qphix_vec_traits.h"
#include "qphix/blas_new_c.h"
using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace MdagMSysSolverQPhiXCloverEnv
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
  class MdagMSysSolverQPhiXClover : public MdagMSystemSolver<T>
  {
  public:
    typedef multi1d<U> Q;
    typedef typename WordType<T>::Type_t REALT;
    typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::FourSpinorBlock QPhiX_Spinor;
    typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::SU3MatrixBlock QPhiX_Gauge;
    typedef typename QPhiX::Geometry<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>::CloverBlock QPhiX_Clover;

    
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMSysSolverQPhiXClover(Handle< LinearOperator<T> > A_,
			      Handle< FermState<T,Q,Q> > state_,
			      const SysSolverQPhiXCloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new CloverTermT<T, U>()), invclov(new CloverTermT<T, U>())
    {
      QDPIO::cout << "MdagMSysSolverQPhiXClover:" << std::endl;
      QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << std::endl;
      

      QDPIO::cout << "Veclen is " << VecTraits<REALT>::Vec << std::endl;
      QDPIO::cout << "Soalen is " << VecTraits<REALT>::Soa << std::endl;
      if ( VecTraits<REALT>::Soa > VecTraits<REALT>::Vec ) { 
	QDPIO::cerr << "PROBLEM: Soalen > Veclen. Please set soalen appropriately (<=VECLEN) at compile time" << std::endl;
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
      
      cbsize_in_blocks = rb[0].numSiteTable()/VecTraits<REALT>::Soa;

      const QPhiX::QPhiXCLIArgs& QPhiXParams = TheQPhiXParams::Instance();

#ifdef QDP_IS_QDPJIT
      int pad_xy=0;
      int pad_xyz=0;
#else
      int pad_xy = QPhiXParams.getPxy();
      int pad_xyz= QPhiXParams.getPxyz();
#endif
      n_blas_simt = QPhiXParams.getSy()*QPhiXParams.getSz();

      QDPIO::cout << "About to grap a Dslash" << std::endl;
      geom = new QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>(Layout::subgridLattSize().slice(),
														   QPhiXParams.getBy(),
														   QPhiXParams.getBz(),
														   QPhiXParams.getNCores(),
														   QPhiXParams.getSy(),
														   QPhiXParams.getSz(),
														   pad_xy,
														   pad_xyz,
														   QPhiXParams.getMinCt());
      
      QDPIO::cout << " Allocating p and c" << std::endl << std::flush ;

#ifndef QPD_IS_QDPJIT
      psi_qphix=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      chi_qphix=(QPhiX_Spinor *)geom->allocCBFourSpinor();

      if( invParam.SolverType == BICGSTAB ) {
    	  // Two step solve need a temporary
    	  tmp_qphix=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      }
#else
      psi_qphix = nullptr;
      chi_qphix = nullptr;
      tmp_qphix = nullptr;

#endif
      
      
      QDPIO::cout << " Allocating Clover" << std::endl << std::flush ;
      QPhiX_Clover* A_cb0=(QPhiX_Clover *)geom->allocCBClov();
      QPhiX_Clover* A_cb1=(QPhiX_Clover *)geom->allocCBClov();
      clov_packed[0] = A_cb0;
      clov_packed[1] = A_cb1;
      
      QDPIO::cout << " Allocating CloverInv " << std::endl << std::flush ;
      QPhiX_Clover* A_inv_cb0=(QPhiX_Clover *)geom->allocCBClov();
      QPhiX_Clover* A_inv_cb1=(QPhiX_Clover *)geom->allocCBClov();
      invclov_packed[0] = A_inv_cb0;
      invclov_packed[1] = A_inv_cb1;
      
      // Pack the gauge field
      QDPIO::cout << "Packing gauge field..."  << std::endl << std::flush ;
      QPhiX_Gauge* packed_gauge_cb0=(QPhiX_Gauge *)geom->allocCBGauge();
      QPhiX_Gauge* packed_gauge_cb1=(QPhiX_Gauge *)geom->allocCBGauge();
      
      QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, *geom);
      u_packed[0] = packed_gauge_cb0;
      u_packed[1] = packed_gauge_cb1;
      
      
      QDPIO::cout << "Creating Clover Term" << std::endl;
      CloverTerm clov_qdp;
      clov->create(state_, invParam.CloverParams);
      QDPIO::cout << "Inverting Clover Term" << std::endl;
      invclov->create(state_, invParam.CloverParams, (*clov));
      for(int cb=0; cb < 2; cb++) { 
	invclov->choles(cb);
      }
      QDPIO::cout << "Done" << std::endl;
      QDPIO::cout << "Packing Clover term..." << std::endl;
      
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*invclov), invclov_packed[cb], *geom, cb);
      }
      
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*clov), clov_packed[cb], *geom, cb);
      }
      QDPIO::cout << "Done" << std::endl;
      
      QDPIO::cout << "Creating the Even Odd Operator" << std::endl;
      M=new QPhiX::EvenOddCloverOperator<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>(u_packed,  
															 clov_packed[1], 
															 invclov_packed[0],  
															 geom,
															 t_boundary,
															 toDouble(aniso_coeffs[0]),
															 toDouble(aniso_coeffs[3]));
      
      
      switch( invParam.SolverType ) { 
      case CG:
	{
	  QDPIO::cout << "Creating the CG Solver" << std::endl;
	  cg_solver = new QPhiX::InvCG<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>((*M), invParam.MaxIter);
	}
      case BICGSTAB:
	{
	  QDPIO::cout << "Creating the BiCGStab Solver" << std::endl;
	  bicgstab_solver = new QPhiX::InvBiCGStab<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>((*M), invParam.MaxIter);
	}
	break;
      default:
	QDPIO::cerr << "UNKNOWN Solver Type" << std::endl;
	QDP_abort(1);
      }
    }

    

    //! Destructor is automatic
    ~MdagMSysSolverQPhiXClover() 
    {
      
      // Need to unalloc all the memory...
      QDPIO::cout << "Destructing" << std::endl;

#ifndef QDP_IS_QDPJIT
      geom->free(psi_qphix);
      geom->free(chi_qphix);
      if( invParam.SolverType == BICGSTAB ) {
      	geom->free(tmp_qphix);
      }

#endif
      psi_qphix = nullptr;
      chi_qphix = nullptr;
      tmp_qphix = nullptr;
      
      geom->free(invclov_packed[0]);
      geom->free(invclov_packed[1]);
      invclov_packed[0] = nullptr;
      invclov_packed[1] = nullptr;
      
      geom->free(clov_packed[0]);
      geom->free(clov_packed[1]);
      clov_packed[0] = nullptr;
      clov_packed[1] = nullptr;
      
      geom->free(u_packed[0]);
      geom->free(u_packed[1]);
      u_packed[0] = nullptr;
      u_packed[1] = nullptr;

      delete geom;

      // Evertyhing else freed as handles freed.
    }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}
    

    SystemSolverResults_t operator()(T& psi, const T& chi, 
				     Chroma::AbsChronologicalPredictor4D<T>& predictor) const 
    {
      
      START_CODE();
      StopWatch swatch;
      swatch.reset(); swatch.start();
      /* Factories here later? */
      SystemSolverResults_t res;
      switch( invParam.SolverType ) { 
      case CG:
	{
	  res = cgSolve(psi,chi,predictor);
	}
	break;
      case BICGSTAB:
	{
	  res = biCGStabSolve(psi,chi,predictor);
	}
	break;
      default:
	QDPIO::cout << "Unknown Solver " << std::endl;
	break;
      }

      swatch.stop();
      QDPIO::cout << "QPHIX_MDAGM_SOLVER: total time: " << swatch.getTimeInSeconds() << " (sec)" << std::endl;
      END_CODE();
      return res;
    }

    
    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */

    // NEED THE CHRONO Predictor argument here... 
    SystemSolverResults_t operator()(T& psi, const T& chi) const
    {
       START_CODE();
       SystemSolverResults_t res;
       Null4DChronoPredictor not_predicting;
       res=(*this)(psi,chi, not_predicting);
     
       END_CODE();
       return res;
    }





  private:
    // Hide default constructor
    MdagMSysSolverQPhiXClover() {}
    


    Handle< LinearOperator<T> > A;
    const SysSolverQPhiXCloverParams invParam;
    Handle< CloverTermT<T, U> > clov;
    Handle< CloverTermT<T, U> > invclov;

    QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>* geom;
    
    Handle< QPhiX::EvenOddCloverOperator<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > M;

    Handle< QPhiX::InvCG<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > cg_solver;

    Handle< QPhiX::InvBiCGStab<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>  > bicgstab_solver;
    
    QPhiX_Clover* invclov_packed[2];
    QPhiX_Clover* clov_packed[2];
    QPhiX_Gauge* u_packed[2];
    
    mutable QPhiX_Spinor* psi_qphix;
    mutable QPhiX_Spinor* chi_qphix;
    mutable QPhiX_Spinor* tmp_qphix;

    size_t cbsize_in_blocks;
    int n_blas_simt;

  SystemSolverResults_t cgSolve(T& psi, const T& chi, AbsChronologicalPredictor4D<T>& predictor) const
  {
      SystemSolverResults_t res;
      Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
      predictor(psi, (*MdagM), chi);

#ifndef QDP_IS_QDPJIT
      // Pack Spinors psi and chi
      QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom,1);
      QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom,1);
#else
      psi_qphix = (QPhiX_Spinor*)(psi.getFjit()) + cbsize_in_blocks;
      chi_qphix = (QPhiX_Spinor*)(chi.getFjit()) + cbsize_in_blocks;
#endif
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      double start = omp_get_wtime();
      int my_isign=1;
      (*cg_solver)(psi_qphix, chi_qphix, toDouble(invParam.RsdTarget), res.n_count, rsd_final, site_flops, mv_apps, my_isign, invParam.VerboseP);
      double end = omp_get_wtime();

#ifndef QDP_IS_QDPJIT
      QPhiX::qdp_unpack_cb_spinor<>(psi_qphix, psi, *geom,1);
#endif

      predictor.newVector(psi);

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      {
	T r = chi;
	T tmp,tmp2;
	(*A)(tmp, psi, PLUS);
	(*A)(tmp2, tmp, MINUS);
	r[ A->subset() ] -= tmp2;
	
	Double r2 = norm2(r,A->subset());
	Double b2 = norm2(chi, A->subset());
	Double rel_resid = sqrt(r2/b2);
	res.resid = rel_resid;
      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rel_resid << std::endl;      

	QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: || r || / || b || = " << rel_resid << std::endl;
	
#if 0
	if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	  QDPIO::cout << "SOLVE FAILED" << std::endl;
	  QDP_abort(1);
	}
#endif
      }  

      size_t num_cb_sites = Layout::vol()/2;
      double total_flops = (static_cast<double>(site_flops) + static_cast<double>(1320+504+1320+504+48)
								*static_cast<double>(mv_apps))*static_cast<double>(num_cb_sites);
      double gflops = total_flops/1.0e9;

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

      END_CODE();
      return res;

    }


  SystemSolverResults_t biCGStabSolve(T& psi, const T& chi,AbsChronologicalPredictor4D<T>& predictor) const
    {
	  START_CODE();
	  StopWatch swatch;
	  swatch.reset(); swatch.start();

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
      QPhiX::qdp_pack_cb_spinor<>(Y, tmp_qphix, *geom,1); // Initial Guess for Y
      QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom,1); // Initial Guess for X
      QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom,1); // RHS
#else
      tmp_qphix = (QPhiX_Spinor *)(Y.getFjit()) + cbsize_in_blocks;
      psi_qphix = (QPhiX_Spinor *)(psi.getFjit()) + cbsize_in_blocks;
      chi_qphix = (QPhiX_Spinor *)(chi.getFjit()) + cbsize_in_blocks;
#endif
      QDPIO::cout << "Done" << std::endl << std::flush;
      double rsd_final;
      size_t num_cb_sites = Layout::vol()/2;

      unsigned long site_flops1=0;
      unsigned long mv_apps1=0;
      int n_count1=0;

      unsigned long site_flops2=0;
      unsigned long mv_apps2=0;
      int n_count2=0;

      QDPIO::cout << "Starting Y solve" << std::endl << std::flush ;
      double start = omp_get_wtime();
      (*bicgstab_solver)(tmp_qphix,chi_qphix, toDouble(invParam.RsdTarget), n_count1, rsd_final, site_flops1, mv_apps1, -1, invParam.VerboseP);
      double end = omp_get_wtime();


      double total_flops = (static_cast<double>(site_flops1) + static_cast<double>(1320+504+1320+504+48)*
							static_cast<double>(mv_apps1))*static_cast<double(num_cb_sites);
      double gflops = total_flops/1.0e9;
      double total_time = end - start;

      Double r_final = sqrt(toDouble(rsd_final)/norm2(chi,A->subset()));
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: " << n_count1 << " iters,  rsd_sq_final=" << rsd_final << " ||r||/||b|| (acc) = " << r_final <<std::endl;
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

      QDPIO::cout << "Starting X solve" << std::endl << std::flush ;
      start = omp_get_wtime();
      (*bicgstab_solver)(psi_qphix,tmp_qphix, toDouble(invParam.RsdTarget), n_count2, rsd_final, site_flops2, mv_apps2, +1, invParam.VerboseP);
      end = omp_get_wtime();
      total_flops = (static_cast<double>(site_flops2) + static_cast<double>(1320+504+1320+504+48)
										*static_cast<double>(mv_apps2))*static_cast<double>(num_cb_sites);
      gflops = total_flops/1.0e9;
      total_time = end - start;

#ifndef QDP_IS_QDPJIT
      // Want Norm of Y
      QPhiX::qdp_unpack_cb_spinor<>(tmp_qphix, Y, *geom,1);
#endif
      double norm2Y;
      QPhiX::norm2Spinor(norm2Y, tmp_qphix, *geom, n_blas_simt);
      r_final = sqrt(toDouble(rsd_final)/norm2Y);

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: " << n_count2 << " iters,  rsd_sq_final=" << rsd_final << " ||r||/||b|| (acc) = " << r_final << std::endl;
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

#ifndef QDP_IS_QDPJIT
      QPhiX::qdp_unpack_cb_spinor<>(psi_qphix, psi, *geom,1);
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
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: total_iters="<<res.n_count<<" || r || / || b || = " << res.resid << std::endl;


#if 0
      if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	QDPIO::cout << "SOLVE FAILED" << std::endl;
	QDP_abort(1);
      }
#endif
      
      swatch.stop();
      QDPIO::cout << "QPHIX_MDAGM_SOLVER: total time: " << swatch.getTimeInSeconds() << " (sec)" << std::endl;
      END_CODE();
      return res;
    }

  };


} // End namespace

#endif // BUILD_QPHIX
#endif 

