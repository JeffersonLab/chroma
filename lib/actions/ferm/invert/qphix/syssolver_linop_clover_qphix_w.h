// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_clover_qphix_h__
#define __syssolver_linop_clover_qphix_h__

#include "chroma_config.h"

#ifdef BUILD_QPHIX
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
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

#include "actions/ferm/invert/qphix/qphix_vec_traits.h"

using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQPhiXCloverEnv
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
  class LinOpSysSolverQPhiXClover : public LinOpSystemSolver<T>
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
    LinOpSysSolverQPhiXClover(Handle< LinearOperator<T> > A_,
			      Handle< FermState<T,Q,Q> > state_,
			      const SysSolverQPhiXCloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new CloverTermT<T, U>()), invclov(new CloverTermT<T, U>())
    {
      QDPIO::cout << "LinOpSysSolverQPhiXClover:" << std::endl;
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

      QDPIO::cout << "About to grap a Dslash" << std::endl;
      geom = new QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>(Layout::subgridLattSize().slice(),
												      invParam.By, 
												      invParam.Bz, 
												      invParam.NCores, 
												      invParam.Sy,
												      invParam.Sz,
												      invParam.PadXY,
												      invParam.PadXYZ,
														   invParam.MinCt);
      
      QDPIO::cout << " Allocating p and c" << std::endl << std::flush ;
      psi_qphix=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      chi_qphix=(QPhiX_Spinor *)geom->allocCBFourSpinor();

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
	  if( invParam.TuneP ) cg_solver->tune();

	}
	break;

      case BICGSTAB:
	{
	  QDPIO::cout << "Creating the BiCGStab Solver" << std::endl;
	  bicgstab_solver = new QPhiX::InvBiCGStab<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>((*M), invParam.MaxIter);
	  if( invParam.TuneP ) bicgstab_solver->tune();

	}
	break;
      default:
	QDPIO::cerr << "UNKNOWN Solver Type" << std::endl;
	QDP_abort(1);
      }
    }

    

    //! Destructor is automatic
    ~LinOpSysSolverQPhiXClover() 
    {
      
      // Need to unalloc all the memory...
      QDPIO::cout << "Destructing" << std::endl;

      geom->free(psi_qphix);
      geom->free(chi_qphix);
      psi_qphix = 0x0;
      chi_qphix = 0x0;
      
      geom->free(invclov_packed[0]);
      geom->free(invclov_packed[1]);
      invclov_packed[0] = 0x0;
      invclov_packed[1] = 0x0;
      
      geom->free(clov_packed[0]);
      geom->free(clov_packed[1]);
      clov_packed[0] = 0x0;
      clov_packed[1] = 0x0;
      
      geom->free(u_packed[0]);
      geom->free(u_packed[1]);
      u_packed[0] = 0x0;
      u_packed[1] = 0x0;

      delete geom;

      // Evertyhing else freed as handles freed.
    }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}
    

    
    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator()(T& psi, const T& chi) const
    {
      /* Factories here later? */
      SystemSolverResults_t res;
      switch( invParam.SolverType ) { 
      case CG:
	{
	  res = cgnrSolve(psi,chi);
	}
	break;
      case BICGSTAB:
	{
	  res = biCGStabSolve(psi,chi);
	}
	break;
      default:
	QDPIO::cout << "Unknown Solver " << std::endl;
	break;
      }
      
      return res;
    }





  private:
    // Hide default constructor
    LinOpSysSolverQPhiXClover() {}
    


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
    
    QPhiX_Spinor* psi_qphix;
    QPhiX_Spinor* chi_qphix;

    SystemSolverResults_t cgnrSolve(T& psi, const T& chi) const
    {

      SystemSolverResults_t res;
      psi = zero;

        // This is LinOpSolver so CGNE
      // Step 1 tranform RHS: chi -> M^\dagger chi
      // Move this into library later.
      T mdag_chi;
      (*A)(mdag_chi, chi, MINUS);
      
      
      //      QDPIO::cout << "Allocating Spinor fields" << std::endl;
      // Pack Spinors psi and chi
      QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom,1);
      QPhiX::qdp_pack_cb_spinor<>(mdag_chi, chi_qphix, *geom,1);
      
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      int my_isign=1;

      double start = omp_get_wtime();
      (*cg_solver)(psi_qphix,chi_qphix, toDouble(invParam.RsdTarget), res.n_count, rsd_final, site_flops, mv_apps, my_isign,invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << std::endl;      
      QPhiX::qdp_unpack_cb_spinor<>(psi_qphix, psi,*geom,1);

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: || r || / || b || = " << rel_resid << std::endl;


      if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	QDPIO::cout << "SOLVE FAILED" << std::endl;
	QDP_abort(1);
      }


      int num_cb_sites = Layout::vol()/2;
      unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
      double gflops = (double)(total_flops)/(1.0e9);

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

      END_CODE();
      return res;

    }

    SystemSolverResults_t biCGStabSolve(T& psi, const T& chi) const
    {
      SystemSolverResults_t res;
      psi=zero;

      // Pack Spinors psi and chi
      QDPIO::cout << "Packing" << std::endl << std::flush ;
      QPhiX::qdp_pack_cb_spinor<>(psi, psi_qphix, *geom,1);
      QPhiX::qdp_pack_cb_spinor<>(chi, chi_qphix, *geom,1);
      QDPIO::cout << "Done" << std::endl << std::flush;
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      QDPIO::cout << "Starting solve" << std::endl << std::flush ;
      double start = omp_get_wtime();
      (*bicgstab_solver)(psi_qphix,chi_qphix, toDouble(invParam.RsdTarget), res.n_count, rsd_final, site_flops, mv_apps, 1, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << std::endl;      
      //      QPhiX::qdp_unpack_spinor<>(psi_s[0], psi_qphix, psi, *geom);
      QPhiX::qdp_unpack_cb_spinor<>( psi_qphix, psi, *geom,1);

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: || r || / || b || = " << rel_resid << std::endl;

#if 0
      if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	QDPIO::cout << "SOLVE FAILED" << std::endl;
	QDP_abort(1);
      }
#endif

      int num_cb_sites = Layout::vol()/2;
      unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
      double gflops = (double)(total_flops)/(1.0e9);

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

      END_CODE();
      return res;
    }

  };


} // End namespace

#endif // BUILD_QPHIX
#endif 

