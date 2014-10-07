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
#include "actions/ferm/linop/clover_term_qdp_w.h"
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

using namespace std;
using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQPhiXCloverEnv
  {
    //! Register the syssolver
    bool registerAll();

    template<typename T>
    struct VecTraits { 
      static const int Vec=1;
      static const int Soa=1;
      static const bool compress12=false;
    };

    // Templates
#if defined CHROMA_QPHIX_ARCH_AVX
#warning QPHIX for AVX
    // AVX Traits:
    template<>
    struct VecTraits<float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    template<>
    struct VecTraits<double> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 

    };
#endif

#if defined CHROMA_QPHIX_ARCH_MIC
#warning QPhiX for MIC
    // MIC Traits
    template<>
    struct VecTraits<float> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };

    template<>
    struct VecTraits<double> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
    };
#endif


  }



  //! Solve a Clover Fermion System using the QPhiX inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
  using namespace LinOpSysSolverQPhiXCloverEnv;
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
      A(A_), invParam(invParam_), clov(new QDPCloverTermT<T, U>()), invclov(new QDPCloverTermT<T, U>())
    {
      QDPIO::cout << "LinOpSysSolverQPhiXClover:" << endl;
      QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << endl;


      QDPIO::cout << "Veclen is " << VecTraits<REALT>::Vec << endl;
      QDPIO::cout << "Soalen is " << VecTraits<REALT>::Soa << endl;
      if ( VecTraits<REALT>::Soa > VecTraits<REALT>::Vec ) { 
	QDPIO::cerr << "PROBLEM: Soalen > Veclen. Please set soalen appropriately (<=VECLEN) at compile time" << endl;
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

      QDPIO::cout << "About to grap a Dslash" << endl;
      geom = new QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>(Layout::subgridLattSize().slice(),
												      invParam.By, 
												      invParam.Bz, 
												      invParam.NCores, 
												      invParam.Sy,
												      invParam.Sz,
												      invParam.PadXY,
												      invParam.PadXYZ,
														   invParam.MinCt);
      
      QDPIO::cout << " Allocating p and c" << endl << flush ;
      p_even=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      p_odd=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      c_even=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      c_odd=(QPhiX_Spinor *)geom->allocCBFourSpinor();
      psi_s[0]=p_even;
      psi_s[1]=p_odd;
      chi_s[0]=c_even;
      chi_s[1]=c_odd;

      QDPIO::cout << " Allocating Clover" << endl << flush ;
      QPhiX_Clover* A_cb0=(QPhiX_Clover *)geom->allocCBClov();
      QPhiX_Clover* A_cb1=(QPhiX_Clover *)geom->allocCBClov();
      clov_packed[0] = A_cb0;
      clov_packed[1] = A_cb1;

      QDPIO::cout << " Allocating CloverInv " << endl << flush ;
      QPhiX_Clover* A_inv_cb0=(QPhiX_Clover *)geom->allocCBClov();
      QPhiX_Clover* A_inv_cb1=(QPhiX_Clover *)geom->allocCBClov();
      invclov_packed[0] = A_inv_cb0;
      invclov_packed[1] = A_inv_cb1;

      // Pack the gauge field
      QDPIO::cout << "Packing gauge field..."  << endl << flush ;
      QPhiX_Gauge* packed_gauge_cb0=(QPhiX_Gauge *)geom->allocCBGauge();
      QPhiX_Gauge* packed_gauge_cb1=(QPhiX_Gauge *)geom->allocCBGauge();

      QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, *geom);
      u_packed[0] = packed_gauge_cb0;
      u_packed[1] = packed_gauge_cb1;

      
      QDPIO::cout << "Creating Clover Term" << endl;
      QDPCloverTerm clov_qdp;
      clov->create(state_, invParam.CloverParams);
      QDPIO::cout << "Inverting Clover Term" << endl;
      invclov->create(state_, invParam.CloverParams, (*clov));
      for(int cb=0; cb < 2; cb++) { 
	invclov->choles(cb);
      }
      QDPIO::cout << "Done" << endl;
      QDPIO::cout << "Packing Clover term..." << endl;
      
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*invclov).getTriBuffer(), invclov_packed[cb], *geom, cb);
      }
    
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*clov).getTriBuffer(), clov_packed[cb], *geom, cb);
      }
      QDPIO::cout << "Done" << endl;

      QDPIO::cout << "Creating the Even Odd Operator" << endl;
      M=new QPhiX::EvenOddCloverOperator<REALT,VecTraits<REALT>::Vec,VecTraits<REALT>::Soa,VecTraits<REALT>::compress12>(u_packed,  
															 clov_packed[1], 
															 invclov_packed[0],  
															 geom,
															 t_boundary,
															 toDouble(aniso_coeffs[0]),
															 toDouble(aniso_coeffs[3]));
      
      
      switch( invParam.SolverType ) { 
      case BICGSTAB:
	{
	  QDPIO::cout << "Creating the BiCGStab Solver" << endl;
	  bicgstab_solver = new QPhiX::InvBiCGStab<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>((*M), invParam.MaxIter,1);

#if 1
	  if( invParam.TuneP ) bicgstab_solver->tune();
#endif

	}
	break;
      default:
	QDPIO::cerr << "UNKNOWN Solver Type" << endl;
	QDP_abort(1);
      }
    }

    

    //! Destructor is automatic
    ~LinOpSysSolverQPhiXClover() 
    {
      
      // Need to unalloc all the memory...
      QDPIO::cout << "Destructing" << endl;

      geom->free(p_even);
      geom->free(p_odd);
      geom->free(c_even);
      geom->free(c_odd);
      psi_s[0] = 0x0;
      psi_s[1] = 0x0;
      chi_s[0] = 0x0;
      chi_s[1] = 0x0;
      
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
	  res = cgSolve(psi,chi);
	}
	break;
      case BICGSTAB:
	{
	  res = biCGStabSolve(psi,chi);
	}
	break;
      default:
	QDPIO::cout << "Unknown Solver " << endl;
	break;
      }
      
      return res;
    }





  private:
    // Hide default constructor
    LinOpSysSolverQPhiXClover() {}
    


    Handle< LinearOperator<T> > A;
    const SysSolverQPhiXCloverParams invParam;
    Handle< QDPCloverTermT<T, U> > clov;
    Handle< QDPCloverTermT<T, U> > invclov;

    QPhiX::Geometry<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>* geom;
    
    Handle< QPhiX::EvenOddCloverOperator<REALT, VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > M;

    Handle< QPhiX::InvCG<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12> > cg_solver;

    Handle< QPhiX::InvBiCGStab<REALT,VecTraits<REALT>::Vec, VecTraits<REALT>::Soa, VecTraits<REALT>::compress12>  > bicgstab_solver;
    
    QPhiX_Clover* invclov_packed[2];
    QPhiX_Clover* clov_packed[2];
    QPhiX_Gauge* u_packed[2];
    
    QPhiX_Spinor* p_even;
    QPhiX_Spinor* p_odd;
    QPhiX_Spinor* c_even;
    QPhiX_Spinor* c_odd;
    QPhiX_Spinor* psi_s[2];
    QPhiX_Spinor* chi_s[2];
    
    SystemSolverResults_t cgSolve(T& psi, const T& chi) const
    {
#if 0
      SystemSolverResults_t res;
      psi = zero;

        // This is LinOpSolver so CGNE
      // Step 1 tranform RHS: chi -> M^\dagger chi
      // Move this into library later.
      T mdag_chi;
      (*A)(mdag_chi, chi, MINUS);
      
      
      //      QDPIO::cout << "Allocating Spinor fields" << endl;
      // Pack Spinors psi and chi
      QPhiX::qdp_pack_spinor<>(psi, psi_s[0], psi_s[1], *geom);
      QPhiX::qdp_pack_spinor<>(mdag_chi, chi_s[0], chi_s[1], *geom);
      
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      double start = omp_get_wtime();
      (*cg_solver)(psi_s[1],chi_s[1], res.n_count, rsd_final, site_flops, mv_apps, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << endl;      
      QPhiX::qdp_unpack_spinor<>(psi_s[0], psi_s[1], psi, (*M).getGeometry());

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: || r || / || b || = " << rel_resid << endl;

#if 0
      if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	QDPIO::cout << "SOLVE FAILED" << endl;
	QDP_abort(1);
      }
#endif

      int num_cb_sites = Layout::vol()/2;
      unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
      double gflops = (double)(total_flops)/(1.0e9);

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_CG_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << endl;

      END_CODE();
      return res;
#endif
    }

    SystemSolverResults_t biCGStabSolve(T& psi, const T& chi) const
    {
      SystemSolverResults_t res;
      psi=zero;

      // Pack Spinors psi and chi
      QDPIO::cout << "Packing" << endl << flush ;
      QPhiX::qdp_pack_spinor<>(psi, psi_s[0], psi_s[1], *geom);
      QPhiX::qdp_pack_spinor<>(chi, chi_s[0], chi_s[1], *geom);
      QDPIO::cout << "Done" << endl << flush;
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      QDPIO::cout << "Starting solve" << endl << flush ;
      double start = omp_get_wtime();
      (*bicgstab_solver)(psi_s[1],chi_s[1], toDouble(invParam.RsdTarget), res.n_count, rsd_final, site_flops, mv_apps, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << endl;      
      QPhiX::qdp_unpack_spinor<>(psi_s[0], psi_s[1], psi, *geom);

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: || r || / || b || = " << rel_resid << endl;

#if 0
      if ( !toBool (  rel_resid < invParam.RsdTarget*invParam.RsdToleranceFactor ) ) {
	QDPIO::cout << "SOLVE FAILED" << endl;
	QDP_abort(1);
      }
#endif

      int num_cb_sites = Layout::vol()/2;
      unsigned long total_flops = (site_flops + (1320+504+1320+504+48)*mv_apps)*num_cb_sites;
      double gflops = (double)(total_flops)/(1.0e9);

      double total_time = end - start;
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << endl;

      END_CODE();
      return res;
    }

  };


} // End namespace

#endif // BUILD_QPHIX
#endif 

