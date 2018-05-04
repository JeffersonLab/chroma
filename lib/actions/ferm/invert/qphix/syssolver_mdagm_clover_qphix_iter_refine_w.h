// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_clover_qphix_iter_refine_h__
#define __syssolver_linop_clover_qphix_iter_refine_h__

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
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"

using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQPhiXCloverIterRefineEnv
 {
    //! Register the syssolver
    bool registerAll();

    template<typename TOuter,typename TInner>
    struct MixedVecTraits { 
      static const int Vec=1;
      static const int Soa=1;
      static const bool compress12=false;
      static const int VecInner=1;
      static const int SoaInner=1;

    };

    // Templates
#if defined CHROMA_QPHIX_ARCH_AVX
#warning QPhix Solver AVX
    // AVX Traits:
    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=4;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=4;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };

   template<>
    struct MixedVecTraits<float,float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };


#endif

#if defined CHROMA_QPHIX_ARCH_MIC
#warning QPhix solver MIC
    // MIC Traits
    template<>
    struct MixedVecTraits<double,double> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=8;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };
    template<>
    struct MixedVecTraits<double,float> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };
    template<>
    struct MixedVecTraits<double,QPhiX::half> { 
      static const int Vec=8;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };
    template<>
    struct MixedVecTraits<float,float> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;

    };
    template<>
    struct MixedVecTraits<float,QPhiX::half> { 
      static const int Vec=16;
      static const int Soa=CHROMA_QPHIX_SOALEN;
      static const bool compress12=CHROMA_QPHIX_COMPRESS12; 
      static const int VecInner=16;
      static const int SoaInner=CHROMA_QPHIX_INNER_SOALEN;
    };
#endif


  }



  //! Solve a Clover Fermion System using the QPhiX inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
  using namespace LinOpSysSolverQPhiXCloverIterRefineEnv;
  template<typename T, typename U>  
  class LinOpSysSolverQPhiXCloverIterRefine : public LinOpSystemSolver<T>
  {
  public:
    typedef multi1d<U> Q;
    typedef typename WordType<T>::Type_t REALT;
    typedef CHROMA_QPHIX_INNER_TYPE  InnerReal;

    typedef typename QPhiX::Geometry<REALT,MixedVecTraits<REALT,InnerReal>::Vec,MixedVecTraits<REALT,InnerReal>::Soa,MixedVecTraits<REALT,InnerReal>::compress12>::FourSpinorBlock QPhiX_Spinor;
    typedef typename QPhiX::Geometry<REALT,MixedVecTraits<REALT,InnerReal>::Vec,MixedVecTraits<REALT,InnerReal>::Soa,MixedVecTraits<REALT,InnerReal>::compress12>::SU3MatrixBlock QPhiX_Gauge;
    typedef typename QPhiX::Geometry<REALT,MixedVecTraits<REALT,InnerReal>::Vec,MixedVecTraits<REALT,InnerReal>::Soa,MixedVecTraits<REALT,InnerReal>::compress12>::CloverBlock  QPhiX_Clover;

    typedef typename QPhiX::Geometry<InnerReal,MixedVecTraits<REALT,InnerReal>::VecInner,
			                 MixedVecTraits<REALT,InnerReal>::SoaInner,
                                         MixedVecTraits<REALT,InnerReal>::compress12>::FourSpinorBlock QPhiX_InnerSpinor;

  typedef typename QPhiX::Geometry<InnerReal,MixedVecTraits<REALT,InnerReal>::VecInner,
			                 MixedVecTraits<REALT,InnerReal>::SoaInner,
                                         MixedVecTraits<REALT,InnerReal>::compress12>::SU3MatrixBlock QPhiX_InnerGauge;

    typedef typename QPhiX::Geometry<InnerReal,MixedVecTraits<REALT,InnerReal>::VecInner,
			                 MixedVecTraits<REALT,InnerReal>::SoaInner,
                                         MixedVecTraits<REALT,InnerReal>::compress12>::CloverBlock QPhiX_InnerClover;

    
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverQPhiXCloverIterRefine(Handle< LinearOperator<T> > A_,
			      Handle< FermState<T,Q,Q> > state_,
			      const SysSolverQPhiXCloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new QDPCloverTermT<T, U>()), invclov(new QDPCloverTermT<T, U>())
    {

      
      QDPIO::cout << "LinOpSysSolverQPhiXCloverIterRefine:" << std::endl;

      if ( toBool( invParam.Delta < 0 ) ) { 
	QDPIO::cout << "Error: Delta parameter for solve should be set" << std::endl;
	QDP_abort(1);
      }

      QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << std::endl;


      QDPIO::cout << "Veclen is " << MixedVecTraits<REALT,InnerReal>::Vec << std::endl;
      QDPIO::cout << "Soalen is " << MixedVecTraits<REALT,InnerReal>::Soa << std::endl;
      QDPIO::cout << "Inner Veclen is " << MixedVecTraits<REALT,InnerReal>::VecInner << std::endl;
      QDPIO::cout << "Inner Soalen is " << MixedVecTraits<REALT,InnerReal>::SoaInner << std::endl;
      
      if ( MixedVecTraits<REALT,InnerReal>::Soa > MixedVecTraits<REALT,InnerReal>::Vec ) { 
	QDPIO::cerr << "PROBLEM: Soalen > Veclen. Please set soalen appropriately (<=VECLEN) at compile time" << std::endl;
	QDP_abort(1);
      }

      if ( MixedVecTraits<REALT,InnerReal>::SoaInner > MixedVecTraits<REALT,InnerReal>::VecInner ) { 
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


      // Grab a dslash from which we will get geometry.
      geom_outer = new QPhiX::Geometry<REALT, 
	MixedVecTraits<REALT,InnerReal>::Vec, 
	MixedVecTraits<REALT,InnerReal>::Soa,
	MixedVecTraits<REALT,InnerReal>::compress12>(Layout::subgridLattSize().slice(),
						     invParam.By, 
						     invParam.Bz, 
						     invParam.NCores, 
						     invParam.Sy,
						     invParam.Sz,
						     invParam.PadXY,
						     invParam.PadXYZ,
						     invParam.MinCt);
      
      geom_inner = new QPhiX::Geometry<InnerReal, 
	MixedVecTraits<REALT,InnerReal>::VecInner, 
	MixedVecTraits<REALT,InnerReal>::SoaInner,
	MixedVecTraits<REALT,InnerReal>::compress12>(Layout::subgridLattSize().slice(),
						     invParam.By, 
						     invParam.Bz, 
						     invParam.NCores, 
						     invParam.Sy,
						     invParam.Sz,
						     invParam.PadXY,
						     invParam.PadXYZ,
						     invParam.MinCt);
						     


      p_even=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      p_odd=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      c_even=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      c_odd=(QPhiX_Spinor *)geom_outer->allocCBFourSpinor();
      psi_s[0]=p_even;
      psi_s[1]=p_odd;
      chi_s[0]=c_even;
      chi_s[1]=c_odd;


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
      QDPCloverTerm clov_qdp;
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
	QPhiX::qdp_pack_clover<>((*invclov).getTriBuffer(), invclov_packed[cb], *geom_outer, cb);
      }
    
      // Pack outer clover
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*clov).getTriBuffer(), clov_packed[cb], *geom_outer, cb);
      }

      // Pack inner clover inverse
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*invclov).getTriBuffer(), invclov_packed_i[cb], *geom_inner, cb);
      }
    
      // Pack inner clover
      for(int cb=0; cb < 2; cb++) { 
	QPhiX::qdp_pack_clover<>((*clov).getTriBuffer(), clov_packed_i[cb], *geom_inner, cb);
      }
      QDPIO::cout << "Done" << std::endl;


      M_outer=new QPhiX::EvenOddCloverOperator<REALT,
	MixedVecTraits<REALT,InnerReal>::Vec,
	MixedVecTraits<REALT,InnerReal>::Soa,
	MixedVecTraits<REALT,InnerReal>::compress12>(u_packed,  
						     clov_packed[1], 
						     invclov_packed[0],  
						     geom_outer,
						     t_boundary,
						     toDouble(aniso_coeffs[0]),
						     toDouble(aniso_coeffs[3]));
    
    M_inner=new QPhiX::EvenOddCloverOperator<InnerReal,
	MixedVecTraits<REALT,InnerReal>::VecInner,
	MixedVecTraits<REALT,InnerReal>::SoaInner,
	MixedVecTraits<REALT,InnerReal>::compress12>(u_packed_i,  
						     clov_packed_i[1], 
						     invclov_packed_i[0],  
						     geom_inner,
						     t_boundary,
						     toDouble(aniso_coeffs[0]),
						     toDouble(aniso_coeffs[3]));
    
    
    bicgstab_inner_solver = new QPhiX::InvBiCGStab<InnerReal,MixedVecTraits<REALT,InnerReal>::VecInner, MixedVecTraits<REALT,InnerReal>::SoaInner, MixedVecTraits<REALT,InnerReal>::compress12>((*M_inner), invParam.MaxIter,1);
    
    mixed_solver =   new QPhiX::InvRichardsonMultiPrec<REALT,
	MixedVecTraits<REALT,InnerReal>::Vec,
	MixedVecTraits<REALT,InnerReal>::Soa,
	MixedVecTraits<REALT,InnerReal>::compress12,
	InnerReal,
	MixedVecTraits<REALT,InnerReal>::VecInner,
	MixedVecTraits<REALT,InnerReal>::SoaInner,
	// NB: &(*) is yukcy. 
	MixedVecTraits<REALT,InnerReal>::compress12>((*M_outer), (*bicgstab_inner_solver), toDouble(invParam.Delta), 1, invParam.MaxIter);
    }
    
    

    //! Destructor 
    ~LinOpSysSolverQPhiXCloverIterRefine() 
    {
      
      // Need to unalloc all the memory...
      QDPIO::cout << "Destructing" << std::endl;
      geom_outer->free(p_even);
      geom_outer->free(p_odd);
      geom_outer->free(c_even);
      geom_outer->free(c_odd);
      psi_s[0] = 0x0;
      psi_s[1] = 0x0;
      chi_s[0] = 0x0;
      chi_s[1] = 0x0;

      geom_outer->free(invclov_packed[0]);
      geom_outer->free(invclov_packed[1]);
      invclov_packed[0] = 0x0;
      invclov_packed[1] = 0x0;


      geom_outer->free(clov_packed[0]);
      geom_outer->free(clov_packed[1]);
      clov_packed[0] = 0x0;
      clov_packed[1] = 0x0;


      geom_outer->free(u_packed[0]);
      geom_outer->free(u_packed[1]);
      u_packed[0] = 0x0;
      u_packed[1] = 0x0;


      geom_inner->free(invclov_packed_i[0]);
      geom_inner->free(invclov_packed_i[1]);
      invclov_packed_i[0] = 0x0;
      invclov_packed_i[1] = 0x0;


      geom_inner->free(clov_packed_i[0]);
      geom_inner->free(clov_packed_i[1]);
      clov_packed_i[0] = 0x0;
      clov_packed_i[1] = 0x0;


      geom_inner->free(u_packed_i[0]);
      geom_inner->free(u_packed_i[1]);
      u_packed_i[0] = 0x0;
      u_packed_i[1] = 0x0;

      delete geom_inner;      
      delete geom_outer;


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
      QDPIO::cout << "Packing Spinors" << std::endl << std::flush ;

      QPhiX::qdp_pack_spinor<>(psi, psi_s[0], psi_s[1], (*M_outer).getGeometry());
      QPhiX::qdp_pack_spinor<>(chi, chi_s[0], chi_s[1], (*M_outer).getGeometry());

      QDPIO::cout << "Done" << std::endl << std::flush;

      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      QDPIO::cout << "Starting solve" << std::endl << std::flush ;
      double start = omp_get_wtime();
      (*mixed_solver)(psi_s[1],chi_s[1], toDouble(invParam.RsdTarget), res.n_count, rsd_final, site_flops, mv_apps, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_ITER_REFINE_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << std::endl;      
      QPhiX::qdp_unpack_spinor<>(psi_s[0], psi_s[1], psi, (*M_outer).getGeometry());
      
#if 1
      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_ITER_REFINE_SOLVER: || r || / || b || = " << rel_resid << std::endl;
#endif
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
      QDPIO::cout << "QPHIX_CLOVER_BICGSTAB_ITER_REFINE_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << std::endl;

      END_CODE();
      return res;
    }





  private:
    // Hide default constructor
    LinOpSysSolverQPhiXCloverIterRefine() {}
    


    Handle< LinearOperator<T> > A;
    const SysSolverQPhiXCloverParams invParam;
    Handle< QDPCloverTermT<T, U> > clov;
    Handle< QDPCloverTermT<T, U> > invclov;

    QPhiX::Geometry<REALT, 
		    MixedVecTraits<REALT,InnerReal>::Vec, 
		    MixedVecTraits<REALT,InnerReal>::Soa, 
		    MixedVecTraits<REALT,InnerReal>::compress12>* geom_outer;

    QPhiX::Geometry<InnerReal, 
		    MixedVecTraits<REALT,InnerReal>::VecInner, 
		    MixedVecTraits<REALT,InnerReal>::SoaInner, 
		    MixedVecTraits<REALT,InnerReal>::compress12>* geom_inner;

    Handle< QPhiX::EvenOddCloverOperator<REALT, 
					 MixedVecTraits<REALT,InnerReal>::Vec, 
					 MixedVecTraits<REALT,InnerReal>::Soa, 
					 MixedVecTraits<REALT,InnerReal>::compress12> > M_outer;

Handle< QPhiX::EvenOddCloverOperator<InnerReal, 
				     MixedVecTraits<REALT,InnerReal>::VecInner, 
				     MixedVecTraits<REALT,InnerReal>::SoaInner, 
				     MixedVecTraits<REALT,InnerReal>::compress12> > M_inner;
    
    
    Handle< QPhiX::InvBiCGStab<InnerReal,
			       MixedVecTraits<REALT,InnerReal>::VecInner, 
			       MixedVecTraits<REALT,InnerReal>::SoaInner, 
			       MixedVecTraits<REALT,InnerReal>::compress12>  > bicgstab_inner_solver;

Handle< QPhiX::InvRichardsonMultiPrec<REALT,
				      MixedVecTraits<REALT,InnerReal>::Vec,
				      MixedVecTraits<REALT,InnerReal>::Soa,
				      MixedVecTraits<REALT,InnerReal>::compress12,
				      InnerReal,
				      MixedVecTraits<REALT,InnerReal>::VecInner,
				      MixedVecTraits<REALT,InnerReal>::SoaInner,
				      MixedVecTraits<REALT,InnerReal>::compress12> > mixed_solver;
    
    QPhiX_Clover* invclov_packed[2];
    QPhiX_Clover* clov_packed[2];
    QPhiX_Gauge* u_packed[2];
    
    QPhiX_InnerClover* invclov_packed_i[2];
    QPhiX_InnerClover* clov_packed_i[2];
    QPhiX_InnerGauge* u_packed_i[2];
    
    QPhiX_Spinor* p_even;
    QPhiX_Spinor* p_odd;
    QPhiX_Spinor* c_even;
    QPhiX_Spinor* c_odd;
    QPhiX_Spinor* psi_s[2];
    QPhiX_Spinor* chi_s[2];
    
  };


} // End namespace

#endif // BUILD_QPHIX
#endif 

