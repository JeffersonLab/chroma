// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_clover_intel_h__
#define __syssolver_linop_clover_intel_h__

#include "chroma_config.h"

#ifdef BUILD_ISOLVER
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/intel_solvers/syssolver_intel_clover_params.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "chroma_config.h"
#include "util/gauge/reunit.h"
#include "io/aniso_io.h"

#include "geometry.h"
#include "cpp_dslash_qdp_packer.h"
#include "clover.h"
#include "invcg.h"
#include "invbicgstab.h"

using namespace std;
using namespace CPlusPlusWilsonDslash;
using namespace QDP;

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverIntelCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }



  //! Solve a Clover Fermion System using the Intel inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */
  template<typename T, typename U>  
  class LinOpSysSolverIntelClover : public LinOpSystemSolver<T>
  {
  public:
    typedef multi1d<U> Q;
    typedef typename WordType<T>::Type_t REALT;
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverIntelClover(Handle< LinearOperator<T> > A_,
			      Handle< FermState<T,Q,Q> > state_,
			      const SysSolverIntelCloverParams& invParam_) : 
      A(A_), invParam(invParam_), clov(new QDPCloverTermT<T, U>()), invclov(new QDPCloverTermT<T, U>())
    {
      QDPIO::cout << "LinOpSysSolverIntelClover:" << endl;
      QDPIO::cout << "Compression is: " << invParam.CompressP << endl;
      QDPIO::cout << "AntiPeriodicT is: " << invParam.AntiPeriodicT << endl;


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
      
      float t_boundary=(float)(1);
      // NB: In this case, the state will have boundaries applied.
      // So we only need to apply our own boundaries if Compression is enabled
      //
      if (invParam.AntiPeriodicT) {
	t_boundary=(float)(-1);
	// Flip off the boundaries -- Dslash expects them off..
	u[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));
      }

      // Create the Clover Op
      QDPIO::cout << "Creating the Clover Op" << endl;
      M = new EvenOddCloverOperator<float, Veclen, Soalen>(Layout::subgridLattSize().slice(), 
							   invParam.By, 
							   invParam.Bz, 
							   invParam.NCores, 
							   invParam.Sy,
							   invParam.Sz,
							   invParam.PadXY,
							   invParam.PadXYZ,
							   invParam.MinCt, 
							   invParam.CompressP,
							   t_boundary,
							   toFloat(aniso_coeffs[0]),
							   toFloat(aniso_coeffs[3])
							   );
      
      QDPIO::cout << "Allocating Packed Gauge Field" << endl;
      ClovDslash<float,Veclen,Soalen>::SU3MatrixBlockF* packed_gauge_cb0 = (ClovDslash<float,Veclen,Soalen>::SU3MatrixBlockF*)(*M).allocCBGauge();
      ClovDslash<float,Veclen,Soalen>::SU3MatrixBlockF* packed_gauge_cb1 = (ClovDslash<float,Veclen,Soalen>::SU3MatrixBlockF*)(*M).allocCBGauge();

      printf("packed_gauge_cb0 = %lx\n", packed_gauge_cb0);
      printf("packed_gauge_cb1 = %lx\n", packed_gauge_cb1);


      QDPIO::cout << "Allocating Psi and Chi spinors" << endl;
      p_even=(ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF*)(*M).allocCBFourSpinor();
      p_odd=(ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF*)(*M).allocCBFourSpinor();
      c_even=(ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF*)(*M).allocCBFourSpinor();
      c_odd=(ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF*)(*M).allocCBFourSpinor();

      psi_s[0] = p_even + 1;
      psi_s[1] = p_odd + 1;

      chi_s[0] = c_even + 1;
      chi_s[1] = c_odd + 1;
      
      QDPIO::cout << "Allocate Packed Clover Term" << endl;
      ClovDslash<float, Veclen, Soalen>::CloverBlockF* A_cb0=(ClovDslash<float,Veclen,Soalen>::CloverBlockF*)(*M).allocCBClover();
      ClovDslash<float, Veclen, Soalen>::CloverBlockF* A_cb1=(ClovDslash<float,Veclen,Soalen>::CloverBlockF*)(*M).allocCBClover();
      ClovDslash<float, Veclen, Soalen>::CloverBlockF* A_inv_cb0=(ClovDslash<float,Veclen,Soalen>::CloverBlockF*)(*M).allocCBClover();
      ClovDslash<float, Veclen, Soalen>::CloverBlockF* A_inv_cb1=(ClovDslash<float,Veclen,Soalen>::CloverBlockF*)(*M).allocCBClover();
      
      invclov_packed[0] = A_inv_cb0;
      invclov_packed[1] = A_inv_cb1;
      clov_packed[0] = A_cb0;
      clov_packed[1] = A_cb1;
          
      // Clover stuff
      // Clover term uses u from state (with potentially antiperiodic BCs on)
      // However I am not doing anything regarding anisotropy, as QDP++ takes care of that.
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
	qdp_pack_clover<float,Veclen, Soalen, multi1d<PrimitiveClovTriang<float> > >((*invclov).getTriBuffer(), invclov_packed[cb], (*M).getGeometry(), cb);
      }
    
      for(int cb=0; cb < 2; cb++) { 
	qdp_pack_clover<float,Veclen, Soalen, multi1d<PrimitiveClovTriang<float> > >((*clov).getTriBuffer(), clov_packed[cb], (*M).getGeometry(), cb);
      }
      QDPIO::cout << "Done" << endl;


    
      QDPIO::cout << "Packing Gauge Field" << endl;
      qdp_pack_gauge<float,Veclen,Soalen, Q >(u,
	packed_gauge_cb0,
	packed_gauge_cb1, 
	(*M).getGeometry());
;
      QDPIO::cout << "Done" << endl;
      
      u_packed[0] = packed_gauge_cb0;
      u_packed[1] = packed_gauge_cb1;
      QDPIO::cout << "Done";


      (*M).setFields(u_packed, clov_packed[1], invclov_packed[0]);
      
      switch( invParam.SolverType ) { 
      case CG: 
	{
	  QDPIO::cout << "Creating the CG Solver" << endl;
	  cg_solver = new InvCG<float,Veclen, Soalen, ClovDslash<float,Veclen, Soalen>::FourSpinorBlockF>((*M), toFloat(invParam.RsdTarget), invParam.MaxIter);
	  if( invParam.TuneP ) cg_solver->tune();
	}
	break;
      case BICGSTAB:
	{
	  QDPIO::cout << "Creating the BiCGStab Solver" << endl;
	  bicgstab_solver = new CPlusPlusWilsonDslash::InvBiCGStab<float,Veclen, Soalen, ClovDslash<float,Veclen, Soalen>::FourSpinorBlockF>((*M), toFloat(invParam.RsdTarget), invParam.MaxIter);

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
    ~LinOpSysSolverIntelClover() 
    {
      // Need to unalloc all the memory...
      QDPIO::cout << "Destructing" << endl;
      (*M).free(p_even);
      (*M).free(p_odd);
      (*M).free(c_even);
      (*M).free(c_odd);
      psi_s[0] = 0x0;
      psi_s[1] = 0x0;
      chi_s[0] = 0x0;
      chi_s[1] = 0x0;
      
      (*M).free(invclov_packed[0]);
      (*M).free(invclov_packed[1]);
      invclov_packed[0] = 0x0;
      invclov_packed[1] = 0x0;
      
      (*M).free(clov_packed[0]);
      (*M).free(clov_packed[1]);
      clov_packed[0] = 0x0;
      clov_packed[1] = 0x0;
      
      (*M).free(u_packed[0]);
      (*M).free(u_packed[1]);
      u_packed[0] = 0x0;
      u_packed[1] = 0x0;

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
    LinOpSysSolverIntelClover() {}
    
    // Somehow, I need to take care of this...
    static const int Veclen=ISOLVER_VECLEN;
    static const int Soalen=ISOLVER_SOALEN;


    Handle< LinearOperator<T> > A;
    const SysSolverIntelCloverParams invParam;
    Handle< QDPCloverTermT<T, U> > clov;
    Handle< QDPCloverTermT<T, U> > invclov;

    Handle< EvenOddCloverOperator<float, Veclen, Soalen> > M;
    Handle< InvCG<float,Veclen, Soalen, ClovDslash<float,Veclen, Soalen>::FourSpinorBlockF> > cg_solver;
    Handle< CPlusPlusWilsonDslash::InvBiCGStab<float,Veclen, Soalen, ClovDslash<float,Veclen, Soalen>::FourSpinorBlockF > > bicgstab_solver;

    ClovDslash<float, Veclen, Soalen>::CloverBlockF* invclov_packed[2];
    ClovDslash<float, Veclen, Soalen>::CloverBlockF* clov_packed[2];
    ClovDslash<float,Veclen,Soalen>::SU3MatrixBlockF* u_packed[2];

    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* p_even;
    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* p_odd;
    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* c_even;
    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* c_odd;
    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* psi_s[2];
    ClovDslash<float,Veclen,Soalen>::FourSpinorBlockF* chi_s[2];

    SystemSolverResults_t cgSolve(T& psi, const T& chi) const
    {
      SystemSolverResults_t res;
      psi = zero;

        // This is LinOpSolver so CGNE
      // Step 1 tranform RHS: chi -> M^\dagger chi
      // Move this into library later.
      T mdag_chi;
      (*A)(mdag_chi, chi, MINUS);
      
      
      //      QDPIO::cout << "Allocating Spinor fields" << endl;
      // Pack Spinors psi and chi
      qdp_pack_spinor<float,Veclen,Soalen, T>(psi, psi_s[0], psi_s[1], (*M).getGeometry());
      qdp_pack_spinor< float, Veclen, Soalen,T >(mdag_chi, chi_s[0], chi_s[1], (*M).getGeometry());
      
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      double start = omp_get_wtime();
      (*cg_solver)(psi_s[1],chi_s[1], res.n_count, rsd_final, site_flops, mv_apps, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "INTEL_CLOVER_CG_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << endl;      
      qdp_unpack_spinor<float, Veclen, Soalen, T >(psi_s[0], psi_s[1], psi, (*M).getGeometry());

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "INTEL_CLOVER_CG_SOLVER: || r || / || b || = " << rel_resid << endl;

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
      QDPIO::cout << "INTEL_CLOVER_CG_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << endl;

      END_CODE();
      return res;
    }

    SystemSolverResults_t biCGStabSolve(T& psi, const T& chi) const
    {
      SystemSolverResults_t res;
      psi=zero;

      // Pack Spinors psi and chi
      QDPIO::cout << "PAcking" << endl << flush ;
      qdp_pack_spinor<float,Veclen,Soalen,T>(psi, psi_s[0], psi_s[1], (*M).getGeometry());
      qdp_pack_spinor<float,Veclen,Soalen,T>(chi, chi_s[0], chi_s[1], (*M).getGeometry());
      QDPIO::cout << "Done" << endl << flush;
      double rsd_final;
      unsigned long site_flops=0;
      unsigned long mv_apps=0;
      
      QDPIO::cout << "Starting solve" << endl << flush ;
      double start = omp_get_wtime();
      (*bicgstab_solver)(psi_s[1],chi_s[1], 1, res.n_count, rsd_final, site_flops, mv_apps, invParam.VerboseP);
      double end = omp_get_wtime();

      QDPIO::cout << "INTEL_CLOVER_BICGSTAB_SOLVER: " << res.n_count << " iters,  rsd_sq_final=" << rsd_final << endl;      
      qdp_unpack_spinor<float, Veclen, Soalen, T >(psi_s[0], psi_s[1], psi, (*M).getGeometry());

      // Chi Should now hold the result spinor 
      // Check it against chroma.
      T r = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[ A->subset() ] -= tmp;

      Double r2 = norm2(r,A->subset());
      Double b2 = norm2(chi, A->subset());
      Double rel_resid = sqrt(r2/b2);

      QDPIO::cout << "INTEL_CLOVER_BICGSTAB_SOLVER: || r || / || b || = " << rel_resid << endl;

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
      QDPIO::cout << "INTEL_CLOVER_BICGTAB_SOLVER: Solver Time="<< total_time <<" (sec)  Performance=" << gflops / total_time << " GFLOPS" << endl;

      END_CODE();
      return res;
    }

  };


} // End namespace

#endif // BUILD_Intel
#endif 

