// -*- C++ -*-
// $Id: eoprec_dwf_qprop_array_cg_dwf_lowmem_w.h,v 3.2 2007-03-05 19:36:32 bjoo Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chroma_config.h"
#if defined(CHROMA_USE_CG_DWF_LOWMEM)


#ifndef __prec_dwf_qprop_array_cg_dwf_lowmem_w_h__
#define __prec_dwf_qprop_array_cg_dwf_lowmem_w_h__


#include "eoprec_constdet_wilstype_fermact_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


namespace Chroma
{

  //! SSE Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  template< typename SinglePrecSolver, typename DoublePrecSolver >
  class CGDWFQpropT : public SystemSolverArray<LatticeFermion>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;


    //! Alternative constructor for compatibility
    /*!
     * \param m_q_       quark mass ( Read )
     */
    CGDWFQpropT(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
		 Handle< LinOpSystemSolverArray<T> > invA_,   // throw away
		 Handle< FermState<T,P,Q> > state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 const AnisoParam_t& anisoParam_,
		 const GroupXML_t& invParam_) : 
      A(A_), OverMass(OverMass_), Mass(Mass_), 
      N5(A->size()), anisoParam(anisoParam_)
      {init(state_, invParam_);}

    //! Need a real destructor
    ~CGDWFQpropT() {fini();}

    //! Expected length of array index
    int size() const {return N5;}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const
    {
      
      QDPIO::cout << "entering CGDWFQpropT::operator()" << endl;
      
      START_CODE();
      
      SystemSolverResults_t res;
      
      //    init();   // only needed because 2 qpropT might be active - SSE CG does not allow this
      
      Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
      
      // Apply SSE inverter
      Real   a5  = 1;
      double M5  = toDouble(1 + a5*(1 + (Nd-1)*ff - OverMass));
      double m_f = toDouble(Mass);
      double out_eps;
      int single_count = 0;
      int double_count = 0;
      res.n_count = 0;

      StopWatch swatch, swatch_total, swatch_solver;
      {
#ifdef SINGLE_PREC_SOLVER
	double rsd = toDouble(invParam.RsdCG);
	double rsd_sq = rsd * rsd;
	int    max_iter = invParam.MaxCG;

	swatch_total.reset(); swatch_solver.reset(); 
	swatch_total.start();

	// Init solver
	std::string dwf_error_str = "single prec.";

	QDPIO::cout << "CGDWFQpropT: Initializing Single Precision Solver " << endl;
	swatch.reset(); swatch.start();
	int stat = single_prec_solver.init(lattice_size.slice(), NULL, NULL);
	if ( stat != 0) {
	  
	  QDPIO::cerr << __func__ << ": error in SP solver init: " << dwf_error_str << " error code is " << stat << endl;
	  QDP_abort(1);
	}
	swatch.stop();
	QDPIO::cout << "CGDWFQpropT: Single Prec Solver init took: " << swatch.getTimeInSeconds() << " seconds " << endl;
	
	// load the gauge
	swatch.reset(); swatch.start();
	single_prec_solver.loadGauge(&u, &v);
	swatch.stop();
	QDPIO::cout << "CGDWFQpropT: Single Precision Gauge Field Import took: " << swatch.getTimeInSeconds() << " seconds " << endl;

	
	// Call the solver
	QDPIO::cout << "CGDWFQpropT: Beginning Single Precision Solve" << endl;
	swatch_solver.start();
	single_prec_solver.cgSolver(psi, M5, m_f, 
				    chi, psi, rsd_sq, max_iter, out_eps, single_count);
	swatch_solver.stop();
	// Delete the gauge
	QDPIO::cout << "CGDWFQpropT: Deleting Single Precision Solver Gauge Field " << endl;
	single_prec_solver.deleteGauge();
	
	QDPIO::cout << "CGDWFQPropT: Finalizing Single Precision Solver " << endl;
	// Destroy the solver
	single_prec_solver.fini();
	swatch_total.stop();
	double tot_time = swatch_total.getTimeInSeconds();
	double sol_time = swatch_solver.getTimeInSeconds();

	QDPIO::cout << "CGDWFQPropT: Total solution time: " << tot_time << " Time in Single Prec Solver: " << sol_time << " Single Prec Overhead: " << 100.0*(tot_time-sol_time)/sol_time << "%" << endl;

	res.n_count += single_count;
      }
#endif
#ifdef DOUBLE_PREC_SOLVER
      {

	double rsd = toDouble(invParam.RsdCGRestart);
	double rsd_sq = rsd * rsd;
	int    max_iter = invParam.MaxCGRestart;

	swatch_total.reset(); swatch_solver.reset(); 
	swatch_total.start();
	
	// Init the solver
	std::string dwf_error_str2 = "double prec.";
	
	QDPIO::cout << "CGDWFQpropT: Initializing Double Precision Solver " << endl;
	swatch.reset(); swatch.start();
	int stat2 = double_prec_solver.init(lattice_size.slice(), NULL, NULL);
	if ( stat2 != 0) {
	  QDPIO::cerr << __func__ << ": error in DP solver init: " << dwf_error_str2 << " error code is " << stat2 << endl;
	  QDP_abort(1);
	}
	swatch.stop();
	QDPIO::cout << "CGDWFQpropT: Double Prec Solver init took: " << swatch.getTimeInSeconds() << " seconds " << endl;

	// Load the gauge field
	swatch.reset(); swatch.start();
	double_prec_solver.loadGauge(&u, &v);
	swatch.stop();
	QDPIO::cout << "Double Precision Gauge Field Import took: " << swatch.getTimeInSeconds() << " seconds " << endl;
	
	// Do the solve
	QDPIO::cout << "CGDWFQpropT: Beginning Double Precision Solve" << endl;
	swatch_solver.start();
	double_prec_solver.cgSolver(psi, M5, m_f, 
				    chi, psi, rsd_sq, max_iter, out_eps, double_count);
	swatch_solver.stop();
	
	// Delete the gauge field
	double_prec_solver.deleteGauge();
	
	// Kill the solver
	double_prec_solver.fini();
	swatch_total.stop();
	double tot_time = swatch_total.getTimeInSeconds();
	double sol_time = swatch_solver.getTimeInSeconds();

	QDPIO::cout << "CGDWFQPropT: Total solution time: " << tot_time << " Time in Double Prec Solver: " << sol_time << " Double Prec Overhead: " << 100.0*(tot_time-sol_time)/sol_time << "%" << endl;
      
	
	res.n_count += double_count;
      }
#endif
      QDPIO::cout << "CGDWFQpropT: Single Prec. Iters = " << single_count << " Double Prec. Iters = " << double_count << " Total Iters = " << res.n_count << endl;
      
      {
	multi1d<LatticeFermion>  r(N5);
	A->unprecLinOp(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }
      
      QDPIO::cout << "exiting CGDWFQpropT::operator()" << endl;
      
      END_CODE();
      
      return res;
    }

  protected:
    //! Private internal initializer
    void init(Handle< FermState<T,P,Q> > state, const GroupXML_t& inv)
    {
      QDPIO::cout << "entering CGDWFQpropT::init" << endl;
      
      if (Nd != 4 || Nc != 3) {
	
	QDPIO::cerr << "CGDWFQpropT: only supports Nd=4 and Nc=3" << endl;
	QDP_abort(1);
	
      }

      // Read the XML for the CG params
      try {
	
	std::istringstream  is(inv.xml);
	XMLReader  paramtop(is);
	
	read(paramtop, inv.path, invParam);
      }
      catch (const std::string& e) {
	
	QDPIO::cerr << "CGDWFQpropT: only support a CG inverter" << endl;
	QDP_abort(1);
      }

      lattice_size.resize(Nd+1);
      lattice_size[Nd] = N5;

      for(int i=0; i < Nd; ++i)
	lattice_size[i] = Layout::lattSize()[i];
      
      if (N5 % 4 != 0) {
	
	QDPIO::cerr << "SSE qpropT only N5 that is multiple of 2" << endl;
	QDP_abort(1);
      }
   
      u.resize(Nd);
      u = state->getLinks();
      
      if (anisoParam.anisoP) {
	
	Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
	
	// Rescale the u fields by the anisotropy
	for(int mu=0; mu < u.size(); ++mu) {
	  
	  if (mu != anisoParam.t_dir)
	    u[mu] *= ff;
	}
      }
      
      // Construct shifted gauge field
      v.resize(Nd);
      for (int i = 0; i < Nd; i++)
	v[i] = shift(u[i], -1, i); // as viewed from the destination

      QDPIO::cout << "exiting CGDWFQpropT::init" << endl;
    }

    //! Private internal destructor
    void fini() {
      QDPIO::cout << "CGDWFQpropT: calling destructor" << endl;
    }
      
  private:
#ifdef SINGLE_PREC_SOLVER
    mutable SinglePrecSolver  single_prec_solver;
#endif

#ifdef DOUBLE_PREC_SOLVER
    mutable DoublePrecSolver double_prec_solver;
#endif

    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Real OverMass;
    Real Mass;
    int  N5;
    AnisoParam_t anisoParam;
    SysSolverCGParams invParam;
    multi1d<int> lattice_size;
    multi1d<LatticeColorMatrix> u;
    multi1d<LatticeColorMatrix> v;
  };

}
#endif

#endif // CHROMA_USE_CG_DWF_LOWMEM
