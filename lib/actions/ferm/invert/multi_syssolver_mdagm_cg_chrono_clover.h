// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2 using CG
 */

#ifndef __multi_syssolver_mdagm_cg_chrono_clover_h__
#define __multi_syssolver_mdagm_cg_chrono_clover_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/linop/lopishift.h"
#include "update/molecdyn/predictor/mre_shifted_predictor.h"
#include "init/chroma_init.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/invert/reliable_cg.h"
#include "actions/ferm/invert/minvcg2.h"
#include "lmdagm.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGChronoCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }

  struct MultiSysSolverCGChronoCloverParams { 
    MultiSysSolverCGChronoCloverParams(XMLReader& xml, const std::string& path);
    MultiSysSolverCGChronoCloverParams() {};
    MultiSysSolverCGChronoCloverParams( const MultiSysSolverCGChronoCloverParams& p) {
      clovParams = p.clovParams;
      MaxIter = p.MaxIter;
      MaxChrono = p.MaxChrono;
      Delta = p.Delta;
      CutoffRsd = p.CutoffRsd;
      RsdTarget = p.RsdTarget; // Array Copy
    }
    CloverFermActParams clovParams;
    int MaxIter;
    int MaxChrono;

    Real Delta;
    Real CutoffRsd;
    multi1d<Real> RsdTarget;

  };
  
#if 1
  void read(XMLReader& xml, const std::string& path, MultiSysSolverCGChronoCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverCGChronoCloverParams& param);
#endif

  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  class MdagMMultiSysSolverCGChronoClover : public MdagMMultiSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;
    typedef multi1d<LatticeColorMatrix> P;
 
    typedef LatticeFermionF TF;
    typedef LatticeColorMatrixF UF;
    typedef multi1d<LatticeColorMatrixF> QF;
    typedef multi1d<LatticeColorMatrixF> PF;

    typedef LatticeFermionD TD;
    typedef LatticeColorMatrixD UD;
    typedef multi1d<LatticeColorMatrixD> QD;
    typedef multi1d<LatticeColorMatrixD> PD;

    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGChronoClover(Handle< LinearOperator<T> > M_,
				      Handle< FermState<T,P,Q> > state_,
				      const MultiSysSolverCGChronoCloverParams& invParam_) : 
      M(M_), invParam(invParam_) 
    {
      QF links_single; links_single.resize(Nd);
      QD links_double; links_double.resize(Nd);
      
      const Q& links = state_->getLinks();
      for(int mu=0; mu < Nd; mu++) { 
	links_single[mu] = links[mu];
	links_double[mu] = links[mu];
      }
      
      // Links single hold the possibly stouted links
      // with gaugeBCs applied... 
      // Now I need to create a single prec state...
      fstate_single = new PeriodicFermState<TF,QF,QF>(links_single);
      fstate_double = new PeriodicFermState<TD,QD,QD>(links_double);
      
    }

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGChronoClover() {}
    
    //! Return the subset on which the operator acts
    const Subset& subset() const {return M->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<Real>& shifts, const T& chi) const
    {
      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      SystemSolverResults_t res;
      res.n_count = 0;

      multi1d<Real> RsdCG(shifts.size());
      if (invParam.RsdTarget.size() == 1) {
	RsdCG = invParam.RsdTarget[0];
      }
      else if (invParam.RsdTarget.size() == RsdCG.size()) {
	
	RsdCG = invParam.RsdTarget;
      }
      else {
	
	QDPIO::cerr << "MdagMMultiSysSolverCGChronoClover: shifts incompatible" << endl;
	QDP_abort(1);
      }


      // Setup residua to Single Prec CG Solver
      multi1d<RealF> modRsdCG(shifts.size());
      for(int i=0; i < shifts.size(); i++) { 
	if( toBool(  RsdCG[i] < invParam.CutoffRsd ) ) { 
	  modRsdCG[i] = invParam.CutoffRsd;
	}
	else {
	  modRsdCG[i] = RsdCG[i];
	}
      }
	

      
      
      Handle< LinearOperator<TF> > M_single(new EvenOddPrecDumbCloverFLinOp( fstate_single, invParam.clovParams ));
      multi1d<TF> psi_f(shifts.size());
      {

	TF chi_f;
	
	// Single Prec RHC
	chi_f[M->subset()] = chi;
	for(int i=0; i < shifts.size(); i++) { 
	  psi_f[i][M->subset()] = zero;
	}
	
	// Setup single prec shifts
	multi1d<RealF> shifts_r(shifts.size());
	for(int i=0; i < shifts.size(); i++ ) {
	  shifts_r[i] = shifts[i];
	}
	
	MInvCG2(*M_single,
		chi_f, 
		psi_f,
		shifts_r, 
		modRsdCG,
		invParam.MaxIter,
		res.n_count);
	
      }


      psi.resize(shifts.size());

      Handle< LinearOperator<TD> > M_double(new EvenOddPrecDumbCloverDLinOp( fstate_double, invParam.clovParams ));

      Handle<LinearOperator<TD> > MM( new MdagMLinOp<TD>(M_double));
      
      MinimalResidualExtrapolationShifted4DChronoPredictor<TD,RealD> chrono(invParam.MaxChrono,*MM);

      chrono.reset();
      TD chi_d; chi_d[M->subset()] = chi;
      TD psi_d;
      for(int i=shifts.size()-1; i >= 0; i--) { 
	RealD rshift = sqrt(shifts[i]);
	RealF rshift_s = rshift;
	RealD dshift = shifts[i];
	SystemSolverResults_t res_tmp;
	Handle< LinearOperator<TD> > Ms(new lopishift<TD,RealD>(M_double, rshift));
	Handle< LinearOperator<TF> > Ms_single( new lopishift<TF,RealF>(M_single, rshift_s) );
	

	psi_d[M->subset()] = psi_f[i];

	chrono.predictX(psi_d, dshift, chi_d);
	
	res_tmp= InvCGReliable(*(Ms),*(Ms_single), chi_d, psi_d, RsdCG[i],invParam.Delta, invParam.MaxIter);
	
	chrono.newXVector(psi_d);
	psi[i][M->subset()] = psi_d;
	QDPIO::cout << "n_count = " << res_tmp.n_count << endl;
	res.n_count += res_tmp.n_count;
      }
      

#if 0
      XMLFileWriter& log = Chroma::getXMLLogInstance();
      push(log, "MultiCG");
      write(log, "shifts", shifts);
      write(log, "RsdCG", RsdCG);
      write(log, "n_count", res.n_count);
      Double chinorm=norm2(chi, M->subset());
      multi1d<Double> r_rel(shifts.size());
      
      for(int i=0; i < shifts.size(); i++) { 
	T tmp1,tmp2;
	(*M)(tmp1, psi[i], PLUS);
	(*M)(tmp2, tmp1, MINUS);  // tmp2 = A^\dagger A psi
	tmp2[ M->subset() ] +=  shifts[i]* psi[i]; // tmp2 = ( A^\dagger A + shift_i ) psi
	T r;
	r[ M->subset() ] = chi - tmp2;
	r_rel[i] = sqrt(norm2(r, M->subset())/chinorm );
#if 1
	QDPIO::cout << "r[" <<i <<"] = " << r_rel[i] << endl;
#endif

      }
      write(log, "ResidRel", r_rel);
      pop(log);
#endif
      swatch.stop();
      double time = swatch.getTimeInSeconds();
      QDPIO::cout << "MULTI_CG_CHRONO_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << endl;
 QDPIO::cout << "MULTI_CG_CHRONO_CLOVER_SOLVER: "<<time<< " sec" << endl;
      END_CODE();
      
      return res;
    }

    
  private:
    // Hide default constructor
    MdagMMultiSysSolverCGChronoClover() {}

    Handle< LinearOperator<T> > M;
    const MultiSysSolverCGChronoCloverParams invParam;
    Handle< FermState<TF, QF, QF> > fstate_single;
    Handle< FermState<TD, QD, QD> > fstate_double;
    
  };


} // End namespace

#endif 

