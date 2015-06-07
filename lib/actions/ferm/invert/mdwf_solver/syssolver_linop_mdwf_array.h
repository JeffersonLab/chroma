// -*- C++ -*-
/*! \file
 *  \brief DWF/SSE double-prec solver
 */

#ifndef SYSSOLVER_LINOP_MDWF_ARRAY_H
#define SYSSOLVER_LINOP_MDWF_ARRAY_H


extern "C" {
  struct QOP_MDWF_State;
  struct QOP_MDWF_Parameters;  
};

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/mdwf_solver/syssolver_mdwf_params.h"


using namespace QDP;
namespace Chroma 
{ 

  //! CG1 system solver namespace
  namespace LinOpSysSolverMDWFArrayEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! AVP's DWF Solver interface
  /*!
   * \ingroup invert
   *
   * @{
   */
  class LinOpSysSolverMDWFArray : public LinOpSystemSolverArray<LatticeFermion>  { 
    
  public:
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    /* Constructor */
    LinOpSysSolverMDWFArray(
			    Handle< LinearOperatorArray<T> > A_,
			    Handle< FermState<T,P,Q> > fs_, 
			    const SysSolverMDWFParams& invParam_ )
      : A(A_), invParam(invParam_) {
      
      // Resize arrays
      b5_in.resize(invParam.N5);
      c5_in.resize(invParam.N5);
      
      // Copy params
      for(int s =0 ; s < invParam.N5; s++ ) {
	b5_in[s] = toDouble(invParam.b5);
	c5_in[s] = toDouble(invParam.c5);
      }

      // Call Init Funcs
      init(fs_);
    } 
    
    
    /* Destructor */  
    ~LinOpSysSolverMDWFArray() { 
      fini(); state = NULL; 
    }
      
    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<LatticeFermion>& psi, 
				      const multi1d<LatticeFermion>& chi) const;


      
    
    int size() const { return invParam.N5; }
    const Subset& subset() const { return all; }

    
  private:
  //    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Handle<LinearOperatorArray<T> > A;
    SysSolverMDWFParams invParam;
    multi1d<LatticeColorMatrix> u; // The gauge field suitably prepared
    multi1d<double> b5_in;
    multi1d<double> c5_in;

    // Internal Pointers
    QOP_MDWF_State *state;
    QOP_MDWF_Parameters *params;

    /* INIT Function */
    void init(Handle< FermState<T,P,Q> > fermstate);

    /* Finalize -- Cleanup */
    void fini(void);
  };
  

  /*! @} */   // end of group qprop
}


#endif
