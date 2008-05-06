// -*- C++ -*-
// $Id: mdwf_solver.h,v 1.4 2008-05-06 13:25:28 bjoo Exp $
/*! \file
 *  \brief DWF/SSE double-prec solver
 */

#ifndef MDWF_SOLVER_H
#define MDWF_SOLVER_H


extern "C" {
  struct QOP_MDWF_State;
  struct QOP_MDWF_Parameters;  
};

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


using namespace QDP;
namespace Chroma 
{ 
  //! AVP's DWF Solver interface
  /*!
   * \ingroup qprop
   *
   * @{
   */
  class MDWFQpropT : public SystemSolverArray<LatticeFermion>  { 
    
  public:
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    /* Constructor */
    MDWFQpropT(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
	       Handle< LinOpSystemSolverArray<T> > invA_,   // throw away
	       Handle< FermState<T,P,Q> > state_, 
	       const Real& OverMass_,
	       const Real& Mass_,
	       const AnisoParam_t& anisoParam_,
	       const GroupXML_t& invParam_) : A(A_), 
					      OverMass(OverMass_), 
					      Mass(Mass_),       
					      N5(A->size()), 
					      anisoParam(anisoParam_)  {
           init(state_, invParam_);
    } 
    
      
    /* Destructor */  
    ~MDWFQpropT() { 
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


      
    
    int size() const { return N5; }
    const Subset& subset() const { return all; }

    
  private:
    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Real OverMass;              // M5
    Real Mass;                  // m_f
    int N5;                     // The 5th Dimension
    AnisoParam_t anisoParam;    // Anisotropy
    SysSolverCGParams invParam; // Inverter Parameters
    multi1d<LatticeColorMatrix> u; // The gauge field suitably prepared
    
    // Internal Pointers
    QOP_MDWF_State *state;
    QOP_MDWF_Parameters *params;

    /* INIT Function */
    void init(Handle< FermState<T,P,Q> > fermstate, const GroupXML_t& inv);

    /* Finalize -- Cleanup */
    void fini(void);
  };
  

  /*! @} */   // end of group qprop
}


#endif
