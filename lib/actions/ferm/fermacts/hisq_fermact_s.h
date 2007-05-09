// -*- C++ -*-
// $Id: hisq_fermact_s.h,v 1.1 2007-05-09 12:43:20 mcneile Exp $
/*! \file
 *  \brief Hisq staggered fermion action
 */

#ifndef __hisq_fermact_s_h__
#define __hisq_fermact_s_h__

#include "stagtype_fermact_s.h"
#include "state.h"
#include "actions/ferm/fermstates/asqtad_state.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermacts/hisq_fermact_params_s.h"


namespace Chroma 
{ 
  //! Name and registration
  namespace HisqFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Hisq staggered fermion action
  /*! \ingroup fermacts
   *
   */
  class HisqFermAct : public EvenOddStaggeredTypeFermAct<
    LatticeStaggeredFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General CreateFermState
    HisqFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		  const HisqFermActParams& p) :
      cfs(cfs_), param(p) {}
  
    //! Copy constructor
    HisqFermAct(const HisqFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Create state should apply the BC
    AsqtadConnectStateBase* createState(const Q& u_) const;

    //! Return the fermion BC object for this action
    const FermBC<T,P,Q>& getFermBC() const {return cfs->getBC();}

    //! Produce a linear operator for this action
    EvenOddLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state_) const;

    //! Produce a linear operator M^dag.M for this action
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state_) const;

    //! Return quark prop solver, solution of unpreconditioned system
    SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam) const;

    //! accessors 
    const Real getQuarkMass() const {return param.Mass;}
    Real getU0() {return param.u0;}

    //! Destructor is automatic
    ~HisqFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    HisqFermAct() {} //hide default constructor
    void operator=(const HisqFermAct& a) {} // Assignment
  
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    HisqFermActParams  param;
  };


}; // End Namespace Chroma

#endif
