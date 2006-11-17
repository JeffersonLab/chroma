// -*- C++ -*-
// $Id: asqtad_fermact_s.h,v 3.4 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief Asqtad staggered fermion action
 */

#ifndef __asqtad_fermact_s_h__
#define __asqtad_fermact_s_h__

#include "stagtype_fermact_s.h"
#include "state.h"
#include "actions/ferm/fermstates/asqtad_state.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermacts/asqtad_fermact_params_s.h"


namespace Chroma 
{ 
  //! Name and registration
  namespace AsqtadFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Asqtad staggered fermion action
  /*! \ingroup fermacts
   *
   */
  class AsqtadFermAct : public EvenOddStaggeredTypeFermAct<
    LatticeStaggeredFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General CreateFermState
    AsqtadFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		  const AsqtadFermActParams& p) :
      cfs(cfs_), param(p) {}
  
    //! Copy constructor
    AsqtadFermAct(const AsqtadFermAct& a) : 
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
    ~AsqtadFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    AsqtadFermAct() {} //hide default constructor
    void operator=(const AsqtadFermAct& a) {} // Assignment
  
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    AsqtadFermActParams  param;
  };


}; // End Namespace Chroma

#endif
