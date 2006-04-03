// -*- C++ -*-
// $Id: asqtad_fermact_s.h,v 3.0 2006-04-03 04:58:44 edwards Exp $
/*! \file
 *  \brief Asqtad staggered fermion action
 */

#ifndef __asqtad_fermact_s_h__
#define __asqtad_fermact_s_h__

#include "fermact.h"
#include "state.h"
#include "actions/ferm/fermacts/asqtad_state.h"
#include "actions/ferm/fermacts/simple_fermstate.h"


namespace Chroma 
{ 
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
    AsqtadFermAct(Handle< FermBC<T,P,Q> > fbc_, 
		  const Real Mass_, const Real u0_) : 
      cfs(new CreateSimpleFermState<T,P,Q>(fbc_)), Mass(Mass_), u0(u0_) {}
  
    //! Copy constructor
    AsqtadFermAct(const AsqtadFermAct& a) : 
      cfs(a.cfs), Mass(a.Mass), u0(a.u0) {}

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
			   const InvertParam_t& invParam) const;

    //! accessors 
    const Real getQuarkMass() const {return Mass;}
    Real getU0() {return u0;}

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
    Real Mass;
    Real u0;
  };


}; // End Namespace Chroma

#endif
