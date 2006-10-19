// -*- C++ -*-
// $Id: tprec_wilson_fermact_w.h,v 3.1 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Time-preconditioned (no even-odd) Wilson fermion action
 */

#ifndef __tprec_wilson_fermact_w_h__
#define __tprec_wilson_fermact_w_h__

#include "tprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace TimePrecWilsonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Even-odd preconditioned Wilson fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned wilson fermion action. 
   * Only defined on odd subset.
   */
  class TimePrecWilsonFermAct : public TimePrecConstDetWilsonTypeFermAct<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    TimePrecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const Real& Mass_) : 
      cfs(cfs_) {param.Mass=Mass_;}

    //! General FermBC with Anisotropy
    TimePrecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const WilsonFermActParams& param_) :
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    TimePrecWilsonFermAct(const TimePrecWilsonFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    TimePrecConstDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<LatticeFermion>(linOp(state));
      }

    //! Return a linear operator solver for this action to solve M*psi=chi 
    LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~TimePrecWilsonFermAct() {}

    Double getQuarkMass(void) const { 
      return param.Mass;
    }

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    TimePrecWilsonFermAct() {}
    //! Assignment
    void operator=(const TimePrecWilsonFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    WilsonFermActParams param;
  };

}

#endif
