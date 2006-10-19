// -*- C++ -*-
// $Id: eoprec_wilson_fermact_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecWilsonFermActEnv
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
  class EvenOddPrecWilsonFermAct : public EvenOddPrecConstDetWilsonTypeFermAct<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const Real& Mass_) : 
      cfs(cfs_) {param.Mass=Mass_;}

    //! General FermBC with Anisotropy
    EvenOddPrecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const WilsonFermActParams& param_) :
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    EvenOddPrecWilsonFermAct(const EvenOddPrecWilsonFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

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
    ~EvenOddPrecWilsonFermAct() {}

    Double getQuarkMass(void) const { 
      return param.Mass;
    }

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    EvenOddPrecWilsonFermAct() {}
    //! Assignment
    void operator=(const EvenOddPrecWilsonFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    WilsonFermActParams param;
  };

}

#endif
