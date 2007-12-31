// -*- C++ -*-
// $Id: unprec_graphene_fermact_w.h,v 3.1 2007-12-31 23:24:26 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Graphene fermion action.
 *
 * This formulation follows Borici's variant of Creutz's graphene
 * fermion construction. Borici's variant is described in
 * arXiv:0712.4401 and Cruetz's original construction is described
 * in arXiv:0712.1201
 */

#ifndef __unprec_graphene_fermact_w_h__
#define __unprec_graphene_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecGrapheneFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Unpreconditioned Graphene fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecGrapheneFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
			      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecGrapheneFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			  const Real& Mass_) : 
      cfs(cfs_) {param.Mass=Mass_;}

    //! General FermBC
    UnprecGrapheneFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			  const WilsonFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecGrapheneFermAct(const UnprecGrapheneFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Assignment
    UnprecGrapheneFermAct& operator=(const UnprecGrapheneFermAct& a)
      {cfs=a.cfs; param=a.param; return *this;}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    //! Destructor is automatic
    ~UnprecGrapheneFermAct() {}

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecGrapheneFermAct() {} //hide default constructor
   
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    WilsonFermActParams param;
  };

}

#endif
