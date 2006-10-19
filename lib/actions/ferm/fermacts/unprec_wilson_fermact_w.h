// -*- C++ -*-
// $Id: unprec_wilson_fermact_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_wilson_fermact_w_h__
#define __unprec_wilson_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecWilsonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Unpreconditioned Wilson fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecWilsonFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
			      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			const Real& Mass_) : 
      cfs(cfs_) {param.Mass=Mass_;}

    //! General FermBC
    UnprecWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			const WilsonFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecWilsonFermAct(const UnprecWilsonFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Assignment
    UnprecWilsonFermAct& operator=(const UnprecWilsonFermAct& a)
      {cfs=a.cfs; param=a.param; return *this;}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    //! Destructor is automatic
    ~UnprecWilsonFermAct() {}

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecWilsonFermAct() {} //hide default constructor
   
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    WilsonFermActParams param;
  };

}

#endif
