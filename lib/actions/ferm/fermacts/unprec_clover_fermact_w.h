// -*- C++ -*-
// $Id: unprec_clover_fermact_w.h,v 2.2 2005-12-29 05:36:42 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_clover_fermact_w_h__
#define __unprec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace UnprecCloverFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Unpreconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Unpreconditioned clover fermion action
   */
  class UnprecCloverFermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    /*! Isotropic action */
    UnprecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			const CloverFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecCloverFermAct(const UnprecCloverFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecCloverFermAct& operator=(const UnprecCloverFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecCloverFermAct() {}

  private:
    // Hide partial constructor
    UnprecCloverFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    CloverFermActParams param;
  };

}

#endif
