// -*- C++ -*-
// $Id: unprec_clover_fermact_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_clover_fermact_w_h__
#define __unprec_clover_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace UnprecCloverFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Unpreconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Unpreconditioned clover fermion action
   */
  class UnprecCloverFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
			      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    /*! Isotropic action */
    UnprecCloverFermAct(Handle< CreateFermState<T,P,Q> > cfs_,
			const CloverFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecCloverFermAct(const UnprecCloverFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    //! Destructor is automatic
    ~UnprecCloverFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    // Hide partial constructor
    UnprecCloverFermAct() {}

    //! Assignment
    void operator=(const UnprecCloverFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> > cfs;
    CloverFermActParams param;
  };

}

#endif
