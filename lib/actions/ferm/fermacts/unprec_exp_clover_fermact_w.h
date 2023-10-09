// -*- C++ -*-
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_exp_clover_fermact_w_h__
#define __unprec_exp_clover_fermact_w_h__

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "unprec_wilstype_fermact_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace UnprecExpCloverFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Unpreconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Unpreconditioned clover fermion action
   */
  class UnprecExpCloverFermAct
    : public UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix>,
				     multi1d<LatticeColorMatrix>>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! General FermBC
    /*! Isotropic action */
    UnprecExpCloverFermAct(Handle<CreateFermState<T, P, Q>> cfs_, const CloverFermActParams& param_)
      : cfs(cfs_), param(param_)
    {
    }

    //! Copy constructor
    UnprecExpCloverFermAct(const UnprecExpCloverFermAct& a) : cfs(a.cfs), param(a.param)
    {
    }

    //! Produce a linear operator for this action
    UnprecLinearOperator<T, P, Q>* linOp(Handle<FermState<T, P, Q>> state) const;

    LinearOperator<T>* hermitianLinOp(Handle<FermState<T, P, Q>> state) const
    {
      return new lgherm<T>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecExpCloverFermAct()
    {
    }

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T, P, Q>& getCreateState() const
    {
      return *cfs;
    }

    // Hide partial constructor
    UnprecExpCloverFermAct()
    {
    }

    //! Assignment
    void operator=(const UnprecExpCloverFermAct& a)
    {
    }

  private:
    Handle<CreateFermState<T, P, Q>> cfs;
    CloverFermActParams param;
  };

}

#endif
