// -*- C++ -*-
// $Id: unprec_w12_fermact_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned W12 fermion action
 */

#ifndef __unprec_w12_fermact_w_h__
#define __unprec_w12_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecW12FermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Unpreconditioned W12 fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecW12FermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
			   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecW12FermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		     const CloverFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecW12FermAct(const UnprecW12FermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    //! Destructor is automatic
    ~UnprecW12FermAct() {}

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecW12FermAct() {} //hide default constructor
    void operator=(const UnprecW12FermAct& a) {} // Hide =
  
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverFermActParams param;
  };

}

#endif
