// -*- C++ -*-
// $Id: unprec_w12_fermact_w.h,v 2.1 2006-02-09 02:23:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned W12 fermion action
 */

#ifndef __unprec_w12_fermact_w_h__
#define __unprec_w12_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecW12FermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Unpreconditioned W12 fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecW12FermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecW12FermAct(Handle< FermBC<LatticeFermion> > fbc_, 
		     const CloverFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecW12FermAct(const UnprecW12FermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecW12FermAct& operator=(const UnprecW12FermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecW12FermAct() {}

  private:
    UnprecW12FermAct() {} //hide default constructor
   
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    CloverFermActParams param;
  };

}

#endif
