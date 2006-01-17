// -*- C++ -*-
// $Id: prec_wilson_fermact_w.h,v 2.4 2006-01-17 16:01:46 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"


namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecWilsonFermActEnv
  {
    extern const std::string name;
    extern const bool registered;   /*!< Name to be used */
  }
  

  //! Even-odd preconditioned Wilson fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned wilson fermion action. 
   * Only defined on odd subset.
   */
  class EvenOddPrecWilsonFermAct : public EvenOddPrecConstDetWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const Real& Mass_) : 
      fbc(fbc_) {param.Mass=Mass_;}

    //! General FermBC with Anisotropy
    EvenOddPrecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const WilsonFermActParams& param_) :
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    EvenOddPrecWilsonFermAct(const EvenOddPrecWilsonFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    EvenOddPrecWilsonFermAct& operator=(const EvenOddPrecWilsonFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddPrecConstDetLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~EvenOddPrecWilsonFermAct() {}

    Double getQuarkMass(void) const { 
      return param.Mass;
    }

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    WilsonFermActParams param;
  };

}

#endif
