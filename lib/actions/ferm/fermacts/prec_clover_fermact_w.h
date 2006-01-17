// -*- C++ -*-
// $Id: prec_clover_fermact_w.h,v 2.4 2006-01-17 16:01:46 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_clover_fermact_w_h__
#define __prec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

namespace Chroma 
{ 
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecCloverFermActEnv
  {
    extern const std::string name;
    extern const bool registered;   /*!< Name to be used */
  }
  

  //! Even-odd preconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned clover fermion action. 
   * Only defined on odd subset.
   */

  class EvenOddPrecCloverFermAct : public EvenOddPrecLogDetWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const CloverFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    EvenOddPrecCloverFermAct(const EvenOddPrecCloverFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    EvenOddPrecCloverFermAct& operator=(const EvenOddPrecCloverFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddPrecLogDetLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~EvenOddPrecCloverFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    CloverFermActParams param;
  };

}; // End Namespace Chroma


#endif
