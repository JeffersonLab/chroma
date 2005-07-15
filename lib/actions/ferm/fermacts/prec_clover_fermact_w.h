// -*- C++ -*-
// $Id: prec_clover_fermact_w.h,v 1.9 2005-07-15 11:06:10 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_clover_fermact_w_h__
#define __prec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"


namespace Chroma 
{ 
  //! Even-odd preconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned clover fermion action. 
   * Only defined on odd subset.
   */

  class EvenOddPrecCloverFermAct : public EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const Real& Mass_, const Real& ClovCoeff_, const Real& u0_) : 
      fbc(fbc_), Mass(Mass_), ClovCoeff(ClovCoeff_), u0(u0_) {}

    //! Copy constructor
    EvenOddPrecCloverFermAct(const EvenOddPrecCloverFermAct& a) : 
      fbc(a.fbc), Mass(a.Mass), ClovCoeff(a.ClovCoeff), u0(a.u0) {}

    //! Assignment
    EvenOddPrecCloverFermAct& operator=(const EvenOddPrecCloverFermAct& a)
    {fbc=a.fbc; Mass=a.Mass; ClovCoeff=a.ClovCoeff; u0=a.u0; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddPrecLinearOperator<LatticeFermion> linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion> lMdagM(Handle<const ConnectState> state) const;

    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~EvenOddPrecCloverFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    Real Mass;
    Real ClovCoeff;
    Real u0;
  };

}; // End Namespace Chroma


#endif
