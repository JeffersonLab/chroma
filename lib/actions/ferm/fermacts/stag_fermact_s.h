// -*- C++ -*-
// $Id: stag_fermact_s.h,v 1.4 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Staggered fermion action
 */

#ifndef __stag_fermact_w_h__
#define __stag_fermact_w_h__

#include "fermact.h"

using namespace QDP;

namespace Chroma 
{ 
  //! Staggered fermion action
  /*! \ingroup fermact
   *
   */
  class StagFermAct : public EvenOddStaggeredTypeFermAct< LatticeStaggeredFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    StagFermAct(Handle< FermBC<LatticeStaggeredFermion> > fbc_, 
		const Real& Mass_) : 
      fbc(fbc_), Mass(Mass_) {}

    //! Copy constructor
    StagFermAct(const StagFermAct& a) : 
      fbc(a.fbc), Mass(a.Mass) {}

    //! Assignment
    StagFermAct& operator=(const StagFermAct& a)
      {fbc=a.fbc; Mass=a.Mass; return *this;}

    //! Return the quark mass
    const Real getQuarkMass() const {return Mass;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeStaggeredFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddLinearOperator< LatticeStaggeredFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeStaggeredFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~StagFermAct() {}

  private:
    StagFermAct() {} //hide default constructor
  
  private:
    Handle< FermBC<LatticeStaggeredFermion> >  fbc;
    Real Mass;
  };

} // End Namespace Chroma

using namespace Chroma;

#endif
