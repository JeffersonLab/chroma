// -*- C++ -*-
// $Id: stag_fermact_s.h,v 3.1 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Staggered fermion action
 */

#ifndef __stag_fermact_w_h__
#define __stag_fermact_w_h__

#include "stagtype_fermact_s.h"


namespace Chroma 
{ 
  //! Staggered fermion action
  /*! \ingroup fermacts
   *
   */
  class StagFermAct : public EvenOddStaggeredTypeFermAct<LatticeStaggeredFermion, 
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    StagFermAct(Handle< FermBC<T,P,Q> > fbc_, 
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

    //! Produce a linear operator for this action
    EvenOddLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const;

    //! Destructor is automatic
    ~StagFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    StagFermAct() {} //hide default constructor
  
  private:
    Handle< FermBC<T,P,Q> >  fbc;
    Real Mass;
  };

} // End Namespace Chroma


#endif
