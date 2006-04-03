// -*- C++ -*-
// $Id: asqtad_linop_s.h,v 3.0 2006-04-03 04:58:49 edwards Exp $
//! Asqtad Staggered-Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Staggered fermions!
 */

#ifndef __asqtad_linop_s_h__
#define __asqtad_linop_s_h__

#include "linearop.h"
#include "actions/ferm/linop/asqtad_dslash.h"


namespace Chroma 
{ 
  class AsqtadLinOp : public EvenOddLinearOperator< LatticeStaggeredFermion, 
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor - Must use create later
    AsqtadLinOp() {}

    //! Full constructor
    AsqtadLinOp(Handle<AsqtadConnectStateBase> state_, const Real& Mass_) 
    {
      create(state_, Mass_);
    }

    //! Creation routine
    void create(Handle<AsqtadConnectStateBase> state_, const Real& Mass_) 
    {
      //u_fat = u_fat_;
      // u_triple = u_triple_;
      Mass = Mass_;
      D.create(state_);
    };

    //! Destructor is automatic
    ~AsqtadLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<LatticeStaggeredFermion,
		 multi1d<LatticeColorMatrix>,
		 multi1d<LatticeColorMatrix> >& getFermBC() const {return D.getFermBC();}

    //! Apply the the even-even block onto a source vector
    inline void evenEvenLinOp(LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, 
			      enum PlusMinus isign) const 
    {
      chi[ rb[0] ] = 2*Mass*psi;
    }
  
    //! Apply the the even-odd block onto a source vector
    void evenOddLinOp(LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-even block onto a source vector
    void oddEvenLinOp(LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, 
		      enum PlusMinus isign) const;

    //! Apply the the odd-odd block onto a source vector
    inline void oddOddLinOp(LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, 
			    enum PlusMinus isign) const
    {
      chi[ rb[1] ] = 2*Mass*psi;
    }

  private:
    Real Mass;

    // These are really only needed for D. I bring back D here where
    // Steve originally had it.  We don't need u_fat and u_triple here
    // they are kept in the action now.
    // multi1d<LatticeColorMatrix> u_fat;
    // multi1d<LatticeColorMatrix> u_triple;
    AsqtadDslash D;
  };


} // End Namespace Chroma


#endif
