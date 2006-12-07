// -*- C++ -*-
// $Id: klein_gordon_linop_s.h,v 1.1 2006-12-07 18:26:18 edwards Exp $
//! Klein-Gordon operator
/*! \file
 *  \ingroup linop
 *
 *  \brief Klein-Gordon boson action masquerading action as a staggered action
 */

#ifndef __klein_gordon_linop_s_h__
#define __klein_gordon_linop_s_h__

#include "linearop.h"
#include "state.h"
#include "io/aniso_io.h"

namespace Chroma 
{ 
  class KleinGordonLinOp : public UnprecLinearOperator< LatticeStaggeredFermion, 
			   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor - Must use create later
    KleinGordonLinOp() {}

    //! Full constructor
    KleinGordonLinOp(Handle< FermState<T,P,Q> > state_,
		     const Real& Mass_,
		     const AnisoParam_t& aniso);

    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const AnisoParam_t& aniso);

    //! Destructor is automatic
    ~KleinGordonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Apply the operator onto a source vector
    void operator() (T& chi, const T& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned operator
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const T& chi, const T& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    Real                          fact;  // tmp holding  Nd+Mass
    Real                          Mass;
    AnisoParam_t                  anisoParam;
    Handle< FermBC<T,P,Q> >       fbc;
    multi1d<LatticeColorMatrix>   u; 
 };


} // End Namespace Chroma


#endif
