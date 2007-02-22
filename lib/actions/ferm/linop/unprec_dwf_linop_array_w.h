// -*- C++ -*-
// $Id: unprec_dwf_linop_array_w.h,v 3.1 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion linear operator
 */

#ifndef __unprec_dwf_linop_array_w_h__
#define __unprec_dwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/unprec_dwflike_linop_base_array_w.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Unpreconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class UnprecDWLinOpArray : public UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			     multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    UnprecDWLinOpArray(Handle< FermState<T,P,Q> > fs,
		       const Real& WilsonMass_, const Real& m_q, int N5_,
		       const AnisoParam_t& aniso_);

    //! Destructor is automatic
    ~UnprecDWLinOpArray() {}

    //! Length of DW flavor index/space
    inline int size() const {return N5;}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Apply the Dminus operator on a lattice fermion.
    inline
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const
    {
      QDPIO::cerr << "Dminus not implemented" << endl;
      QDP_abort(1);
    }

    //! Derivative
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
	       enum PlusMinus isign) const;


  protected:
    //! Partial constructor
    UnprecDWLinOpArray() {}
    //! Hide =
    void operator=(const UnprecDWLinOpArray&) {}

  private:
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;

    Real fact1;
    Real fact2;

    WilsonDslash  D;
    Handle< FermBC<T,P,Q> > fbc;
  };

} // End Namespace Chroma



#endif
