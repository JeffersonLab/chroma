// -*- C++ -*-
// $Id: unprec_ovdwf_linop_array_w.h,v 3.1 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) linear operator
 */

#ifndef __unprec_ovdwf_linop_array_w_h__
#define __unprec_ovdwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/unprec_dwflike_linop_base_array_w.h"


namespace Chroma
{
  //! Unpreconditioned Overlap-DWF (Borici) linear operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class UnprecOvDWLinOpArray : public UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			       multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    UnprecOvDWLinOpArray(Handle< FermState<T,P,Q> > fs, 
			 const Real& WilsonMass_, const Real& m_q, int N5_)
      {create(fs,WilsonMass_,m_q,N5_);}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs, 
		const Real& WilsonMass_, const Real& m_q_, int N5_);

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvDWLinOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;

  protected:
    //! Partial constructor
    UnprecOvDWLinOpArray() {}
    //! Hide =
    void operator=(const UnprecOvDWLinOpArray&) {}

  private:
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
    WilsonDslash  D;
    Handle< FermBC<T,P,Q> > fbc;
  };

} // End Namespace Chroma


#endif
