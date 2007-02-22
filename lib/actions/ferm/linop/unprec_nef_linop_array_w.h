// -*- C++ -*-
// $Id: unprec_nef_linop_array_w.h,v 3.1 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion linear operator
 */

#ifndef __unprec_nef_linop_array_w_h__
#define __unprec_nef_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/unprec_dwflike_linop_base_array_w.h"


namespace Chroma
{
  //! Unpreconditioned domain-wall Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class UnprecNEFDWLinOpArray : public UnprecDWLikeLinOpBaseArray<LatticeFermion,
				multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    /*!
      Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
      Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
    */
    UnprecNEFDWLinOpArray(Handle< FermState<T,P,Q> > fs,
			  const Real& WilsonMass_, 
			  const multi1d<Real>& b5_, const multi1d<Real>& c5_, 
			  const Real& m_q_, int N5_)
      {create(fs,WilsonMass_,b5_,c5_,m_q_,N5_);}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& WilsonMass_, 
		const multi1d<Real>& b5_, const multi1d<Real>& c5_,
		const Real& m_q_, int N5_);

    //! Length of DW flavor index/space
    inline int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecNEFDWLinOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;

    //! Derivative
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
	       enum PlusMinus isign) const;

  protected:
    //! Partial constructor
    UnprecNEFDWLinOpArray() {}
    //! Hide =
    void operator=(const UnprecNEFDWLinOpArray&) {}

  private:
    Real WilsonMass;
    multi1d<Real> b5;
    multi1d<Real> c5;
    Real m_q;
    int  N5;
    WilsonDslash  D;
    Handle< FermBC<T,P,Q> > fbc;

    multi1d<Real> fb5;
    multi1d<Real> fc5;
  };

} // End Namespace Chroma



#endif
