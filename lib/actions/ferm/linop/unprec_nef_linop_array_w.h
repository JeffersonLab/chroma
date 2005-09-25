// -*- C++ -*-
// $Id: unprec_nef_linop_array_w.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
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
  class UnprecNEFDWLinOpArray : public UnprecDWLikeLinOpBaseArray<LatticeFermion,multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecNEFDWLinOpArray() {}

    //! Full constructor
    /*!
      Set b5 = 1.0 and c5=0.0 to get Shamir DWF with a5=1.
      Set b5 = 1.0 and c5=1.0 to get Borichi DWF.
    */
    UnprecNEFDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
			  const Real& WilsonMass_, 
			  const multi1d<Real>& b5_, const multi1d<Real>& c5_, 
			  const Real& m_q_, int N5_)
      {create(u_,WilsonMass_,b5_,c5_,m_q_,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& WilsonMass_, 
		const multi1d<Real>& b5_, const multi1d<Real>& c5_,
		const Real& m_q_, int N5_);

    //! Length of DW flavor index/space
    inline int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecNEFDWLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

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

  private:
    Real WilsonMass;
    multi1d<Real> b5;
    multi1d<Real> c5;
    Real m_q;
    int  N5;
    WilsonDslash  D;

    multi1d<Real> fb5;
    multi1d<Real> fc5;
  };

}; // End Namespace Chroma



#endif
