// -*- C++ -*-
// $Id: unprec_ovdwf_linop_array_w.h,v 1.10 2005-01-02 05:21:10 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) linear operator
 */

#ifndef __unprec_ovdwf_linop_array_w_h__
#define __unprec_ovdwf_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_base_array_w.h"

using namespace QDP;

namespace Chroma
{
  //! Unpreconditioned Overlap-DWF (Borici) linear operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  class UnprecOvDWLinOpArray : public UnprecDWLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecOvDWLinOpArray() {}

    //! Full constructor
    UnprecOvDWLinOpArray(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q, int N5_)
      {create(u_,WilsonMass_,m_q,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, const Real& WilsonMass_, const Real& m_q_, int N5_);

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvDWLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
    void Dminus(LatticeFermion& chi,
		const LatticeFermion& psi,
		enum PlusMinus isign,
		int s5) const;


  private:
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
    WilsonDslash  D;
  };

}; // End Namespace Chroma

using namespace Chroma;

#endif
