// -*- C++ -*-
// $Id: unprec_ovext_linop_array_w.h,v 1.5 2005-01-14 20:13:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __unprec_ovext_linop_array_w_h__
#define __unprec_ovext_linop_array_w_h__

#include "linearop.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * This operator implements  hep-lat/0005004
   */

  class UnprecOvExtLinOpArray : public UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    UnprecOvExtLinOpArray() {}

    //! Full constructor
    UnprecOvExtLinOpArray(const multi1d<LatticeColorMatrix>& u_, 
			  const Real& WilsonMass_, const Real& m_q, int N5_)
    {create(u_,WilsonMass_,m_q,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const Real& WilsonMass_, const Real& m_q_, int N5_);

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvExtLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

  private:
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;    // total number of 4D fields
    int  NN5;   // number of poles

    multi1d<Real> cc;
    multi1d<Real> ss;
    Real          fact1;
    Real          fact2;

    UnprecWilsonLinOp  D;
  };

}; // End Namespace Chroma


#endif
