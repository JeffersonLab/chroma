// -*- C++ -*-
// $Id: unprec_ovext_linop_array_w.h,v 1.7 2005-04-21 14:04:38 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __unprec_ovext_linop_array_w_h__
#define __unprec_ovext_linop_array_w_h__

#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP;

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
			  const int Npoles_,
			  const Real& coeffP_,
			  const multi1d<Real>& resP_,
			  const multi1d<Real>& rootQ_,
			  const multi1d<Real>& beta_,
			  const Real& OverMass_,
			  const Real& Mass_,
			  const Real& b5_,
			  const Real& c5_)

    {create(u_,Npoles_, coeffP_, resP_, rootQ_, beta_, 
	    OverMass_,Mass_,b5_,c5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		const int Npoles_,
		const Real& coeffP_,
		const multi1d<Real>& resP_,
		const multi1d<Real>& rootQ_,
		const multi1d<Real>& beta_,
		const Real& OverMass_,
		const Real& m_q_,
		const Real& b5_,
		const Real& c5_);

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvExtLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Derivative
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
	       enum PlusMinus isign) const;

  private:
    Handle< const UnprecWilsonLinOp > Dw;
    int Npoles;
    int N5;
    Real R;
    Real alpha;
    Real a5;
    Real coeffP;
    multi1d<Real> p_by_beta_sqrt;
    multi1d<Real> q_sqrt;
    multi1d<Real> beta;
  };

}; // End Namespace Chroma


#endif
