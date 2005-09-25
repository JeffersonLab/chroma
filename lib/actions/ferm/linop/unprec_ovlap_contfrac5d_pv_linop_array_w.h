// -*- C++ -*-
// $Id: unprec_ovlap_contfrac5d_pv_linop_array_w.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Pauli-Villars Continued Fraction 5D
 */

#ifndef __unprec_ovlap_contfrac5d_pv_linop_array_w_h__
#define __unprec_ovlap_contfrac5d_pv_linop_array_w_h__

#include "linearop.h"
#include "fermact.h"
#include "state.h"


namespace Chroma 
{ 
  //! Unpreconditioned Pauli-Villars Continued Fraction 5D
  /*!
   * \ingroup linop
   *
   * Paulli-Villars cont. frac.
   */

  class UnprecOvlapContFrac5DPVLinOpArray : public UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:

    //! Full constructor
    UnprecOvlapContFrac5DPVLinOpArray(const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& S_aux,
				      Handle<const ConnectState> state,
				      const Real& _m_q,
				      int _N5,
				      const Real& _scale_fac,
				      const multi1d<Real>& _alpha,
				      const multi1d<Real>& _beta,
				      int _NEig,
				      const multi1d<Real>& _EigValFunc,
				      const multi1d<LatticeFermion>& _EigVec,
				      const bool _isLastZeroP ) :
      M(S_aux.linOp(state)), m_q(_m_q), N5(_N5), scale_fac(_scale_fac), alpha(_alpha),
      beta(_beta), NEig(_NEig), EigVec(_EigVec), EigValFunc(_EigValFunc),
      isLastZeroP(_isLastZeroP) 
    {
      QDPIO::cout << "LinOpPV isLastZeroP = " << isLastZeroP << endl;
    }

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvlapContFrac5DPVLinOpArray() {}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

    //! Derivative
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
	       enum PlusMinus isign) const;

  private:
    Handle< const DiffLinearOperator<LatticeFermion, multi1d<LatticeColorMatrix> > > M;  // Wilson(esque) Op
    const Real m_q;
    const int  N5;    // Size of the 5th dimension
    const Real scale_fac;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const int NEig;
    const multi1d<LatticeFermion> EigVec;
    const multi1d<Real> EigValFunc;
    const bool isLastZeroP;
  };

}; // End Namespace Chroma


#endif
