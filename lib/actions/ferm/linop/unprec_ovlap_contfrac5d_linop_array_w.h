// -*- C++ -*-
// $Id: unprec_ovlap_contfrac5d_linop_array_w.h,v 1.1 2004-09-29 21:48:34 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __unprec_ovlap_contfrac5d_linop_array_w_h__
#define __unprec_ovlap_contfrac5d_linop_array_w_h__

#include "linearop.h"
#include "fermact.h"
#include "state.h"

using namespace QDP;

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

class UnprecOvlapContFrac5DLinOpArray : public LinearOperator< multi1d<LatticeFermion> >
{
public:

  //! Full constructor
  UnprecOvlapContFrac5DLinOpArray(const UnprecWilsonTypeFermAct<LatticeFermion>& S_aux,
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
    m_q(_m_q), N5(_N5), scale_fac(_scale_fac), alpha(_alpha),
    beta(_beta), NEig(_NEig), EigValFunc(_EigValFunc), EigVec(_EigVec),
    M(S_aux.linOp(state)), isLastZeroP(_isLastZeroP) 
  {
    QDPIO::cout << "LinOp isLastZeroP = " << isLastZeroP << endl;
  }

  //! Length of DW flavor index/space
  int size() const {return N5;}

  //! Destructor is automatic
  ~UnprecOvlapContFrac5DLinOpArray() {}

  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (multi1d<LatticeFermion>& chi, 
		   const multi1d<LatticeFermion>& psi, 
		   enum PlusMinus isign) const;

private:
  Handle< const LinearOperator<LatticeFermion> > M;  // Wilson(esque) Op
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

#endif
