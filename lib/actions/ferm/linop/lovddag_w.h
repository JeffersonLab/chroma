// -*- C++ -*-
// $Id: lovddag_w.h,v 1.8 2004-04-16 14:58:29 bjoo Exp $
/*! \file
 *  \brief Internal Overlap-pole operator
 */

#ifndef __lovddag_w_h__
#define __lovddag_w_h__

#include "linearop.h"
#include "fermact.h" 
#include "meas/eig/eig_w.h"
#include <iostream>
using namespace std;

using namespace QDP;

//! Internal Overlap-pole operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Wilson fermions!
 *
 *   Chi  =   (1/2)*((1+m_q) + (1-m_q) * gamma_5 * B) . Psi 
 *  where  B  is the pole approx. to eps(H(m)) 
 *
 * Internally, it computes
 *   Chi  =   ((1+m_q)/(1-m_q) + gamma_5 * B) . Psi 
 * and then rescales at the end to the correct normalization
 *
 *  NOTE: B is hermitian, so       
 *     (1 + gamma_5 * B)^dag = (1 + B * gamma_5) 
 *                           = gamma_5 * (1 + gamma_5 * B) * gamma_5 
 */

class lovddag : public LinearOperator<LatticeFermion>
{
public:
  //! Creation routine
  /*!
   * \ingroup linop
   *
   * \param _MdagM          M^dag.M of underlying linop M      (Read)
   * \param _M              Underlying linop M	               (Read)
   * \param _m_q            quark mass                         (Read)
   * \param _numroot 	    number of poles in expansion       (Read)
   * \param _constP         constant coeff                     (Read)
   * \param _resP           numerator                          (Read)
   * \param _rootQ          denom                              (Read)
   * \param _OperEigVec     eigenvectors      	               (Read)
   * \param _EigValFunc     eigenvalues      	               (Read)
   * \param _NEig           number of eigenvalues              (Read)
   * \param _MaxCG          MaxCG inner CG                     (Read)
   * \param _RsdCG          residual for inner CG              (Read)
   */
  lovddag(const UnprecWilsonTypeFermAct<LatticeFermion>& S_aux,
	   Handle<const ConnectState> state,
	   const Real& _m_q, int _numroot, 
	   const Real& _constP, 
	   const multi1d<Real>& _resP,
	   const multi1d<Real>& _rootQ, 
	   int _NEig,
	   const multi1d<Real>& _EigValFunc,
	   const multi1d<LatticeFermion>& _EigVec,
	   int _MaxCG,
	   const Real& _RsdCG,
	  const int _ReorthFreq,
	   const Chirality _ichiral) :
    M(S_aux.linOp(state)), MdagM(S_aux.lMdagM(state)), m_q(_m_q), 
    numroot(_numroot), constP(_constP),
    resP(_resP), rootQ(_rootQ), EigVec(_EigVec), EigValFunc(_EigValFunc),
    NEig(_NEig), MaxCG(_MaxCG), RsdCG(_RsdCG), ReorthFreq(_ReorthFreq), ichiral(_ichiral) {};

  //! Destructor is automatic
  ~lovddag() {}
 
  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Handle<const LinearOperator<LatticeFermion> > M;
  Handle<const LinearOperator<LatticeFermion> > MdagM;

  // Copy all of these rather than reference them.
  const Real m_q;
  int numroot;
  const Real constP;
  const multi1d<Real> rootQ;
  const multi1d<Real> resP;
  const multi1d<LatticeFermion> EigVec;
  const multi1d<Real> EigValFunc;
  int NEig;
  int MaxCG;
  const Real RsdCG;
  const int  ReorthFreq;
  Chirality ichiral;
};

#endif
