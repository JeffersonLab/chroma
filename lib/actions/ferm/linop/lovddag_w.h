// -*- C++ -*-
// $Id: lovddag_w.h,v 1.1 2003-04-09 05:57:15 edwards Exp $
/*! \file
 *  \brief Internal Overlap-pole operator for D^dag.D
 */

#ifndef __lovddag_w_h__
#define __lovddag_w_h__

#include "linearop.h"

using namespace QDP;

//! Apply overlap operator to a chiral source
/*!
 * \ingroup linop
 *
 *   Chi  =   (1/4)*(2*(1+m_q^2) + (1-m_q^2)*(gamma_5 +- 1) * B ) . Psi_+- 
 *
 *  where  B  is the pole approx. to eps(H(m)) 
 *
 * Internally, it computes
 *   Chi  =   (2*(1+m_q^2)/(1-m_q^2) + (gamma_5 +- 1) * B ) . Psi_+- 
 * and then rescales at the end to the correct normalization
 *
 *  NOTE: B is hermitian, so       
 *     (1 + gamma_5 * B)^dag = (1 + B * gamma_5) 
 *                           = gamma_5 * (1 + gamma_5 * B) * gamma_5
 *
 */
class lovddag : public LinearOperator
{
public:
  //! Creation routine
  /*!
   * \ingroup linop
   *
   * \param _A              A lovlapms linear operator          (Read)
   * \param _ichiral        Chirality of expected source vector (Read)
   */
  lovddag(const lovlapms& _A, int _ichiral) :
    A(_A), ichiral(_ichiral) {}

  //! Destructor is automatic
  ~lovddag() {}

  //! Only defined on the entire lattice
  const Subset& subset() const {return all;}

  //! Apply the operator onto a chiral source vector
  LatticeFermion operator() (const LatticeFermion& psi, enum LinOpSign isign) const;

private:
  const lovlapms& A;
  const int ichiral;
};

#endif
