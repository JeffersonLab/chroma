// -*- C++ -*-
// $Id: overlapbu_linop_w.h,v 2.0 2005-09-25 21:04:29 edwards Exp $
/*! \file
 *  \brief A variant of the 4D overlap operator
 */

#ifndef __overlapbu_linop_w_h__
#define __overlapbu_linop_w_h__

#include "linearop.h"


namespace Chroma 
{ 
  //! Unpreconditioned Wilson-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *                                                      ~      ~+
   * This subroutine applies the unpreconditioned matrix  M  or  M   the vector
   * Psi,
   *
   *      	       	   {   ~
   *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
   *      	   Chi  =  {
   *      	       	   {   ~   +
   *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS

   * Algorithm:

   * The kernel for Wilson fermions is

   *      M  =  1 - k D'
   *
   */

  class OverlapBULinOp : public LinearOperator<LatticeFermion>
  {
  public:
    //! Partial constructor
    OverlapBULinOp() {}

    //! Full constructor
    OverlapBULinOp(const multi1d<LatticeColorMatrix>& _u, const Real& _OverMass, const Real& _m_q)
      {create(_u,_OverMass,_m_q);}

    //! Destructor is automatic
    ~OverlapBULinOp() {}

    //! Only defined on the odd subset
    const OrderedSubset& subset() const {return all;}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& _u, const Real& _OverMass, const Real& _m_q);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

  private:
    Real OverMass;
    Real m_q;
    multi1d<LatticeColorMatrix> u;
  };


}; // End Namespace Chroma


#endif
