// -*- C++ -*-
// $Id: prec_asqtad_linop_s.h,v 1.1 2003-12-10 14:27:08 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */
// NEW $Id: asqtad_linop_s.h 2003/11/13 steve 
// Asqtad Staggered fermion linear operator

#ifndef __asqtad_linop_s_h__
#define __asqtad_linop_s_h__

#include "linearop.h"
// #include "actions/ferm/linop/dslash_s.h"
#include "actions/ferm/dslash/dslash_s.h"

using namespace QDP;

//! Asqtad Staggered-Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Staggered fermions!
 *
 *                                     ~      ~+
 * This subroutine applies the matrix  M  or  M   to the vector
 * Psi,
 *
 *      	       	   {   ~
 *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
 *      	   Chi  =  {
 *      	       	   {   ~   +
 *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS

 * Algorithm:

 * The kernel for Staggered fermions is

 *      M  =  2m + D'
 *
 *  The KS phases are already included in the fat and triple links.
 *  This is done when creating the "connectState" where the fat
 *  and triple links are themselves computed.
 *
 */

class AsqtadLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  AsqtadLinOp() {}

  //! Full constructor
  AsqtadLinOp(const multi1d<LatticeColorMatrix>& u_fat_, const multi1d<LatticeColorMatrix>& u_triple_, const Real& Mass_) : u_fat(u_fat_), u_triple(u_triple_), Mass(Mass_) {};

  //! Destructor is automatic
  ~AsqtadLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Real Mass;
  multi1d<LatticeColorMatrix>& u_fat;
  multi1d<LatticeColorMatrix>& u_triple;
  //AsqtadDslash D;
};

#endif
