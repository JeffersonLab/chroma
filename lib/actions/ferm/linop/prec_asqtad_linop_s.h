// -*- C++ -*-
// $Id: prec_asqtad_linop_s.h,v 1.2 2003-12-10 16:21:00 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */
// NEW $Id: asqtad_linop_s.h 2003/11/13 steve 
// Asqtad Staggered fermion linear operator

#ifndef __asqtad_linop_s_h__
#define __asqtad_linop_s_h__

#include "linearop.h"

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

class EvenOddPrecAsqtadLinOp : public EvenOddPrecLinearOperator<LatticeFermion>
{
public:
  //! Partial constructor - Must use create later
  EvenOddPrecAsqtadLinOp() {}

  //! Full constructor
  EvenOddPrecAsqtadLinOp(const multi1d<LatticeColorMatrix>& u_fat_, const multi1d<LatticeColorMatrix>& u_triple_, const Real& Mass_) 
  {
    create(u_fat_, u_triple_, Mass_);
  }


  void create(const multi1d<LatticeColorMatrix>& u_fat_, const multi1d<LatticeColorMatrix>& u_triple_, const Real& Mass_) {
    u_fat = u_fat_;
    u_triple = u_triple_;
    Mass = Mass_;

  };
  //! Destructor is automatic
  ~EvenOddPrecAsqtadLinOp() {}

  //! Apply the the even-even block onto a source vector
  void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const;
  
  //! Apply the the even-odd block onto a source vector
  void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

  //! Apply the the odd-even block onto a source vector
  void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

  //! Apply the the odd-odd block onto a source vector
  void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		   enum PlusMinus isign) const;

  void evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
                   enum PlusMinus isign) const { };

private:
  Real Mass;
  multi1d<LatticeColorMatrix> u_fat;
  multi1d<LatticeColorMatrix> u_triple;
};

#endif
