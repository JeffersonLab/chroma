// -*- C++ -*-
// $Id: asqtad_linop_s.h,v 1.3 2003-12-12 13:56:40 bjoo Exp $
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
 */

class AsqtadLinOp : public EvenOddLinearOperator<LatticeFermion>
{
public:
  //! Partial constructor - Must use create later
  AsqtadLinOp() {}

  //! Full constructor
  AsqtadLinOp(const multi1d<LatticeColorMatrix>& u_fat_, const multi1d<LatticeColorMatrix>& u_triple_, const Real& Mass_) 
  {
    create(u_fat_, u_triple_, Mass_);
  }


  void create(const multi1d<LatticeColorMatrix>& u_fat_, const multi1d<LatticeColorMatrix>& u_triple_, const Real& Mass_) {
    u_fat = u_fat_;
    u_triple = u_triple_;
    Mass = Mass_;

  };
  //! Destructor is automatic
  ~AsqtadLinOp() {}

  //! Apply the the even-even block onto a source vector
  inline void evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		     enum PlusMinus isign) const 
  {
    chi[ rb[0] ] = 2*Mass*psi;
  }
  
  //! Apply the the even-odd block onto a source vector
  void evenOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

  //! Apply the the odd-even block onto a source vector
  void oddEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		    enum PlusMinus isign) const;

  //! Apply the the odd-odd block onto a source vector
  inline void oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
		   enum PlusMinus isign) const
  {
    chi[ rb[1] ] = 2*Mass*psi;
  }

private:
  Real Mass;
  multi1d<LatticeColorMatrix> u_fat;
  multi1d<LatticeColorMatrix> u_triple;
};

#endif
