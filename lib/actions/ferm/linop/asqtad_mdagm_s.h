// -*- C++ -*-
// $Id: asqtad_mdagm_s.h,v 2.0 2005-09-25 21:04:28 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */
// NEW $Id: asqtad_linop_s.h 2003/11/13 steve 
// Asqtad Staggered fermion linear operator

#ifndef __asqtad_mdagm_s_h__
#define __asqtad_mdagm_s_h__

#include "linearop.h"
#include "actions/ferm/linop/asqtad_dslash.h"



namespace Chroma 
{ 
//! Asqtad Staggered-Dirac operator
/*!
 * \ingroup linop
 *
 * This routine is specific to Staggered fermions!
 *
 *                                           +
 * This subroutine applies the matrix  or  (M M)    to the vector
 *      					E,E
 * Psi,
 *
 *      	       	   {   ~
 *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
 *      	   Chi  =  {
 *      	       	   {   ~   +
 *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS
 * NEED TO THINK MORE ABOUT THE ISIGN HERE AS IS IT NOT REALLY NEEDED!!!
 * FOR NOW JUST CALL THIS ROUTINE WITH A PLUS!!

 * Algorithm:

 * The kernel for Staggered fermions is
 *        +
 *      (M M) =  4m**2  - D  D
 * 	     E             EO OE
 */

class AsqtadMdagM : public LinearOperator<LatticeStaggeredFermion>
{
public:
  //! Partial constructor
  AsqtadMdagM() {}

  //! Full constructor
  AsqtadMdagM(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple, const Real& _Mass)
    {create(_u_fat,_u_triple,_Mass);}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& _u_fat,const multi1d<LatticeColorMatrix>& _u_triple, const Real& _Mass);
    
  //! Destructor is automatic
  ~AsqtadMdagM() {}

  //! Only defined on the even subset
  const OrderedSubset& subset() const {return rb[0];}

  //! Apply the operator onto a source vector
  void operator() (LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, enum PlusMinus isign) const;

private:
  Real Mass;
  // multi1d<LatticeColorMatrix> u_fat;
  // multi1d<LatticeColorMatrix> u_triple;

  AsqtadDslash D;
};


}; // End Namespace Chroma


#endif
