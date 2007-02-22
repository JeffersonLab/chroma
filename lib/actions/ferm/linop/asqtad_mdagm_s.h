// -*- C++ -*-
// $Id: asqtad_mdagm_s.h,v 3.1 2007-02-22 21:11:46 bjoo Exp $
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

  class AsqtadMdagM : public DiffLinearOperator<LatticeStaggeredFermion, 
                             multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Partial constructor
    AsqtadMdagM() {}

    //! Full constructor
    AsqtadMdagM(Handle<AsqtadConnectStateBase> state, const Real& Mass_)
      {create(state, Mass_);}

    //! Creation routine
    void create(Handle<AsqtadConnectStateBase> state, const Real& Mass_);
    
    //! Destructor is automatic
    ~AsqtadMdagM() {}

    //! Only defined on the even subset
    const Subset& subset() const {return rb[0];}

    //! Return the fermion BC object for this linear operator
    const FermBC<LatticeStaggeredFermion,
		 multi1d<LatticeColorMatrix>,
		 multi1d<LatticeColorMatrix> >& getFermBC() const {return D.getFermBC();}

    //! Apply the operator onto a source vector
    void operator() (LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, enum PlusMinus isign) const;

  private:
    Real Mass;
    AsqtadDslash D;
  };


} // End Namespace Chroma


#endif
