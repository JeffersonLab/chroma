// -*- C++ -*-
//  $Id: asq_dsl_s.h,v 3.3 2007-09-01 18:53:55 edwards Exp $
/*! \file
 *  \brief The "asq" or "asqtad" dslash operator D'
 */

#ifndef __asqdslash_h__
#define __asqdslash_h__

#include "linearop.h"
#include "actions/ferm/fermstates/asqtad_state.h"


namespace Chroma 
{ 
  //! The "asq" or "asqtad" dslash operator D'
  /*!
   * \ingroup linop
   *
   * This routine is specific to staggered fermions!
   *
   * Description:
   *
   * This routine applies the "asq" or "asqtad" operator D' to Psi, 
   * putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \                     F                     
   *   chi(x)  :=  >  isign eta  (x) [U  (x) psi(x+mu)
   *	       /             mu	     mu
   *	       ---
   *	       mu=0
   *
   *			+ c_3 U  (x) U  (x+mu) U  (x+2mu) psi(x+3mu) ]
   *                             mu     mu        mu
   *
   *	             Nd-1
   *	             ---
   *	             \                      +F
   *                -    >  isign eta  (x)  [U  (x-mu) psi(x-mu)
   *	             /             mu	    mu
   *	             ---
   *	             mu=0
   *
   *                             +      +          +
   *			+ c_3 U  (x) U  (x-2mu) U  (x-3mu) psi(x-3mu) ]
   *                             mu     mu         mu
   * Note the KS phase factors are already included in the U's!
   */

  class QDPStaggeredDslash : public DslashLinearOperator< 
    LatticeStaggeredFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {  
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Empty constructor. Must use create later
    QDPStaggeredDslash() {}
 
    //! Full constructor
    QDPStaggeredDslash(Handle<AsqtadConnectStateBase> state_)
    {create(state_);}
 
    //! Creation routine  
    void create(Handle<AsqtadConnectStateBase> state_);
 
    //! No real need for cleanup here
    ~QDPStaggeredDslash() {}

    /*! Arguments:
     *
     *  \param u_fat     Fat7 links                	  		(Read)
     *  \param u_triple  triple links					(Read)
     *  \param psi       Pseudofermion field - Source		        (Read)
     *  \param isign     D' or D'^+  ( +1 | -1 ) respectively		(Read)
     *  \param cb	       Checkerboard of OUTPUT vector			(Read) 
     */
    void apply (LatticeStaggeredFermion& chi, const LatticeStaggeredFermion& psi, 
		enum PlusMinus isign, int cb) const;
  
    //! Subset is all here
    const Subset& subset() const {return all;}
    
    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return state->getBC();}

  private:
    Handle<AsqtadConnectStateBase> state;
  };

} // End Namespace Chroma


#endif


