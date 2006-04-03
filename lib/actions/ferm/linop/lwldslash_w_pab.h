// -*- C++ -*-
// $Id: lwldslash_w_pab.h,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_pab_h__
#define __lwldslash_pab_h__

#include "actions/ferm/linop/lwldslash_base_w.h"
#include "state.h"

#include <wfm.h>


namespace Chroma 
{ 
  typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;

  namespace PABDslashEnv { 
    extern int refcount;
  }

  //! General Wilson-Dirac dslash
  /*!
   * \ingroup linop
   *
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */

  class PABWilsonDslash : public WilsonDslashBase
  {
  public:
    //! Empty constructor. Must use create later
    PABWilsonDslash();

    //! Full constructor
    PABWilsonDslash(Handle< FermState<T,P,Q> > state);

    //! Full constructor with anisotropy
    PABWilsonDslash(Handle< FermState<T,P,Q> > state,
		    const AnisoParam_t& aniso_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state);

    //! Creation routine with anisotropy
    void create(Handle< FermState<T,P,Q> > state, 
		const AnisoParam_t& aniso_);

    //! No real need for cleanup here
    ~PABWilsonDslash();

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb	    Checkerboard of OUTPUT vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply(LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const AnisoParam_t& getAnisoParam() const {return anisoParam;}

  private:
    AnisoParam_t  anisoParam;
    PrimitiveSU3Matrix* packed_gauge;  // fold in anisotropy
    WilsonArg wil;
    unsigned long wil_cbsize;
    Handle< FermBC<T,P,Q> > fbc;
  };


} // End Namespace Chroma


#endif
