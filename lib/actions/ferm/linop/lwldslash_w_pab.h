// -*- C++ -*-
// $Id: lwldslash_w_pab.h,v 2.1 2005-12-18 23:53:26 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_pab_h__
#define __lwldslash_pab_h__

#include "actions/ferm/linop/lwldslash_base_w.h"

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
    PABWilsonDslash(const multi1d<LatticeColorMatrix>& u_);

    //! Full constructor with anisotropy
    PABWilsonDslash(const multi1d<LatticeColorMatrix>& u_, 
		    const AnisoParam_t& aniso_);

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_);

    //! Creation routine with anisotropy
    void create(const multi1d<LatticeColorMatrix>& u_, const AnisoParam_t& aniso_);

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
    void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;

  protected:
    //! Get the anisotropy parameters
    const AnisoParam_t& getAnisoParam() const {return anisoParam;}

  private:
    AnisoParam_t  anisoParam;
    PrimitiveSU3Matrix* packed_gauge;  // fold in anisotropy
    WilsonArg wil;
    unsigned long wil_cbsize;
  };


}; // End Namespace Chroma


#endif
