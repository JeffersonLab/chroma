// -*- C++ -*-
// $Id: lwldslash_array_pab_w.h,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator over arrays
 */

#ifndef __lwldslash_array_pab_h__
#define __lwldslash_array_pab_h__

#include "actions/ferm/linop/lwldslash_base_array_w.h"
#include <wfm.h>

namespace Chroma 
{ 

  typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;

  //! General Wilson-Dirac dslash of arrays
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

  class PABWilsonDslashArray : public WilsonDslashBaseArray
  {
  public:
    //! Empty constructor. Must use create later
    PABWilsonDslashArray() {}

    //! Full constructor
    PABWilsonDslashArray(const multi1d<LatticeColorMatrix>& u_, 
			 int N5_);

    //! Full constructor
    PABWilsonDslashArray(const multi1d<LatticeColorMatrix>& u_, 
			 int N5_,
			 const AnisoParam_t& aniso_);

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		int N5_);

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 
		int N5_,
		const AnisoParam_t& aniso_);

    //! Expected length of array index
    int size() const {return N5;}

    //! Destroy / deal with refcounting
    ~PABWilsonDslashArray();

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply (multi1d<LatticeFermion>& chi, 
		const multi1d<LatticeFermion>& psi, 
		enum PlusMinus isign, int cb) const;

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply (LatticeFermion& chi, 
		const LatticeFermion& psi, 
		enum PlusMinus isign, int cb) const;
    
  protected:
    //! Get the anisotropy parameters
    const AnisoParam_t& getAnisoParam() const {return anisoParam;}

  private:
    AnisoParam_t  anisoParam;
    PrimitiveSU3Matrix* packed_gauge;
    WilsonArg wil;
    unsigned long wil_cbsize;
    int N5;
  };


} // End Namespace Chroma


#endif
