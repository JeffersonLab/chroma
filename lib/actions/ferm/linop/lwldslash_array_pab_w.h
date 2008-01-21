// -*- C++ -*-
// $Id: lwldslash_array_pab_w.h,v 3.5 2008-01-21 20:18:50 edwards Exp $
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
    PABWilsonDslashArray();

    //! Full constructor
    PABWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_);

    //! Full constructor
    PABWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_,
			 const AnisoParam_t& aniso_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_,
		const AnisoParam_t& aniso_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_,
		const multi1d<Real>& coeffs_);

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
    
    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    PrimitiveSU3Matrix* packed_gauge;
    WilsonArg wil;
    unsigned long wil_cbsize;
    int N5;
    Handle< FermBC<T,P,Q> > fbc;
  };


} // End Namespace Chroma


#endif
