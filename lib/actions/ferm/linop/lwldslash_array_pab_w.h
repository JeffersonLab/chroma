// -*- C++ -*-
// $Id: lwldslash_array_pab_w.h,v 1.3 2005-06-17 15:17:53 bjoo Exp $
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
    PABWilsonDslashArray(const multi1d<LatticeColorMatrix>& u_, int N5_) {create(u_,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, int N5_);

    //! Expected length of array index
    int size() const {return N5;}

    //! No real need for cleanup here
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
    //! Get the u field
    const multi1d<LatticeColorMatrix>& getU() const {return u;}

  private:
    PrimitiveSU3Matrix* packed_gauge;
    multi1d<LatticeColorMatrix> u; // For derivative
    WilsonArg wil;
    unsigned long wil_cbsize;
    int N5;
  };


}; // End Namespace Chroma


#endif
