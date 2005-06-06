// -*- C++ -*-
// $Id: lwldslash_array_w.h,v 1.1 2005-06-06 03:47:13 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator over arrays
 */

#ifndef __lwldslash_array_h__
#define __lwldslash_array_h__

#include "actions/ferm/linop/lwldslash_base_array_w.h"


namespace Chroma 
{ 
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

  class QDPWilsonDslashArray : public WilsonDslashBaseArray
  {
  public:
    //! Empty constructor. Must use create later
    QDPWilsonDslashArray() {}

    //! Full constructor
    QDPWilsonDslashArray(const multi1d<LatticeColorMatrix>& u_, int N5_) {create(u_,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, int N5_);

    //! Expected length of array index
    int size() const {return N5;}

    //! No real need for cleanup here
    ~QDPWilsonDslashArray() {}

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

  protected:
    //! Get the u field
    const multi1d<LatticeColorMatrix>& getU() const {return u;}

  private:
    multi1d<LatticeColorMatrix> u;
    int N5;
  };


}; // End Namespace Chroma


#endif
