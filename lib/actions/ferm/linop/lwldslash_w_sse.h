// -*- C++ -*-
// $Id: lwldslash_w_sse.h,v 1.14 2005-02-20 02:55:48 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_sse_h__
#define __lwldslash_sse_h__

#include "actions/ferm/linop/lwldslash_base_w.h"


namespace Chroma 
{ 
  typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;

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
                                                                                
  class SSEWilsonDslash : public WilsonDslashBase
  {
  public:
    //! Empty constructor. Must use create later
    SSEWilsonDslash() {init();}

    //! Full constructor
    SSEWilsonDslash(const multi1d<LatticeColorMatrix>& u_) {init();create(u_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_);

    //! No real need for cleanup here
    ~SSEWilsonDslash();

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
    void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;

  protected:
    //! Get the u field
    const multi1d<LatticeColorMatrix>& getU() const {return u;}

    //! Init internals
    void init();

  private:
    multi1d<PrimitiveSU3Matrix> packed_gauge;
  
    multi1d<LatticeColorMatrix> u;   // Needed only for derivative. Should find some alternative
  };


}; // End Namespace Chroma


#endif
