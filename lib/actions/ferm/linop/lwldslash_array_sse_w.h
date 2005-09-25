// -*- C++ -*-
// $Id: lwldslash_array_sse_w.h,v 2.0 2005-09-25 21:04:29 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator array
 */

#ifndef __lwldslash_array_sse_w_h__
#define __lwldslash_array_sse_w_h__

#include "actions/ferm/linop/lwldslash_base_array_w.h"


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
                                                                                
  class SSEWilsonDslashArray : public WilsonDslashBaseArray
  {
  public:
    //! Empty constructor. Must use create later
    SSEWilsonDslashArray() {init();}

    //! Full constructor
    SSEWilsonDslashArray(const multi1d<LatticeColorMatrix>& u_, int N5_) {init(); create(u_,N5_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, int N5_);

    //! Expected length of array index
    int size() const {return N5;}

    //! No real need for cleanup here
    ~SSEWilsonDslashArray();

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

    //! Init internals
    void init();

  private:
    multi1d<PrimitiveSU3Matrix> packed_gauge;
  
    multi1d<LatticeColorMatrix> u;   // Needed only for derivative. Should find some alternative
    int N5;
  };


}; // End Namespace Chroma


#endif
