// -*- C++ -*-
// $Id: lwldslash_w.h,v 1.10 2004-12-20 03:59:31 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_h__
#define __lwldslash_h__

#include "actions/ferm/linop/lwldslash_base_w.h"

using namespace QDP;

namespace Chroma 
{ 
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

  class QDPWilsonDslash : public WilsonDslashBase
  {
  public:
    //! Empty constructor. Must use create later
    QDPWilsonDslash() {}

    //! Full constructor
    QDPWilsonDslash(const multi1d<LatticeColorMatrix>& u_) {create(u_);}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_);

    //! No real need for cleanup here
    ~QDPWilsonDslash() {}

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

  private:
    multi1d<LatticeColorMatrix> u;
// Real CoeffWilsr_s;
  };


}; // End Namespace Chroma

using namespace Chroma;

#endif
