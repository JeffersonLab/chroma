// -*- C++ -*-
// $Id: lwldslash_base_3d_w.h,v 3.2 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief 3D Wilson Dslash linear operator
 */

#ifndef __lwldslash_3d_base_h__
#define __lwldslash_3d_base_h__

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "linearop.h"
#include "io/aniso_io.h"


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

  class WilsonDslash3DBase : public DslashLinearOperator<LatticeFermion, 
           multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! No real need for cleanup here
    virtual ~WilsonDslash3DBase() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const;


    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      3D Checkerboard of chi vector               (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! Get the anisotropy parameters
    virtual const multi1d<Real>& getCoeffs() const = 0;
  };


} // End Namespace Chroma

#endif 
#endif
#endif

#endif
