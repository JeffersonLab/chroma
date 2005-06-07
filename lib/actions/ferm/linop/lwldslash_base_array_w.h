// -*- C++ -*-
// $Id: lwldslash_base_array_w.h,v 1.2 2005-06-07 19:36:48 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator over arrays
 */

#ifndef __lwldslash_base_array_h__
#define __lwldslash_base_array_h__

#include "linearop.h"


namespace Chroma 
{ 
  //! General Wilson-Dirac dslash over arrays
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

  class WilsonDslashBaseArray : public DslashLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! No real need for cleanup here
    virtual ~WilsonDslashBaseArray() {}

    //! Subset is all here
    const OrderedSubset& subset() const {return all;}

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
    virtual void apply (LatticeFermion& chi, 
			const LatticeFermion& psi, 
			enum PlusMinus isign, int cb) const = 0;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(multi1d<LatticeColorMatrix>& ds_u, 
		       const multi1d<LatticeFermion>& chi, 
		       const multi1d<LatticeFermion>& psi, 
		       enum PlusMinus isign) const;


    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(multi1d<LatticeColorMatrix>& ds_u, 
		       const multi1d<LatticeFermion>& chi, 
		       const multi1d<LatticeFermion>& psi, 
		       enum PlusMinus isign, int cb) const;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    virtual void deriv(multi1d<LatticeColorMatrix>& ds_u, 
		       const LatticeFermion& chi, 
		       const LatticeFermion& psi, 
		       enum PlusMinus isign, int cb) const;


  protected:
    //! Get the u field
    virtual const multi1d<LatticeColorMatrix>& getU() const = 0;

  };


}; // End Namespace Chroma


#endif
