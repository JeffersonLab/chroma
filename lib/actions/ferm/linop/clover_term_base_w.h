// -*- C++ -*-
// $Id: clover_term_base_w.h,v 2.2 2006-01-01 05:12:30 edwards Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_base_w_h__
#define __clover_term_base_w_h__

#include "linearop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"


namespace Chroma 
{ 
  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  class CloverTermBase : public DslashLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! No real need for cleanup here
    virtual ~CloverTermBase() {}

    //! Subset is all here
    const OrderedSubset& subset() const {return all;}

    //! Creation routine
    virtual void create(const multi1d<LatticeColorMatrix>& u_, 	
			const CloverFermActParams& param_);

    //! Invert
    /*!
     * Computes the inverse of the term on cb using Cholesky
     */
    virtual void choles(int cb);

    //! Invert
    /*!
     * Computes the inverse of the term on cb using Cholesky
     *
     * \return logarithm of the determinant  
     */
    virtual Double cholesDet(int cb);

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   chi^dag * \dot(D} * psi  
     */
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
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
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign, int cb) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! Calculates Tr_D ( Gamma_mat L )
    virtual void triacntr(LatticeColorMatrix& B, int mat, int cb) const = 0;

    //! Get the u field
    virtual const multi1d<LatticeColorMatrix>& getU() const = 0;
  };


} // End Namespace Chroma


#endif
