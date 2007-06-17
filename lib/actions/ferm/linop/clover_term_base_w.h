// -*- C++ -*-
// $Id: clover_term_base_w.h,v 3.7 2007-06-17 02:25:16 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_base_w_h__
#define __clover_term_base_w_h__

#include "linearop.h"


namespace Chroma 
{ 
  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  class CloverTermBase : public DslashLinearOperator<LatticeFermion, 
			 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! No real need for cleanup here
    virtual ~CloverTermBase() {}

    //! Subset is all here
    const Subset& subset() const {return all;}


    virtual void applySite(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int site) const = 0;

    //! Invert
    /*!
     * Computes the inverse of the term on cb using Cholesky
     */
    virtual void choles(int cb) = 0;

    //! Invert
    /*!
     * Computes the determinant of the term
     *
     * \return logarithm of the determinant  
     */
    virtual Double cholesDet(int cb) const = 0;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
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
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign, int cb) const;

    //! Take derivative of TrLn D
    void derivTrLn(multi1d<LatticeColorMatrix>& ds_u, 
		   enum PlusMinus isign, int cb) const;


    void deriv_loops(const int u, const int mu, const int cb,
		     LatticeColorMatrix& ds_u_mu,
		     LatticeColorMatrix& ds_u_nu,
		     const LatticeColorMatrix& Lambda) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

    //! Calculates Tr_D ( Gamma_mat L )
    virtual void triacntr(LatticeColorMatrix& B, int mat, int cb) const = 0;

  protected:

    //! Get the u field
    virtual const multi1d<LatticeColorMatrix>& getU() const = 0;

    //! get the clover coefficient 
    virtual Real getCloverCoeff(int mu, int nu) const = 0;

  };


} // End Namespace Chroma


#endif
