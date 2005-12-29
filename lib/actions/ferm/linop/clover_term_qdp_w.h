// -*- C++ -*-
// $Id: clover_term_qdp_w.h,v 2.2 2005-12-29 05:37:36 edwards Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_qdp_w_h__
#define __clover_term_qdp_w_h__

#include "actions/ferm/linop/clover_term_base_w.h"

namespace Chroma 
{ 

  //! Special structure used for triangular objects
  struct PrimitiveClovTriang
  {
    RScalar<REAL>   diag[2][2*Nc];
    RComplex<REAL>  offd[2][2*Nc*Nc-Nc];
  };

  // Reader/writers
  /*! \ingroup linop */
  void read(XMLReader& xml, const string& path, PrimitiveClovTriang& param);

  /*! \ingroup linop */
  void write(XMLWriter& xml, const string& path, const PrimitiveClovTriang& param);


  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  class QDPCloverTerm : public CloverTermBase
  {
  public:
    //! No real need for cleanup here
    ~QDPCloverTerm() {}

    //! Creation routine
    void create(const multi1d<LatticeColorMatrix>& u_, 	
		const CloverFermActParams& param_);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     */
    void choles(int cb);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     * \return logarithm of the determinant  
     */
    Double cholesDet(int cb);

    /**
     * Apply a dslash
     *
     * Performs the operation
     *
     *  chi <-   (L + D + L^dag) . psi
     *
     * where
     *   L       is a lower triangular matrix
     *   D       is the real diagonal. (stored together in type TRIANG)
     *
     * Arguments:
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT vector               (Read) 
     */
    void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;

  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<LatticeColorMatrix>& f, const Real& diag_mass);

    //! Invert the clover term on cb
    void chlclovms(bool DetP, Double& logdet, int cb);

    //! Get the u field
    const multi1d<LatticeColorMatrix>& getU() const {return u;}

  private:
    multi1d<LatticeColorMatrix>  u;
    CloverFermActParams          param;

    multi1d<PrimitiveClovTriang>  tri;
  };


} // End Namespace Chroma


#endif
