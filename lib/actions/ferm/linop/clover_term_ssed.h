// -*- C++ -*-
// $Id: clover_term_ssed.h,v 1.4 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_ssed_clover_w_h__
#define __clover_term_ssed_clover_w_h__

#include "qdp_precision.h"
#if BASE_PRECISION == 32
#error "This code only works in double precision"
#endif

#include "state.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/clover_term_base_w.h"

namespace Chroma 
{ 

  //! Special structure used for triangular objects
  //  Pad to cache line length
  typedef  RScalar<REAL64> PrimitiveClovDiag[2][8];
 
  // Of the 16 elements only 15 are used. but I want to pad
  // to 16  
  typedef  RComplex<REAL64> PrimitiveClovOffDiag[2][16];

  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  class SSEDCloverTerm : public CloverTermBase<LatticeFermionD,
					       LatticeColorMatrixD>
  {
  public:
    // Typedefs to save typing
    typedef LatticeDiracFermionD3        T;
    typedef multi1d<LatticeColorMatrixD3>  P;
    typedef multi1d<LatticeColorMatrixD3>  Q;

    //! Empty constructor. Must use create later
    SSEDCloverTerm();
    //! Free the internals
    ~SSEDCloverTerm();

    //! Create from another
    void create(Handle< FermState<T,P,Q> > fs, 
		const CloverFermActParams& param_, 
		const SSEDCloverTerm& from);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
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
    Double cholesDet(int cb) const ;

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
    void apply (LatticeDiracFermionD3& chi, const LatticeDiracFermionD3& psi, enum PlusMinus isign, int cb) const;


    void applySite(LatticeDiracFermionD3& chi, const LatticeDiracFermionD3& psi, enum PlusMinus isign, int site) const;

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(LatticeColorMatrixD3& B, int mat, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<LatticeColorMatrixD3>& f, const Real& diag_mass);

    //! Invert the clover term on cb using LDL^\dagger decomp
    void ldagdlinv(LatticeDouble& tr_log_diag, int cb);

    //! Invert the clover term on cb using Cholesky decomp
    void chlclovms(LatticeDouble& log_diag, int cb);

    //! Get the u field
    const multi1d<LatticeColorMatrixD3>& getU() const {return u;}

    //! Calculates Tr_D ( Gamma_mat L )
    Real getCloverCoeff(int mu, int nu) const;

  private:
    Handle< FermBC<T,P,Q> >      fbc;
    multi1d<LatticeColorMatrixD3>  u;
    CloverFermActParams          param;
    LatticeDouble                  tr_log_diag_; // Fill this out during create
                                                  // but save the global sum until needed.
    multi1d<bool> choles_done;   // Keep note of whether the decomposition has been done
                                 // on a particular checkerboard. 
    PrimitiveClovDiag*           tri_diag;
    PrimitiveClovOffDiag *       tri_off_diag;
  };


  
} // End Namespace Chroma


#endif
