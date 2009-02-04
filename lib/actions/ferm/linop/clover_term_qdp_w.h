// -*- C++ -*-
// $Id: clover_term_qdp_w.h,v 3.7 2009-02-04 15:20:31 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_qdp_w_h__
#define __clover_term_qdp_w_h__

#include "state.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
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
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Empty constructor. Must use create later
    QDPCloverTerm();

    //! No real need for cleanup here
    ~QDPCloverTerm() {}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    virtual void create(Handle< FermState<T,P,Q> > fs,
			const CloverFermActParams& param_,
			const QDPCloverTerm& from_);

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
    void apply (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int cb) const;


    void applySite(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, int site) const;

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(LatticeColorMatrix& B, int mat, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<LatticeColorMatrix>& f, const Real& diag_mass);

    //! Invert the clover term on cb
    void chlclovms(LatticeReal& log_diag, int cb);
    void ldagdlinv(LatticeReal& tr_log_diag, int cb);

    //! Get the u field
    const multi1d<LatticeColorMatrix>& getU() const {return u;}

    //! Calculates Tr_D ( Gamma_mat L )
    Real getCloverCoeff(int mu, int nu) const;

  private:
    Handle< FermBC<T,P,Q> >      fbc;
    multi1d<LatticeColorMatrix>  u;
    CloverFermActParams          param;
    LatticeReal                  tr_log_diag_; // Fill this out during create
                                                  // but save the global sum until needed.
    multi1d<bool> choles_done;   // Keep note of whether the decomposition has been done
                                 // on a particular checkerboard. 
    multi1d<PrimitiveClovTriang>  tri;
    
  };

  struct QDPCloverApplyStruct { 
    LatticeFermion& chi;
    const LatticeFermion& psi;
    const multi1d<PrimitiveClovTriang>& tri;
    int cb;
  };
  
  void QDPCloverDispatchFunction(int lo, int hi, int myId, QDPCloverApplyStruct* a);


 
 
} // End Namespace Chroma


#endif
