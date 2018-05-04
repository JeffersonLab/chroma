// -*- C++ -*-
/*! \file
 *  \brief Wilson Dslash linear operator over arrays
 */

#ifndef __lwldslash_array_h__
#define __lwldslash_array_h__

#include "state.h"
#include "actions/ferm/linop/lwldslash_base_array_w.h"


namespace Chroma 
{ 
  //! General Wilson-Dirac dslash of arrays
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

  class QDPWilsonDslashArray : public WilsonDslashBaseArray
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Empty constructor. Must use create later
    QDPWilsonDslashArray() {}

    //! Full constructor
    QDPWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_)
      {create(state,N5_);}

    //! Full constructor
    QDPWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_,
			 const AnisoParam_t& aniso_)
      {create(state,N5_,aniso_);}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_,
		const AnisoParam_t& aniso_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state,
		int N5_,
		const multi1d<Real>& coeffs_);

    //! Expected length of array index
    int size() const {return N5;}

    //! No real need for cleanup here
    ~QDPWilsonDslashArray() {}

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
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
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply (LatticeFermion& chi, 
		const LatticeFermion& psi, 
		enum PlusMinus isign, int cb) const;
     
    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    int N5;
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    multi1d<LatticeColorMatrix> u;  // fold in anisotropy
    Handle< FermBC<T,P,Q> > fbc;
  };


} // End Namespace Chroma


#endif
