// -*- C++ -*-
// $Id: lwldslash_w_cppf.h,v 3.2 2009-05-20 19:28:11 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_cpp_float_h__
#define __lwldslash_cpp_float_h__

#include "actions/ferm/linop/lwldslash_base_w.h"
#include "io/aniso_io.h"
#include "state.h"
#include "cpp_dslash.h"
#include "cpp_dslash_qdp_packer.h" 

using namespace CPlusPlusWilsonDslash;

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
                                                                                
  class CPPWilsonDslashF : public WilsonDslashBase<LatticeFermionF,
						  multi1d<LatticeColorMatrixF>, 
						  multi1d<LatticeColorMatrixF> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermionF               T;
    typedef multi1d<LatticeColorMatrixF>  P;
    typedef multi1d<LatticeColorMatrixF>  Q;

    //! Empty constructor. Must use create later
    CPPWilsonDslashF();

    //! Full constructor
    CPPWilsonDslashF(Handle< FermState<T,P,Q> > state);

    //! Full constructor with anisotropy
    CPPWilsonDslashF(Handle< FermState<T,P,Q> > state,
		    const AnisoParam_t& aniso_);

    //! Full constructor with general coefficients
    CPPWilsonDslashF(Handle< FermState<T,P,Q> > state,
		    const multi1d<Real>& coeffs_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state);

    //! Creation routine with anisotropy
    void create(Handle< FermState<T,P,Q> > state, 
		const AnisoParam_t& aniso_);

    //! Full constructor with general coefficients
    void create(Handle< FermState<T,P,Q> > state, 
		const multi1d<Real>& coeffs_);

    //! No real need for cleanup here
    ~CPPWilsonDslashF();

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
    void apply(T& chi, const T& psi, 
	       enum PlusMinus isign, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

    //! Init internals
    void init();

  private:
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    multi1d<PrimitiveSU3MatrixF> packed_gauge;  // fold in anisotropy
    Handle< FermBC<T,P,Q> > fbc;
    Handle< Dslash<float> > D;

  };


} // End Namespace Chroma


#endif
