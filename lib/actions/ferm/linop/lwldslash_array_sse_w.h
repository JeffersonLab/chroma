// -*- C++ -*-
// $Id: lwldslash_array_sse_w.h,v 3.2 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator array
 */

#ifndef __lwldslash_array_sse_w_h__
#define __lwldslash_array_sse_w_h__

#include "actions/ferm/linop/lwldslash_base_array_w.h"
#include "state.h"
#include "sse_dslash_qdp_packer.h" 

using namespace SSEDslash;

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
                                                                                
  class SSEWilsonDslashArray : public WilsonDslashBaseArray
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Empty constructor. Must use create later
    SSEWilsonDslashArray();

    //! Full constructor
    SSEWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_);

    //! Full constructor
    SSEWilsonDslashArray(Handle< FermState<T,P,Q> > state,
			 int N5_,
			 const AnisoParam_t& aniso_);

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
    ~SSEWilsonDslashArray();

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
    void apply (multi1d<LatticeFermion>& chi, 
		const multi1d<LatticeFermion>& psi, 
		enum PlusMinus isign, int cb) const;

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
    void apply (LatticeFermion& chi, 
		const LatticeFermion& psi, 
		enum PlusMinus isign, int cb) const;
    
    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Init internals
    void init();

    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    multi1d<PrimitiveSU3Matrix> packed_gauge;
    int N5;
    Handle< FermBC<T,P,Q> > fbc;
  };


} // End Namespace Chroma


#endif
