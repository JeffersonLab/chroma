// -*- C++ -*-
// $Id: lwldslash_qdpopt_w.h,v 3.3 2009-04-17 02:05:34 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_qdpopt_h__
#define __lwldslash_qdpopt_h__

#include "state.h"
#include "io/aniso_io.h"
#include "actions/ferm/linop/lwldslash_base_w.h"


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

  template<typename T, typename P, typename Q>
  class QDPWilsonDslashOptT : public WilsonDslashBase<T,P,Q>
  {
  public:

    //! Empty constructor. Must use create later
    QDPWilsonDslashOptT();

    //! Full constructor
    QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state);

    //! Full constructor with anisotropy
    QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state,
		       const AnisoParam_t& aniso_);

    //! Full constructor with general coefficients
    QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state,
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
    ~QDPWilsonDslashOptT() {}

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
    void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    Handle< FermBC<T,P,Q> >  fbc;
    Q   u;
  };


#if 0
   template<typename T>
  class HalfFermionType{};

  template<>
  class HalfFermionType<LatticeFermionF> {
  public: 
    typedef LatticeHalfFermionF Type_t;
  };

  template<>
  class HalfFermionType<LatticeFermionD> {
  public:
    typedef LatticeHalfFermionD Type_t;
  };
#endif

  //! General Wilson-Dirac dslash
  /*! \ingroup linop
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

  //! Empty constructor
  template<typename T, typename P, typename Q>
  QDPWilsonDslashOptT<T,P,Q>::QDPWilsonDslashOptT() {}
  
  //! Full constructor
  template<typename T, typename P, typename Q>
  QDPWilsonDslashOptT<T,P,Q>::QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state)
  {
    create(state);
  }
  
  //! Full constructor with anisotropy
  template<typename T, typename P, typename Q>
  QDPWilsonDslashOptT<T,P,Q>::QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state,
					 const AnisoParam_t& aniso_) 
  {
    create(state, aniso_);
  }

  //! Full constructor with general coefficients
  template<typename T, typename P, typename Q>
  QDPWilsonDslashOptT<T,P,Q>::QDPWilsonDslashOptT(Handle< FermState<T,P,Q> > state,
					 const multi1d<Real>& coeffs_)
  {
    create(state, coeffs_);
  }

  //! Creation routine
  template<typename T, typename P, typename Q>
  void QDPWilsonDslashOptT<T,P,Q>::create(Handle< FermState<T,P,Q> > state)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, cf);
  }

  //! Creation routine with anisotropy
  template<typename T, typename P, typename Q>
  void QDPWilsonDslashOptT<T,P,Q>::create(Handle< FermState<T,P,Q> > state,
					  const AnisoParam_t& anisoParam) 
  {
    START_CODE();
    
    create(state, makeFermCoeffs(anisoParam));
    
    END_CODE();
  }

  //! Full constructor with general coefficients
  template<typename T, typename P, typename Q>
  void QDPWilsonDslashOptT<T,P,Q>::create(Handle< FermState<T,P,Q> > state,
					  const multi1d<Real>& coeffs_)
  {
    START_CODE();

    // Save a copy of the aniso params original fields and with aniso folded in
    coeffs = coeffs_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "QDPWilsonDslashOpt: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fold in anisotropy
    u.resize(Nd);
    for(int mu=0; mu < u.size(); ++mu) {
      u[mu] = (state->getLinks())[mu];
    }
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
    }
  }


  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result				                (Write)
   *  \param psi	      Pseudofermion field				(Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  template<typename T, typename P, typename Q>
  void 
  QDPWilsonDslashOptT<T,P,Q>::apply (T& chi, const T& psi, 
			     enum PlusMinus isign, int cb) const
  {
    START_CODE();
#if QDP_NC == 3
    /*     F 
     *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
     *     mu          mu                    mu
     */
    /*     B           +
     *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
     *     mu          mu                       mu
     */
    // Recontruct the bottom two spinor components from the top two
    /*                        F           B
     *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
     *                        mu          mu
     */

    /* Why are these lines split? An array syntax would help, but the problem is deeper.
     * The expression templates require NO variable args (like int's) to a function
     * and all args must be known at compile time. Hence, the function names carry
     * (as functions usually do) the meaning (and implicit args) to a function.
     */

    int otherCB = (cb == 0 ? 1 : 0);

    switch (isign)
    {
    case PLUS: 
      {


	typename HalfFermionType<T>::Type_t tmp;
	typename HalfFermionType<T>::Type_t tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Minus(u[0]*tmp2);
	
	
	tmp[rb[otherCB]]  = spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(u[1]*tmp2);

	tmp[rb[otherCB]]  = spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(tmp2);	
	
      }


      break;

    case MINUS:
      {

	typename HalfFermionType<T>::Type_t tmp;
	typename HalfFermionType<T>::Type_t tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(u[1]*tmp2);

	tmp[rb[otherCB]] = spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Minus(psi);    
	tmp2 = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(tmp2);	
	
      }

      break;
    }

    getFermBC().modifyF(chi, QDP::rb[cb]);
#else
    QDPIO::cerr<<"lwldshash_qdpopt_w: not implemented for NC!=3\n";
    QDP_abort(13);
#endif

    END_CODE();
  }

  typedef QDPWilsonDslashOptT< LatticeFermion, 
			       multi1d<LatticeColorMatrix>, 
			       multi1d<LatticeColorMatrix> > QDPWilsonDslashOpt;

  typedef QDPWilsonDslashOptT< LatticeFermionF, 
			       multi1d<LatticeColorMatrixF>, 
			       multi1d<LatticeColorMatrixF> > QDPWilsonDslashOptF;

  typedef QDPWilsonDslashOptT< LatticeFermionD, 
			       multi1d<LatticeColorMatrixD>, 
			       multi1d<LatticeColorMatrixD> > QDPWilsonDslashOptD;


} // End Namespace Chroma


#endif
