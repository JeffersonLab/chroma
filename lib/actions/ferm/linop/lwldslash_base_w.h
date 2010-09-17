// -*- C++ -*-
// $Id: lwldslash_base_w.h,v 3.4 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_base_h__
#define __lwldslash_base_h__

#include "linearop.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "tower_array.h"
#include "pq_traits.h"

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
  template <typename T, typename P, typename Q>
  class WilsonDslashBase : 
    public DslashLinearOperator<T,P,Q>
  {
  public:

    //! No real need for cleanup here
    virtual ~WilsonDslashBase() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const;

    virtual void
    deriv(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
	  const Tower<T>& chi, const Tower<T>& psi, 
	  const P& p,
	  enum PlusMinus isign) const;

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const ;



    virtual void 
    deriv(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
	  const Tower<T>& chi, const Tower<T>& psi, 
	  const P& p,
	  enum PlusMinus isign, int cb) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

    //! Creation routine
    virtual void create(Handle< FermState<T,P,Q> > state) =0;

    //! Apply checkerboarded linear operator
    /*! 
     * To avoid confusion (especially of the compilers!), call the checkerboarded
     * apply instead of operator()
     */

    //! Apply checkerboarded linear operator to tower
    void applyTower( Tower<T>& chi, 
			const Tower<T>& psi,
			const P& mom,
			enum PlusMinus isign, 
			int cb)  {

      QDPIO::cout << "Applying with cb="<<cb<<endl;
      if( psi.getHeight() != chi.getHeight() ) { 
	QDPIO::cout << "chi and psi have incompatible heights" << endl;
	QDP_abort(1);
      }

      int N = psi.getHeight();

      // Always do this for now.
      apply(chi[0], psi[0], isign,cb);

      // Preserve the links 
      Handle<FermState<T,P,Q> >  old_state  = getState();
      
      // Now do chi[i] = Sum_{j=0}^i Binom(i,j)  D^{j} psi{i-j}
      // Where D^{j} is D^j = (-P)^j U

      for(int i=1; i < N; i++) { 

	create(old_state);

	// Start off with 
	// J=0 step. Original state so D(U) = D(old_state), psi[i]
	apply(chi[i], psi[i], isign,cb);

	// Get the links from the old state
	Q links(Nd);
	for(int mu=0; mu < Nd; mu++) { 
	  links[mu] = old_state->getLinks()[mu];
	}

	// J=1...i
	for(int j=1; j <= i; j++) { 
	  
	  // Lift: Make D( (-P)^j U ): multiply links by -P
	  for(int mu=0; mu < Nd; mu++) { 
	    Q tmp_u(1);
	    tmp_u[0]= links[mu];
	    links[mu] = -mom[mu]*tmp_u[0];
	  }

	  Handle< FermState<T,P,Q> > newfs( new PeriodicFermState<T,P,Q>(links));
     
	  create(newfs);

	  // Get D_j psi[ i - j]
	  T result=zero;
	  apply(result, psi[i-j], isign,cb);

	  // Add with Binom Coeff
	  chi[i][rb[cb]] += Real(psi.Choose(i,j))*result;
	} // Next j
      } // Next i

      create(old_state); // Restore Operator 
    }


  protected:
    //! Get the anisotropy parameters
    virtual const multi1d<Real>& getCoeffs() const = 0;
    
    //! Get at the state
    virtual Handle<FermState<T,P,Q> > getState() const = 0;

  };

  template<typename T>
  class HalfFermionType{};

  template<>
  class HalfFermionType<LatticeFermionF>  { 
  public:
    typedef LatticeHalfFermionF Type_t;
  };

  template<>
  class HalfFermionType<LatticeFermionD>  { 
  public:
    typedef LatticeHalfFermionD Type_t;
  };


  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  template<typename T, typename P, typename Q>
  void
  WilsonDslashBase<T,P,Q>::deriv(P& ds_u,
			  const T& chi, const T& psi, 
			  enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    P ds_tmp; 
    deriv(ds_u, chi, psi, isign, 0);
    deriv(ds_tmp, chi, psi, isign, 1);
    ds_u += ds_tmp;

    END_CODE();
  }
  template<typename T, typename P, typename Q>
  void
  WilsonDslashBase<T,P,Q>::deriv(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				 const Tower<T>& chi, const Tower<T>& psi, 
				 const P& p,
				 enum PlusMinus isign) const
  {
    START_CODE();
    

    TowerArray<typename PQTraits<Q>::Base_t> ds_tmp(ds_u.getHeight());
    deriv(ds_u, chi, psi, p, isign, 0);
    deriv(ds_tmp, chi, psi, p, isign, 1);

    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] += ds_tmp[mu];
    }
    END_CODE();
  }


  //! Take deriv of D
  /*! \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$  */
  template<typename T, typename P, typename Q>
  void 
  WilsonDslashBase<T,P,Q>::deriv(P& ds_u,
				 const T& chi, const T& psi, 
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();

    ds_u.resize(Nd);

    const multi1d<Real>& anisoWeights = getCoeffs();

    for(int mu = 0; mu < Nd; ++mu) 
    {
      // Break this up to use fewer expressions:
      T temp_ferm1;
      typename HalfFermionType<T>::Type_t tmp_h;

      switch (isign) 
      {
      case PLUS:
      {
	// Undaggered: Minus Projectors
	switch(mu) 
	{ 
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Minus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Minus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Minus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Minus(tmp_h);
	  break;
	default:
	  break;
	};
	
      }
      break;

      case MINUS:
      {
	// Daggered: Plus Projectors
	typename HalfFermionType<T>::Type_t tmp_h;
	switch(mu) 
	{
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Plus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Plus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Plus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Plus(tmp_h);
	  break;
	default:
	  break;
	};
      }
      break;

      default:
	QDP_error_exit("unknown case");
      }

      // QDP Shifts the whole darn thing anyhow
      T temp_ferm2 = shift(temp_ferm1, FORWARD, mu);
      P temp_mat;
      temp_mat.resize(1);

      // This step supposedly optimised in QDP++
      (temp_mat[0])[rb[cb]] = traceSpin(outerProduct(temp_ferm2,chi));
    
      // Just do the bit we need.
      ds_u[mu][rb[cb]] = anisoWeights[mu] * temp_mat[0];
      ds_u[mu][rb[1-cb]] = zero;    
    }
    (*this).getFermBC().zero(ds_u);

    END_CODE();
  }


  //! Take deriv of D
  /*! \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$  */
  template<typename T, typename P, typename Q>
  void 
  WilsonDslashBase<T,P,Q>::deriv(TowerArray<typename PQTraits<Q>::Base_t>& ds_u,
				 const Tower<T>& chi, const Tower<T>& psi, 
				 const P& p,
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if(chi.getHeight() != psi.getHeight() ) { 
      QDPIO::cout << "chi and psi towers have incompatible heights: chi.getHeight()= " << chi.getHeight() << " psi.getHeight()=" << psi.getHeight() << endl; 
      QDP_abort(1);
    }
    if(ds_u.getHeight() != psi.getHeight()) { 
      QDPIO::cout << "ds_u and psi have incompatible heights: dsu.getHeight() = " << ds_u.getHeight() <<  "  psi.getHeight()=" << psi.getHeight() << endl;
      QDP_abort(1);
    }

    int N=psi.getHeight();
     
    const multi1d<Real>& anisoWeights = getCoeffs();

    


    for(int mu = 0; mu < Nd; ++mu) {

      Tower<T> temp_ferm1(N);
      for(int level =0; level < N; level++ ) { 

	typename HalfFermionType<T>::Type_t tmp_h;

	switch (isign) {

	case PLUS:
	  {
	    // Undaggered: Minus Projectors
	    switch(mu) {

	    case 0:
		tmp_h[rb[1-cb]] = spinProjectDir0Minus(psi[level]);
		temp_ferm1[level][rb[1-cb]] = spinReconstructDir0Minus(tmp_h);
		break;
	    case 1:
	      tmp_h[rb[1-cb]] = spinProjectDir1Minus(psi[level]);
	      temp_ferm1[level][rb[1-cb]] = spinReconstructDir1Minus(tmp_h);
	      break;
	    case 2:
	      tmp_h[rb[1-cb]] = spinProjectDir2Minus(psi[level]);
	      temp_ferm1[level][rb[1-cb]] = spinReconstructDir2Minus(tmp_h);
	      break;
	    case 3:
	      tmp_h[rb[1-cb]] = spinProjectDir3Minus(psi[level]);
	      temp_ferm1[level][rb[1-cb]] = spinReconstructDir3Minus(tmp_h);
	      break;
	    default:
	      break;
	    };
	    
	  }
	  break;

	case MINUS:
	  {
	    // Daggered: Plus Projectors
	    typename HalfFermionType<T>::Type_t tmp_h;
	    switch(mu) 
	      {
	      case 0:
		tmp_h[rb[1-cb]] = spinProjectDir0Plus(psi[level]);
		temp_ferm1[level][rb[1-cb]] = spinReconstructDir0Plus(tmp_h);
		break;
	      case 1:
		tmp_h[rb[1-cb]] = spinProjectDir1Plus(psi[level]);
		temp_ferm1[level][rb[1-cb]] = spinReconstructDir1Plus(tmp_h);
		break;
	      case 2:
		tmp_h[rb[1-cb]] = spinProjectDir2Plus(psi[level]);
		temp_ferm1[level][rb[1-cb]] = spinReconstructDir2Plus(tmp_h);
		break;
	      case 3:
		tmp_h[rb[1-cb]] = spinProjectDir3Plus(psi[level]);
		temp_ferm1[level][rb[1-cb]] = spinReconstructDir3Plus(tmp_h);
		break;
	      default:
		break;
	      };
	  }
	  break;

	default:
	  QDP_error_exit("unknown case");
	}
      } // end loop over levels - spin projections all done

      // QDP Shifts the whole darn thing anyhow
      Tower<T> temp_ferm2(N);
      temp_ferm2 = shiftTower(temp_ferm1, FORWARD, mu);
      
      Tower<typename PQTraits<Q>::Base_t> temp_mat(N);

      // Funky outer product 
      spinTraceOuterProduct(temp_mat, temp_ferm2, chi, rb[cb]);
	
      //      // This step supposedly optimised in QDP++
      //  (temp_mat[0])[rb[cb]] = traceSpin(outerProduct(temp_ferm2,chi));
      
      // Just do the bit we need.
      for(int level=0; level < N; level++) { 
	ds_u[mu][level][rb[cb]] = anisoWeights[mu] * temp_mat[level];
	ds_u[mu][level][rb[1-cb]] = zero;    

      }
    }
    
    (*this).getFermBC().zero(ds_u);

    END_CODE();
  }




  //! Return flops performed by the operator()
  template<typename T, typename P, typename Q>
  unsigned long 
  WilsonDslashBase<T,P,Q>::nFlops() const {return 1320;}

#if 0

 class WilsonDslashBaseF : 
    public DslashLinearOperator<LatticeFermionF,
				multi1d<LatticeColorMatrixF>,
				multi1d<LatticeColorMatrixF> >
  {
  public:
    typedef LatticeFermionF T;
    typedef multi1d<LatticeColorMatrixF> P;
    typedef multi1d<LatticeColorMatrixF> Q;

    //! No real need for cleanup here
    virtual ~WilsonDslashBaseF() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const {
      QDPIO::cout << "Not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const {
      QDPIO::cout << "Not Implemented" << endl;
    }

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! Get the anisotropy parameters
    virtual const multi1d<Real>& getCoeffs() const = 0;
  };

#endif

} // End Namespace Chroma


#endif
