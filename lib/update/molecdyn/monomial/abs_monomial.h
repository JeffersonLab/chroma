// -*- C++ -*-
// $Id: abs_monomial.h,v 1.1 2005-01-13 15:55:04 bjoo Exp $

/*! @file
 * @brief Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __abs_monomial_h__
#define __abs_monomial_h__

#include "fermact.h"
#include "gaugeact.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/chrono_predictor.h"
namespace Chroma
{
  //! An abstract monomial class, for inexact algorithms
  /*! @ingroup actions
   *
   * Inexact in this case means energy computation is not supported,
   * (in an inexact algorithm sense -- obviously it is weird to have
   * a hamiltonian where you cannot compute the energy. We may need
   * to think more about this)
   * 
   * This serves the following purpose. It definees 
   * an interface for computing the total force 
   * and can refresh the momenta,
   * 
   * 
   * We don't specify how the momenta is refreshed. It is "virtual".
   * HMD type algorithms will porbably use gaussian noise. 
   * GHMD type algorithms will mix the previous momenta some. How
   * to do that will be encoded in the derived class, probably 
   * through the constructor.
   * 
   * 
   * For this it needs to know the types of coordinates and the momenta
   * so that it can act on the right kind of state.
   */
  template<typename P, typename Q>
  class Monomial
  {
  public:
    //! virtual destructor:
    virtual ~Monomial() {}

    //! Compute dsdq for the system... 
    /*! Not specified how to actually do this s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermion fields if any
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) =0 ;

    //! Copy pseudofermion fields from another monomial...
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Abstract monomial class, for exact algorithms
  /*! @ingroup actions
   *
   * Now define similar classes for exact algorithms.
   * These are basically the same as before but they can compute
   * energies too. Do these need to inerit?
   * Yes. Reason: We can always give it to an inexact algorithm through
   * a downcast. In that case the energy calculations would be hidden.
   */
  template<typename P, typename Q>
  class ExactMonomial : public Monomial<P, Q> 
  {
  public:
    //! virtual destructor:
    virtual ~ExactMonomial() {}

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    // Compute the energies 

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const  = 0;

    //! Refresh pseudofermion fields if any
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * The fermion hierachy would splits at the very top into
   * inexact and exact monomials. An exact monomial can be used
   * for an inexact algorithm, but not vice-versa.
   */
 
  /* Unfortunately we need to template on the Phi-s because
     we need that template for the FermActs */
  template<typename P, typename Q, typename Phi>
  class FermMonomial : public Monomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~FermMonomial() {}

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) = 0;

    // Refresh all pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0 ;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * The fermion hierachy would splits at the very top into
   * inexact and exact monomials. An exact monomial can be used
   * for an inexact algorithm, but not vice-versa.
   *
   * Unfortunately we need to template on the Phi-s because
   *  we need that template for the FermActs 
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial : public ExactMonomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Exact fermionic monomials with pseudofermions living in 4D
   *
   * We need to template on the Phi-s because of the fermacts
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial4D : public ExactFermMonomial<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial4D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Exact fermionic monomials with pseudofermions living in 4D
   *
   * We need to template on the Phi-s because of the fermacts
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial5D : public ExactFermMonomial<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Wilson-like fermion monomials. Not sure what these really do that
   * is new. There can be a staggered version.
   */
  template<typename P, typename Q, typename Phi>
  class ExactWilsonTypeFermMonomial : public ExactFermMonomial4D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Wilson-like fermion monomials. Not sure what these really do that
   * is new. There can be a staggered version.
   */
  template<typename P, typename Q, typename Phi>
  class ExactWilsonTypeFermMonomial5D : public ExactFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct5D<Phi,P>& getFermAct(void) const = 0;

  };



  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f  chi^dag*(M^dag*M)^(-1)*chi     
      // Here, M is some operator.
      //
      // Need
      // dS_f/dU = - psi^dag * [d(M^dag)*M + M^dag*dM] * psi,  psi = (M^dag*M)^(-1)*chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = -X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X
      // where
      //    X = (M^dag*M)^(-1)*chi   Y = M*X = (M^dag)^(-1)*chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const DiffLinearOperator<Phi,P> > lin(FA.linOp(state));
	
      Phi X, Y;

      // Get X out here
      (getMDSolutionPredictor())(X);
      int n_count = getX(X,s);
      (getMDSolutionPredictor()).newVector(X);

      (*lin)(Y, X, PLUS);

      lin->deriv(F, X, Y, MINUS);
      
      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      lin->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;
 

      for(int mu=0; mu < Nd; mu++) { 
	 F[mu] *= Real(-1);   // IS THIS SIGN CORRECT???
      }

    }
  
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const WilsonTypeFermAct<Phi,P>& S_f = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(S_f.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator<Phi> > M(S_f.linOp(f_state));
      
      Phi eta=zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, M->subset());
      
      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now HIT IT with the ROCK!!!! (Or in this case M^{dagger})
      (*M)(getPhi(), eta, MINUS);
      getMDSolutionPredictor().reset();
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) {
      try {
	const TwoFlavorExactWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const TwoFlavorExactWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }

      // Resetting pseudofermion fields implies resetting the chrono predictor
      getMDSolutionPredictor().reset();
    }
  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(Phi& X, const AbsFieldState<P,Q>& s) const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor fermact monomial in extra dimensions
  /*! @ingroup actions
   *
   * Exact 2 degen flavor fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactWilsonTypeFermMonomial5D : public ExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*V*(M^dag*M)^(-1)*V^dag*chi     
      // Here, M is some 5D operator and V is the Pauli-Villars field
      //
      // Need
      // dS_f/dU =  chi^dag * dV * (M^dag*M)^(-1) * V^dag * chi 
      //         -  chi^dag * V * (M^dag*M)^(-1) * [d(M^dag)*M + M^dag*dM] * V^dag * (M^dag*M)^(-1) * chi
      //         +  chi^dag * V * (M^dag*M)^(-1) * d(V^dag) * chi 
      //
      //         =  chi^dag * dV * psi
      //         -  psi^dag * [d(M^dag)*M + M^dag*dM] * psi
      //         +  psi^dag * d(V^dag) * chi 
      //
      // where  psi = (M^dag*M)^(-1) * V^dag * chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = chi^dag*\dot(V)*X - X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X + X*\dot{V}^dag*chi
      // where
      //    X = (M^dag*M)^(-1)*V^dag*chi   Y = M*X = (M^dag)^(-1)*V^dag*chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Get linear operator
      Handle< const DiffLinearOperator<multi1d<Phi>, P> > M(FA.linOp(state));
	
      // Get Pauli-Villars linear operator
      Handle< const DiffLinearOperator<multi1d<Phi>, P> > PV(FA.linOpPV(state));
	
      // Get/construct the pseudofermion solution
      multi1d<Phi> X(FA.size()), Y(FA.size());

      (getMDSolutionPredictor())(X);
      int n_count = getX(X,s);
      (getMDSolutionPredictor()).newVector(X);

      (*M)(Y, X, PLUS);

      // First PV contribution
      PV->deriv(F, getPhi(), X, PLUS);

      // First interior term
      P F_tmp;
      M->deriv(F_tmp, X, Y, MINUS);
      F -= F_tmp;   // NOTE SIGN
      
      // fold M^dag into X^dag ->  Y  !!
      M->deriv(F_tmp, Y, X, PLUS);
      F -= F_tmp;   // NOTE SIGN
      
      // Last PV contribution
      PV->deriv(F_tmp, X, getPhi(), MINUS);
      F += F_tmp;   // NOTE SIGN

    }
  
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator< multi1d<Phi> > > M(FA.linOp(f_state));
      
      // Get Pauli-Villars linear operator
      Handle< const LinearOperator< multi1d<Phi> > > PV(FA.linOpPV(f_state));
	
      const int N5 = FA.size();
      multi1d<Phi> eta(N5);
      eta = zero;
      
      // Fill the eta field with gaussian noise
      for(int s=0; s < N5; ++s)
	gaussian(eta[s], M->subset());
      
      // Temporary: Move to correct normalisation
      for(int s=0; s < N5; ++s)
	eta[s][M->subset()] *= sqrt(0.5);
      
      // Build  phi = V * (V^dag*V)^(-1) * M^dag * eta
      multi1d<Phi> tmp(N5);
      (*M)(tmp, eta, MINUS);

      // Solve  (V^dag*V)*eta = tmp
      int n_pv_count = getXPV(eta, tmp, field_state);

      // Finally, get phi
      (*PV)(getPhi(), eta, PLUS);

      // Reset the chronological predictor
      getMDSolutionPredictor().reset();
    }				    

    virtual void setInternalFields(const Monomial<P,Q>& m) {
      try {
	const TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>& fm = dynamic_cast< const TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>& >(m);

	// Do a resize here -- otherwise if the fields have not yet
	// been refreshed there may be trouble
	getPhi().resize(fm.getPhi().size());

	for(int i=0 ; i < fm.getPhi().size(); i++) { 
	  (getPhi())[i] = (fm.getPhi())[i];
	}
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactWilsonTypeFermMonomial5D" << endl;
	QDP_abort(1);
      }

      // Reset the chronological predictor
      getMDSolutionPredictor().reset();
    }
  

  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P>& getFermAct(void) const = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s) const = 0;

    //! Get X = (PV^dag*PV)^{-1} eta
    virtual int getXPV(multi1d<Phi>& X, const multi1d<Phi>& eta, const AbsFieldState<P,Q>& s) const = 0;

    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor(void) = 0;
   };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecWilsonTypeFermMonomial : public TwoFlavorExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)   const
    {
      Phi X;
      
      // Energy calc doesnt use Chrono Predictor
      X = zero;

      int n_count = getX(X,s);

      // Action on the entire lattice
      Double action = innerProductReal(getPhi(), X);
      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int  getX(Phi& X, const AbsFieldState<P,Q>& s) const = 0;
    
    //! Get at the chronological predcitor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;
  };


  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial : public TwoFlavorExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) const
    {
      const EvenOddPrecWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      Handle<const ConnectState> bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<Phi,P> > lin(FA.linOp(bc_g_state));
      // Get the X fields
      Phi X;

      // Action calc doesnt use chrono predictor use zero guess
      X[ lin->subset() ] = zero;

      int n_count = getX(X, s);
      Double action = innerProductReal(getPhi(), X, lin->subset());
      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)  const {
      return S_even_even(s) + S_odd_odd(s);
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get (M^dagM)^{-1} phi
    virtual int getX(Phi& X, const AbsFieldState<P,Q>& s) const  = 0;

    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;
  };



  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned fermact monomial living in extra dimensions
  /*! @ingroup actions
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecWilsonTypeFermMonomial5D : public TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const
    {
           // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(FA.createState(s.getQ()));
      Handle< const LinearOperator< multi1d<Phi> > > PV(FA.linOpPV(f_state));
 
      multi1d<Phi> X(FA.size());
      multi1d<Phi> tmp(FA.size());

      // Paranoia -- to deal with subsets.
      tmp = zero; 

      // Energy calc does not use chrono predictor
      X = zero;
      // X is now (M^dagM)^{-1} V^{dag} phi
      int n_count = getX(X,s);

      // tmp is now V (M^dag M)^{-1} V^{dag} phi
      (*PV)(tmp, X, PLUS);

      // Action on the entire lattice
      Double action = zero;
      for(int s=0; s < getFermAct().size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s]);
      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P>& getFermAct(void) const = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s) const  = 0;

    //! Get X = (PV^dag*PV)^{-1} eta
    virtual int getXPV(multi1d<Phi>& X, const multi1d<Phi>& eta, const AbsFieldState<P,Q>& s) const = 0;

    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor(void) = 0;
  };


  //! Exact 2 degen flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup actions
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D : public TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) const
    {
      const EvenOddPrecWilsonTypeFermAct5D<Phi,P>& FA = getFermAct();

      Handle<const ConnectState> bc_g_state(FA.createState(s.getQ()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<multi1d<Phi>,P> > lin(FA.linOp(bc_g_state));

      Handle< const EvenOddPrecLinearOperator<multi1d<Phi>,P> > PV(FA.linOpPV(bc_g_state));
      // Get the X fields
      multi1d<Phi> X(FA.size());

      // X is now (M^dag M)^{-1} V^dag phi

      // Chrono predictor not used in energy calculation
      X = zero;
      int n_count = getX(X, s);

      multi1d<Phi> tmp(FA.size());
      (*PV)(tmp, X, PLUS);

      Double action = zero;
      // Total odd-subset action. NOTE: QDP has norm2(multi1d) but not innerProd
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s], lin->subset());
      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s) const  {
      return S_even_even(s) + S_odd_odd(s);
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P>& getFermAct() const = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s) const  = 0;

    //! Get X = (PV^dag*PV)^{-1} eta
    virtual int getXPV(multi1d<Phi>& X, const multi1d<Phi>& eta, const AbsFieldState<P,Q>& s) const  = 0;
  };



}

using namespace Chroma;

#endif
