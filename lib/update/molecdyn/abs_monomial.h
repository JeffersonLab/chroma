// -*- C++ -*-
// $Id: abs_monomial.h,v 1.3 2004-12-21 18:01:06 bjoo Exp $

/*! @file
 * @brief Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __abs_monomial_h__
#define __abs_monomial_h__

#include "fermact.h"
#include "gaugeact.h"

#include "update/molecdyn/field_state.h"

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
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;
  };


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
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    // Compute the energies 

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;
  };

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
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;


    // Refresh all pseudofermions
    virtual void refresh(const AbsFieldState<P,Q>& field_state) = 0;
  };


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
  class ExactFermMonomial : public ExactMonomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    //! Refresh pseudofermions
    virtual void refresh(const AbsFieldState<P,Q>& field_state) = 0;
  };


  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Wilson-like fermion monomials. Not sure what these really do that
   * is new. There can be a staggered version.
   */
  template<typename P, typename Q, typename Phi>
  class ExactWilsonTypeFermMonomial : public ExactFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactWilsonTypeFermMonomial() {}

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;
    virtual void refresh(const AbsFieldState<P,Q>& field_state) = 0;

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct<Phi>& getFermAct(void) const = 0;

  };


  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Unpreconditioned Wilson-like fermion monomials. 
   */
  template<typename P, typename Q, typename Phi>
  class ExactUnprecWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactUnprecWilsonTypeFermMonomial() {}


    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    virtual void refresh(const AbsFieldState<P,Q>& field_state) = 0;

  protected:

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;
    
  };



  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Even-odd preconditioned Wilson-like fermion monomials. 
   */
  template<typename P, typename Q, typename Phi>
  class ExactEvenOddPrecWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q, Phi>
  {
  public:
     //! virtual destructor:
    ~ExactEvenOddPrecWilsonTypeFermMonomial() {}


    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) const = 0;

    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;


    virtual void refresh(const AbsFieldState<P,Q>& field_state) = 0;

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

  };


  //! Exact 2 degen flavor unpreconditioned fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecWilsonTypeFermMonomial : public ExactUnprecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecWilsonTypeFermMonomial() {}


    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const
    {
      /**** Identical code to even-odd prec case *****/
      
      // S_f  chi^dag*(M^dag*M)*chi     
      // Here, M is unprecond. operator.
      //
      // Need
      // dS_f/dU =  psi^dag * [d(M^dag)*M + M^dag*dM] * psi,  psi = (M^dag*M)^(-1)*chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = -X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X
      // where
      //    X = (M^dag*M)^(-1)*chi   Y = M*X = (M^dag)^(-1)*chi
      // In Robert's notation,  X -> psi .
      //
      const UnprecWilsonTypeFermAct<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const UnprecLinearOperator<Phi,P> > lin(FA.linOp(state));
	
      Phi X, Y;
      getX(X,s);
      
      (*lin)(Y, X, PLUS);

      lin->deriv(F, X, Y, MINUS);
      
      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      lin->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;
      
      for(int mu=0; mu < Nd; mu++) { 
      F[mu] *= Real(-1);   // This is problematic. Need convention on where to put minus
      }
    }
  
    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const {
      Phi X;
      getX(X,s);

      Double action = innerProductReal(getPhi(), X);
      return action;
    }

    
    virtual void refresh(const AbsFieldState<P,Q>& field_state) {
       // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const FermionAction<Phi>& S_f = getFermAct();
      
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
    }				    
  

  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Get (M^dagM)^{-1} phi
    virtual void getX(Phi& X, const AbsFieldState<P,Q>& s) const = 0;

  };


  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial : public ExactEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>
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
      getX(X, s);
      Double action = innerProductReal(getPhi(), X, lin->subset());
      return action;
    }

    Double S(const AbsFieldState<P,Q>& s) const {
      return S_even_even(s) + S_odd_odd(s);
    }

    //! Refresh any pseudo fermions
    /*! Can generically implement this method */
    virtual void refresh(const AbsFieldState<P,Q>& s) { 
      
      // Get at the ferion action for piece i
      const FermionAction<Phi>& S_f = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(S_f.createState(s.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator<Phi> > M(S_f.linOp(f_state));
      
      Phi eta=zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, M->subset());
      
      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now HIT IT with the ROCK!!!! (Or in this case M^{dagger})
      (*M)(getPhi(), eta, MINUS);
    }				    

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const
    {
      /**** Identical code to unprec case *****/
      
      // S_f  chi^dag*(M^dag*M)*chi
      // Here, M is precond. operator such that
      //   M*psi_o = [A_oo - D_eo*A_ee^(-1)*D_oe]*psi_o
      // NOTE: any (-1/2) from D are inside of D_eo.
      //
      // Need
      // dS_f/dU =  psi^dag * [d(M^dag)*M + M^dag*dM] * psi,  psi = (M^dag*M)^(-1)*chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = -X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X
      // where
      //    X = (M^dag*M)^(-1)*chi   Y = M*X = (M^dag)^(-1)*chi
      // In Robert's notation,  X -> psi .
      //
      const EvenOddPrecWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      Handle<const ConnectState> bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<Phi,P> > lin(FA.linOp(bc_g_state));
      // Get the X fields
      Phi X;
      getX(X, s);


      // Assume we have X and chi
      Phi Y;
      (*lin)(Y, X, PLUS);
      lin->deriv(F, X, Y, MINUS);

      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      lin->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;

      for(int mu=0; mu < Nd; mu++) { 
	F[mu] *= Real(-1);   // This is problematic. Need convention on where to put minus
      }
    }
  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get (M^dagM)^{-1} phi
    virtual void getX(Phi& X, const AbsFieldState<P,Q>& s) const = 0;

  };
}

using namespace Chroma;

#endif
