// -*- C++ -*-
// $Id: abs_monomial.h,v 1.1 2004-12-13 22:41:00 edwards Exp $

/*! @file
 * @brief Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __abs_monomial_h__
#define __abs_monomial_h__

#include "fermact.h"
#include "gaugeact.h"

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

    //! clone function for virtual copy constructs
    virtual Monomial<P,Q>* clone(void) const = 0;

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

    //! clone function for virtual copy constructs
    virtual ExactMonomial<P,Q>* clone(void) const = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    // Compute the energies 

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const
    {
      return mesKE(s) + mesPE(s);
    }

    //! The total energy
    virtual void mesE(Double& KE, Double& PE, const AbsFieldState<P,Q>& s) const 
    {
      KE = mesKE(s);
      PE = mesPE(s);
    }

    //! The Kinetic Energy
    virtual Double mesKE(const AbsFieldState<P,Q>& s) const = 0;

    //! The Potential Energy 
    virtual Double mesPE(const AbsFieldState<P,Q>& s) const = 0;
  };



  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * The fermion hierachy would splits at the very top into
   * inexact and exact monomials. An exact monomial can be used
   * for an inexact algorithm, but not vice-versa.
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial : public ExactMonomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial() {}

    //! Refresh any pseudo fermions
    virtual void refresh(const AbsFieldState<P,Q>& s) = 0;

    //! clone function for virtual copy constructs
    virtual ExactFermMonomial<P,Q,Phi>* clone(void) const = 0;
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

    //! clone function for virtual copy constructs
    virtual ExactWilsonTypeFermMonomial<P,Q>* clone(void) const = 0;

    //! Problem here, the dsdq cannot be defined yet - fermact derivs not yet defined
    /*! It is not possible to define a generic function yet */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi>& getFermAct() const = 0;

    // Compute the energies. The fermionic energy can now finally be computed

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const
    {
      return mesKE(s) + mesPE(s) + mesFE(s);
    }

    //! The total energy
    virtual void mesE(Double& KE, Double& GE, Double& FE,
		      const AbsFieldState<P,Q>& s) const 
    {
      KE = mesKE(s);
      GE = mesGE(s);
      FE = mesFE(s);
    }

    //! The Kinetic Energy
    virtual Double mesKE(const AbsFieldState<P,Q>& s) const = 0;

    //! The Potential Energy 
    virtual Double mesGE(const AbsFieldState<P,Q>& s) const = 0;

    //! The Fermionic Energy 
    virtual Double mesFE(const AbsFieldState<P,Q>& s) const = 0;
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

    //! clone function for virtual copy constructs
    virtual ExactUnprecWilsonTypeFermMonomial<P,Q,Phi>* clone(void) const = 0;

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Compute dsdq for the system... 
    /*! Now can define forces in terms of fermacts */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;
  };



  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup actions
   *
   * Even-odd preconditioned Wilson-like fermion monomials. 
   */
  template<typename P, typename Q, typename Phi>
  class ExactEvenOddPrecWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q>
  {
  public:
     //! virtual destructor:
    ~ExactEvenOddPrecWilsonTypeFermMonomial() {}

    //! clone function for virtual copy constructs
    virtual ExactEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>* clone(void) const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Compute dsdq for the system... 
    /*! Now can define forces in terms of fermacts */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const = 0;
  };



  //! Exact 2 degen flavor unpreconditioned fermact monomial
  /*! @ingroup actions
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecWilsonTypeFermMonomial : public ExactUnprecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactUnprecWilsonTypeFermMonomial() {}

    //! clone function for virtual copy constructs
    virtual ExactUnprecWilsonTypeFermMonomial<P,Q,Phi>* clone(void) const = 0;

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Refresh any pseudo fermions
    /*! Can generically implement this method */
    virtual void refresh(const AbsFieldState<P,Q>& s) = 0;

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
      const EvenOddPrecWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<T,P> > lin(FA->linOp(s));

      // .....
      // Do something to get  X  and chi
      // .....
      Phi X, chi;
      chi = X = zero;

      // Assume we have X and chi
      Phi Y;
      lin(Y, X, PLUS);
      lin->deriv(F, X, Y, MINUS);

      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      lin->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;

      F *= -1;   // This is problematic. Need convention on where to put minus
    }
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
    ~ExactEvenOddPrecWilsonTypeFermMonomial() {}

    //! clone function for virtual copy constructs
    virtual ExactEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>* clone(void) const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Refresh any pseudo fermions
    /*! Can generically implement this method */
    virtual void refresh(const AbsFieldState<P,Q>& s) = 0;

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

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<T,P> > lin(FA->linOp(s));

      // .....
      // Do something to get  X  and chi
      // .....
      Phi X, chi;
      chi = X = zero;

      // Assume we have X and chi
      Phi Y;
      lin(Y, X, PLUS);
      lin->deriv(F, X, Y, MINUS);

      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      lin->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;

      F *= -1;   // This is problematic. Need convention on where to put minus
    }
  };

}

using namespace Chroma;

#endif
