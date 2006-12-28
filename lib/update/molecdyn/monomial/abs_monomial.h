// -*- C++ -*-
// $Id: abs_monomial.h,v 3.2 2006-12-28 17:34:00 bjoo Exp $

/*! @file
 * @brief Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __abs_monomial_h__
#define __abs_monomial_h__

#include "wilstype_fermact_w.h"
#include "gaugeact.h"

#include "update/molecdyn/field_state.h"
#include "io/xmllog_io.h"

namespace Chroma
{
  //! An abstract monomial class, for inexact algorithms
  /*! @ingroup monomial
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

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Abstract monomial class, for exact algorithms
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermion fields if any
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };

  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

  };

}


#endif
