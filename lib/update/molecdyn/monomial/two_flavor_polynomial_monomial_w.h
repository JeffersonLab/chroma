// -*- C++ -*-
// $Id: two_flavor_polynomial_monomial_w.h,v 3.6 2007-03-22 17:39:23 bjoo Exp $

/*! @file
 * @brief Two flavor Monomials
 */

#ifndef __two_flavor_polynomial_monomial_w_h__
#define __two_flavor_polynomial_monomial_w_h__

#include "wilstype_polyfermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/predictor/chrono_predictor.h"

#include <typeinfo>
using namespace std;

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactPolynomialWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactPolynomialWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialWilsonTypeFermMonomial");

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
      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.polyLinOp(state));
	
      lin->deriv(F, getPhi(), getPhi(), PLUS);
      
      // F now holds derivative with respect to possibly fat links
      // now derive it with respect to the thin links if needs be
      state->deriv(F);

      monitorForces(xml_out, "Forces", F);

      pop(xml_out);
    
      END_CODE();
    }
 
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      START_CODE();

      // Heatbath all the fields
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialWilsonTypeFermMonomial");
      
      // Get at the ferion action for piece i
      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > MdagM(FA.lMdagM(f_state));
      Handle< PolyLinearOperator<Phi,P,Q> > Poly(FA.polyLinOp(f_state));
      
      Phi eta = zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, MdagM->subset());
      
      // Account for fermion BC by modifying the proposed field
      FA.getFermBC().modifyF(eta);

      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now HIT IT with the ROCK!!!! (Or in this case H)
      Phi tmp1, tmp2;
      Poly->applyA(tmp1, eta, MINUS);
      (*MdagM)(tmp2, tmp1, PLUS);

      // Solve [Q*P(Q^2)*Q]^{-1} tmp2 = phi
      // Get system solver
      Handle< PolyPrecSystemSolver<Phi> > invPolyPrec(FA.invPolyPrec(f_state, getInvParams()));

      // Do the inversion
      SystemSolverResults_t res = (*invPolyPrec)(getPhi(), tmp2);

      write(xml_out, "n_count", res.n_count);
      pop(xml_out);
    
      END_CODE();
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try {
	const TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactPolynomialWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }
    
      END_CODE();
    }
    
  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams(void) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExactPolynomial monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecPolynomialWilsonTypeFermMonomial : public TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactUnprecPolynomialWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialUnprecWilsonTypeFermMonomial");

      // Get at the fermion action
      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > Poly(FA.polyLinOp(state));

      // Energy calc doesnt use Chrono Predictor
      Phi X = zero;

      // Action on the entire lattice
      (*Poly)(X, getPhi(), PLUS);

      Double action = innerProductReal(getPhi(), X);
      
      int n_count = 0;
      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams(void) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecPolynomialWilsonTypeFermMonomial : public TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactEvenOddPrecPolynomialWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();

      Handle< FermState<Phi,P,Q> > state = FA.createState(s.getQ());

      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > Poly(FA.polyLinOp(state));

      // Energy calc doesnt use Chrono Predictor
      Phi X = zero;

      // Action on the entire lattice
      (*Poly)(X, getPhi(), PLUS);

      Double action = innerProductReal(getPhi(), X, Poly->subset());
      
      int n_count = 0;
      write(xml_out, "n_count", n_count);
      write(xml_out, "S_oo", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)  
    {
      START_CODE();

      XMLWriter& xml_out=TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialEvenOddPrecWilsonTypeFermMonomial");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams(void) const = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Constand even even determinant so can supplyt
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecConstDetPolynomialWilsonTypeFermMonomial : public TwoFlavorExactEvenOddPrecPolynomialWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactEvenOddPrecConstDetPolynomialWilsonTypeFermMonomial() {}

    //! Even even contribution (eg For this kind of Monomial it is 0)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

    // Inherit everything from Base Class
  protected:
    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;
  };

}


#endif
