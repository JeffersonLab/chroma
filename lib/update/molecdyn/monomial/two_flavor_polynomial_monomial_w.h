// -*- C++ -*-
// $Id: two_flavor_polynomial_monomial_w.h,v 2.1 2006-02-09 22:26:41 edwards Exp $

/*! @file
 * @brief Two flavor Monomials
 */

#ifndef __two_flavor_polynomial_monomial_w_h__
#define __two_flavor_polynomial_monomial_w_h__

#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "actions/ferm/invert/invert.h"
#include "invtype.h"

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
    ~TwoFlavorExactPolynomialWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
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
      const PolyWilsonTypeFermAct<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const DiffLinearOperator<Phi,P> > lin(FA.polyLinOp(state));
	
      lin->deriv(F, getPhi(), getPhi(), PLUS);
      
      // F now holds derivative with respect to possibly fat links
      // now derive it with respect to the thin links if needs be
      state->deriv(F);

      Double F_sq = norm2(F);

      int n_count = 0;
      write(xml_out, "n_count", n_count);
      write(xml_out, "F_sq", F_sq);
      pop(xml_out);
    }
 
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      // Heatbath all the fields
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialWilsonTypeFermMonomial");
      
      // Get at the ferion action for piece i
      const PolyWilsonTypeFermAct<Phi,P>& S_f = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(S_f.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator<Phi> > MdagM(S_f.lMdagM(f_state));
      Handle< const LinearOperator<Phi> > Poly(S_f.polyLinOp(f_state));
      Handle< const LinearOperator<Phi> > PolyPrec(S_f.polyPrecLinOp(state));
      
      Phi eta=zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, MdagM->subset());
      
      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now HIT IT with the ROCK!!!! (Or in this case H)
      Phi tmp1, tmp2;
      Poly->ApplyA(tmp1, eta, MINUS);
      (*MdagM)(tmp2, tmp1, PLUS);

      // Solve [Q*P(Q^2)*Q]^{-1} tmp2 = phi
      // Do the inversion...
      int n_count = invert(getPhi(), *PolyPrec, tmp2);

      write(xml_out, "n_count", n_count);
      pop(xml_out);
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) {
      try {
	const TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactPolynomialWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }
    }
  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    // Get X = (Q*P(Q^2)*Q)^{-1} eta
    virtual int invert(Phi& X, const LinearOperator<Phi>& M, const Phi& eta) const
    {
      const InvertParam_t& inv_param = getInvParams();

      int n_count = 0;

      // Do the inversion...
      switch( inv_param.invType) {
      case CG_INVERTER:
      {
	// Solve M X = eta
	InvCG1(M, eta, X, inv_param.RsdCG, inv_param.MaxCG, n_count);
	QDPIO::cout << "2Flav::invert,  n_count = " << n_count << endl;
      }
      break;
      default:
      {
	QDPIO::cerr << "Currently only CG Inverter is implemented" << endl;
	QDP_abort(1);
      }
      break;
      };
      
      return n_count;
    }

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
  class TwoFlavorExactPolynomialUnprecWilsonTypeFermMonomial : public TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactPolynomialUnprecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialUnprecWilsonTypeFermMonomial");

      // Get at the fermion action
      const PolyWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > state(FA.createState(s.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator<Phi> > Poly(FA.polyLinOp(state));

      // Energy calc doesnt use Chrono Predictor
      Phi X = zero;

      // Action on the entire lattice
      (*Poly)(X, getPhi(), PLUS);

      Double action = innerProductReal(getPhi(), X);
      
      int n_count = 0;
      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactPolynomialEvenOddPrecWilsonTypeFermMonomial : public TwoFlavorExactPolynomialWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactPolynomialEvenOddPrecWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "S_odd_odd");

      const PolyWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      Handle<const ConnectState> state = FA.createState(s.getQ());

      // Create a linear operator
      Handle< const LinearOperator<Phi> > Poly(FA.polyLinOp(state));

      // Energy calc doesnt use Chrono Predictor
      Phi X = zero;

      // Action on the entire lattice
      (*Poly)(X, getPhi(), PLUS);

      Double action = innerProductReal(getPhi(), X, Poly->subset());
      
      int n_count = 0;
      write(xml_out, "n_count", n_count);
      write(xml_out, "S_oo", action);
      pop(xml_out);

      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)  {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactPolynomialEvenOddPrecWilsonTypeFermMonomial");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
      return action;

    }

  protected:
    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

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
  class TwoFlavorExactPolynomialEvenOddPrecConstDetWilsonTypeFermMonomial : public TwoFlavorExactPolynomialEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactPolynomialEvenOddPrecConstDetWilsonTypeFermMonomial() {}

    //! Even even contribution (eg For this kind of Monomial it is 0)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

    // Inherit everything from Base Class
  protected:
    //! Get at fermion action
    virtual const PolyWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;
  };

}


#endif
