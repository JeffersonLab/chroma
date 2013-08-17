// -*- C++ -*-
// $Id: two_flavor_polyprec_monomial_w.h,v 3.7 2007-03-22 17:39:23 bjoo Exp $

/*! @file
 * @brief Two flavor Monomials
 */

#ifndef __two_flavor_polyprec_monomial_w_h__
#define __two_flavor_polyprec_monomial_w_h__

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
  class TwoFlavorExactPolyPrecWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactPolyPrecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolyPrecWilsonTypeFermMonomial");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f =  chi^dag*(M^dag*M)^(-1)*chi     
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
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.polyPrecLinOp(state));
	
      Phi X;

      // Get X out here
      int n_count = this->getX(X,s);

      lin->deriv(F, X, X, PLUS);
      
      for(int mu=0; mu < F.size(); ++mu)
	F[mu] *= Real(-1);

      // F now holds derivative with respect to possibly fat links
      // now derive it with respect to the thin links if needs be
      state->deriv(F);

      write(xml_out, "n_count", n_count);
      monitorForces(xml_out, "Forces", F);

      pop(xml_out);
    
      END_CODE();
    }
 
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      START_CODE();

      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< LinearOperator<Phi> > H(FA.hermitianLinOp(f_state));
      Handle< PolyLinearOperator<Phi,P,Q> > Poly(FA.polyLinOp(f_state));
      
      Phi eta=zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, H->subset());
      
      // Account for fermion BC by modifying the proposed field
      FA.getFermBC().modifyF(eta);

      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now HIT IT with the ROCK!!!! (Or in this case H)
      Phi tmp;
      (*H)(tmp, eta, PLUS);
      Poly->applyA(getPhi(), tmp, PLUS);
      QDPIO::cout << "TwoFlavPolyWilson4DMonomial: resetting Predictor after field refresh" << endl;
      getMDSolutionPredictor().reset();
    
      END_CODE();
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try {
	const TwoFlavorExactPolyPrecWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const TwoFlavorExactPolyPrecWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactPolyPrecWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }

      END_CODE();
    }

    //! Reset predictors
    virtual void resetPredictors(void) {
      getMDSolutionPredictor().reset();

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

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;

    //! Get (Q*P(Q^2)*Q)^{-1} phi
    virtual int getX(Phi& X, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Grab the fermact
      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Guess is passed in
   
      // Get linop
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.polyPrecLinOp(state));

      // Solve [Q*P(Q^2)*Q]^{-1} X = phi
      // Get system solver
      Handle< PolyPrecSystemSolver<Phi> > invPolyPrec(FA.invPolyPrec(state, getInvParams()));

      // Do the inversion...
      (getMDSolutionPredictor())(X, *M, getPhi());
      SystemSolverResults_t res = (*invPolyPrec)(X, getPhi());
      (getMDSolutionPredictor()).newVector(X);

      END_CODE();

      return res.n_count;
    }

  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExactPolyPrec monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecPolyPrecWilsonTypeFermMonomial : public TwoFlavorExactPolyPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactUnprecPolyPrecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactPolyPrecUnprecWilsonTypeFermMonomial");

      Phi X;
      
      // Energy calc doesnt use Chrono Predictor
      X = zero;
      QDPIO::cout << "TwoFlavPolyWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();
      int n_count = this->getX(X,s);

      // Action on the entire lattice
      Double action = innerProductReal(getPhi(), X);
      
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

    //! Get at the chronological predcitor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecPolyPrecWilsonTypeFermMonomial : public TwoFlavorExactPolyPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactEvenOddPrecPolyPrecWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const PolyWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();

      Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      // Get the X fields
      Phi X;

      // Action calc doesnt use chrono predictor use zero guess
      X[ lin->subset() ] = zero;

      // getX noe always uses chrono predictor. Best to Nuke it therefore
      QDPIO::cout << "TwoFlavPolyWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();
      int n_count = this->getX(X, s);
      Double action = innerProductReal(getPhi(), X, lin->subset());
      
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
      push(xml_out, "TwoFlavorExactPolyPrecEvenOddPrecWilsonTypeFermMonomial");

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

    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor(void) = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Constand even even determinant so can supplyt
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecConstDetPolyPrecWilsonTypeFermMonomial : public TwoFlavorExactEvenOddPrecPolyPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    virtual ~TwoFlavorExactEvenOddPrecConstDetPolyPrecWilsonTypeFermMonomial() {}

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
