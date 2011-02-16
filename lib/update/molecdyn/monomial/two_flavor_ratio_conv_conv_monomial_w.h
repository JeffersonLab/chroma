// -*- C++ -*-
// $Id: two_flavor_ratio_conv_conv_monomial_w.h,v 3.3 2009-06-02 15:56:40 bjoo Exp $

/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __two_flavor_ratio_conv_conv_monomial_w_h__
#define __two_flavor_ratio_conv_conv_monomial_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/predictor/chrono_predictor.h"

#include <typeinfo>
using namespace std;

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 2 flavor RatioConvConv type monomial
  /*! @ingroup monomial
   *
   * Exact 2 flavor RatioConvConv type  monomial. 
   *
   * Can supply default dsdq()
   *                    pseudoferm refresh
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactRatioConvConvWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactRatioConvConvWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*M_prec(M^dag*M)^{-1}M^{dag}_prec*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactRatioConvConvWilsonTypeFermMonomial");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = phi^dag  M_prec (M^dag*M)^(-1) M_prec^dag phi     
      // Here, M is some operator.
      //
      // Need
      // dS_f/dU = phi^dag d(M_prec) phi
      //         - psi^dag * [d(M^dag)*M + M^dag*dM] * psi,  
      //           psi^dag d(M_prec^dag) phi,   psi = (M^dag*M)^(-1)M_prec^dag*phi
      //
      // In Balint's notation, the result is  
      // \dot{S} = chi^dag d(M_prec) X 
      //          - X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X
      //          + X^dag*\dot{M_prec}^\dag chi
      // where
      //    X = (M^dag*M)^(-1)M_prec^\dag*chi   Y = M*X = (M^dag)^(-1)M_prec^\dag *chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();          // for M
      const WilsonTypeFermAct<Phi,P,Q>& FA_prec = getDenomFermAct();  // for M_prec

      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Get system solver
      Handle< MdagMSystemSolver<Phi> > invMdagM(FA.invMdagM(state,getNumerInvParams()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(state));	
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(FA_prec.linOp(state));

      Phi X=zero;
      Phi Y=zero;

      // Need MdagM for CG based predictor
      Handle< DiffLinearOperator<Phi,P,Q> > MdagM(FA.lMdagM(state));
      Phi M_dag_prec_phi;

      // M_dag_prec phi = M^{dag}_prec \phi - the RHS
      (*M_prec)(M_dag_prec_phi, getPhi(), MINUS);

      //(getMDSolutionPredictor())(X, *MdagM, M_dag_prec_phi);

      // Solve MdagM X = eta
      SystemSolverResults_t res = (*invMdagM)(X, M_dag_prec_phi, getMDSolutionPredictor());

      // (getMDSolutionPredictor()).newVector(X);
      
      (*M)(Y, X, PLUS);

      // \phi^{\dagger} \dot(M_prec) X
      M_prec->deriv(F, getPhi(), X, PLUS);
      
      // - X^{\dagger} \dot( M^{\dagger}) Y
      P F_tmp;
      M->deriv(F_tmp, X, Y, MINUS);
      F -= F_tmp;
 
      // - Y^{\dagger} \dot( M ) X
      M->deriv(F_tmp, Y, X, PLUS);
      F -= F_tmp;

      // + X^{\dagger} \dot(M_prec)^dagger \phi
      M_prec->deriv(F_tmp, X, getPhi(), MINUS);
      F += F_tmp;

      // F now holds derivative with respect to possibly fat links
      // now derive it with respect to the thin links if needs be
      state->deriv(F);

      write(xml_out, "n_count", res.n_count);
      monitorForces(xml_out, "Forces", F);

      pop(xml_out);

      END_CODE();
    }
  
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      START_CODE();

      // Heatbath all the fields
      
      // Get at the fermion action for the expensive matrix
      const WilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();

      // Get the fermion action for the preconditioner
      const WilsonTypeFermAct<Phi,P,Q>& FA_prec = getDenomFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > state(FA.createState(field_state.getQ()));
      
      // Create a linear operator for the Expensive op
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(state));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(FA_prec.linOp(state));

      Phi eta = zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, M->subset());
      
      // Account for fermion BC by modifying the proposed field
      FA.getFermBC().modifyF(eta);

      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Now we have to get the refreshment right:
      //
      // We have  \eta^{\dag} \eta = \phi M_prec (M^dag M )^-1 M^dag_prec \phi
      //
      //  so that \eta = (M^\dag)^{-1} M^{dag}_prec \phi
      //
      //  So need to solve M^{dag}_prec \phi = M^{\dag} \eta
      //
      // Which we can solve as 
      //
      //      M^{dag}_prec M_prec ( M_prec^{-1} ) \phi = M^{\dag} \eta
      //
      // I will dedicate a function to this:
      //
      // 
      //
      // prepare the source
      // Zero out everything - to make sure there is no junk
      // in uninitialised places
      Phi eta_tmp = zero;
      Phi phi_tmp = zero;
      getPhi() = zero;

      (*M)(eta_tmp, eta, MINUS);  // M^\dag \eta

      // Get system solver
      Handle< MdagMSystemSolver<Phi> > invMdagM(FA_prec.invMdagM(state, getDenomInvParams()));

      // Solve MdagM_prec X = eta
      SystemSolverResults_t res = (*invMdagM)(phi_tmp, eta_tmp);

      (*M_prec)(getPhi(), phi_tmp, PLUS); // (Now get phi = M_prec (M_prec^{-1}\phi)

      // Now invert M_prec^{dagger} on it
      QDPIO::cout << "TwoFlavRatioConvConvWilson4DMonomial: resetting Predictor at end of field refresh" << endl;
      getMDSolutionPredictor().reset();
      XMLWriter& xml_out = TheXMLLogWriter::Instance();

      push(xml_out, "FieldRefreshment");
      write(xml_out, "n_count", res.n_count);
      pop(xml_out);

      END_CODE();
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try {
	const TwoFlavorExactRatioConvConvWilsonTypeFermMonomial<P,Q,Phi>& fm = 
	  dynamic_cast<  const TwoFlavorExactRatioConvConvWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactRatioConvConvWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }


      END_CODE();
    }

    //! Reset predictors
    virtual void resetPredictors() {
      getMDSolutionPredictor().reset();

    }

  protected:
    //! Accessor for pseudofermion with Pf
    virtual const Phi& getPhi() const = 0;

    //! Mutator for pseudofermion
    virtual Phi& getPhi() = 0;    

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct() const
      {return getNumerFermAct();}

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action for preconditioner
    virtual const WilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Parameters for inverting with the action of the numerator
    virtual const GroupXML_t& getNumerInvParams() const = 0;

    //! Parameters for inverting with the action of the denominator
    // NB: This is needed, because for some optimized solvers, the ferm act params
    // are actually in part of the solver params. Calling the invMdagM factory
    // function of the Preconditioning Ferm Act may not be sufficient
    // (Also in this case num and denom refer to the final determinant numerator/denominator)
    virtual const GroupXML_t& getDenomInvParams() const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned RatioConvConv type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor unpreconditioned RatioConvConv tpye fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecRatioConvConvWilsonTypeFermMonomial : public TwoFlavorExactRatioConvConvWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecRatioConvConvWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactUnprecRatioConvConvWilsonTypeFermMonomial");

      const WilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();          // for M
      const WilsonTypeFermAct<Phi,P,Q>& FA_prec = getDenomFermAct();  // for M_prec

      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Get the fermion action for the preconditioner
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(state));	
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(FA_prec.linOp(state));

      Phi X;
      
      // Energy calc doesnt use Chrono Predictor
      // Now X here is (M^{dag}M)^{-1} M^{dag}_prec \phi
      //
      // We now need to multiply by M_prec afterwords
      X = zero;
      QDPIO::cout << "TwoFlavRatioConvConvWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();

      // Get system solver
      Handle< MdagMSystemSolver<Phi> > invMdagM(FA.invMdagM(state,getNumerInvParams()));

      // M_dag_prec phi = M^{dag}_prec \phi - the RHS
      Phi M_dag_prec_phi;
      (*M_prec)(M_dag_prec_phi, getPhi(), MINUS);

      // Solve MdagM X = eta
      SystemSolverResults_t res = (*invMdagM)(X, M_dag_prec_phi);

      Phi phi_tmp=zero;
      (*M_prec)(phi_tmp, X, PLUS);

      // Action on the entire lattice
      Double action = innerProductReal(getPhi(), phi_tmp);

      write(xml_out, "n_count", res.n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf
    virtual const Phi& getPhi() const = 0;

    //! mutator for pseudofermion with Pf
    virtual Phi& getPhi() = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const UnprecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get parameters for the inverter
    virtual const GroupXML_t& getNumerInvParams() const = 0;

    //! Get at the chronological predcitor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned RatioConvConv type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned RatioConvConv type fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecRatioConvConvWilsonTypeFermMonomial : public TwoFlavorExactRatioConvConvWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecRatioConvConvWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      // Fermion actions
      const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();
      const WilsonTypeFermAct<Phi,P,Q>& FA_prec = getDenomFermAct();  // for M_prec

      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperator<Phi,P,Q> > M(FA.linOp(state));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(FA_prec.linOp(state));      

      // Get the X fields
      Phi X;

      // Action calc doesnt use chrono predictor use zero guess
      X[ M->subset() ] = zero;

      // Reset chrono predictor
      QDPIO::cout << "TwoFlavRatioConvConvWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();

      // Get system solver
      Handle< MdagMSystemSolver<Phi> > invMdagM(FA.invMdagM(state,getNumerInvParams()));

      // M_dag_prec phi = M^{dag}_prec \phi - the RHS
      Phi M_dag_prec_phi;
      (*M_prec)(M_dag_prec_phi, getPhi(), MINUS);

      // Solve MdagM X = eta
      SystemSolverResults_t res = (*invMdagM)(X, M_dag_prec_phi);

      Phi phi_tmp=zero;
      (*M_prec)(phi_tmp, X, PLUS);

      Double action = innerProductReal(getPhi(), phi_tmp, M->subset());
      
      write(xml_out, "n_count", res.n_count);
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
      push(xml_out, "TwoFlavorExactEvenOddPrecRatioConvConvWilsonTypeFermMonomial");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get parameters for the inverter
    virtual const GroupXML_t& getNumerInvParams() const = 0;

    //! Accessor for pseudofermion with Pf
    virtual const Phi& getPhi() const = 0;

    //! mutator for pseudofermion with Pf 
    virtual Phi& getPhi() = 0;    

    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned RatioConvConv type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned RatioConvConv type fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecConstDetRatioConvConvWilsonTypeFermMonomial : public TwoFlavorExactEvenOddPrecRatioConvConvWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecConstDetRatioConvConvWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    };

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Accessor for pseudofermion with Pf
    virtual const Phi& getPhi() const = 0;
  };

}


#endif
