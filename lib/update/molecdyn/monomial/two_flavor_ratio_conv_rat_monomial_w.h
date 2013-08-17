// -*- C++ -*-
/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __two_flavor_ratio_conv_rat_monomial_w_h__
#define __two_flavor_ratio_conv_rat_monomial_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/monomial/remez_coeff.h"
#include "update/molecdyn/predictor/chrono_predictor.h"

#include <typeinfo>
using namespace std;

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 2 flavor RatioConvRat type monomial
  /*! @ingroup monomial
   *
   * Exact 2 flavor RatioConvRat type  monomial. 
   *
   * Can supply default dsdq()
   *                    pseudoferm refresh
   *                    getX() algorithm
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactRatioConvRatWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactRatioConvRatWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*M_prec(M^dag*M)^{-1}M^{dag}_prec*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactRatioConvRatWilsonTypeFermMonomial");

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
      const WilsonTypeFermAct<Phi,P,Q>& FAPrec = getDenomFermAct();  // for M_prec

      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.linOp(state));	
      Handle< DiffLinearOperator<Phi,P,Q> > linPrec(FAPrec.linOp(state));

      Phi X=zero;
      Phi Y=zero;

      // Get X out here
      int n_count = this->getX(X,s);
      
      (*lin)(Y, X, PLUS);

      // \phi^{\dagger} \dot(M_prec) X
      linPrec->deriv(F, getPhi(), X, PLUS);
      

      // - X^{\dagger} \dot( M^{\dagger}) Y
      P F_tmp;
      lin->deriv(F_tmp, X, Y, MINUS);
      F -= F_tmp;
 
      // - Y^{\dagger} \dot( M ) X
      lin->deriv(F_tmp, Y, X, PLUS);
      F -= F_tmp;

      // + X^{\dagger} \dot(M_prec)^dagger \phi
      linPrec->deriv(F_tmp, X, getPhi(), MINUS);
      F += F_tmp;

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
      
      // Get at the fermion action for the expensive matrix
      const WilsonTypeFermAct<Phi,P,Q>& S_f = getNumerFermAct();

      // Get the fermion action for the preconditioner
      const WilsonTypeFermAct<Phi,P,Q>& S_prec = getDenomFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(S_f.createState(field_state.getQ()));
      
      // Create a linear operator for the Expensive op
      Handle< DiffLinearOperator<Phi,P,Q> > M(S_f.linOp(f_state));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(S_prec.linOp(f_state));

      Phi eta = zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, M->subset());
      
      // Account for fermion BC by modifying the proposed field
      S_f.getFermBC().modifyF(eta);

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
      Handle< MdagMSystemSolver<Phi> > invMdagM(S_prec.invMdagM(f_state, getNumerInvParams()));

      // Solve MdagM_prec X = eta
      SystemSolverResults_t res = (*invMdagM)(phi_tmp, eta_tmp);

      (*M_prec)(getPhi(), phi_tmp, PLUS); // (Now get phi = M_prec (M_prec^{-1}\phi)

      // Now invert M_prec^{dagger} on it
      QDPIO::cout << "TwoFlavRatioConvRatWilson4DMonomial: resetting Predictor at end of field refresh" << endl;
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
	const TwoFlavorExactRatioConvRatWilsonTypeFermMonomial<P,Q,Phi>& fm = 
	  dynamic_cast<  const TwoFlavorExactRatioConvRatWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactRatioConvRatWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }


      END_CODE();
    }

    //! Reset predictors
    virtual void resetPredictors() {
      getMDSolutionPredictor().reset();

    }

    // We want to generate X = (M^dag M)^{-1} M^{\dagger}_prec \phi
    // Which is a normal solve on M^dag M X = M^{\dagger}_prec \phi
    virtual int getX( Phi& X, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Grab the fermact
      const WilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();
      const WilsonTypeFermAct<Phi,P,Q>& FA_prec = getDenomFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get linop
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(state));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(FA_prec.linOp(state));

      // Get system solver
      const GroupXML_t& inv_param = getNumerInvParams();
      Handle< MdagMSystemSolver<Phi> > invMdagM(FA.invMdagM(state,inv_param));

      SystemSolverResults_t res;

      // Do the inversion...
//    case CG_INVERTER:
      {
	// Solve MdagM X = M^{dag}_prec \phi
	// Do the inversion...

	// Need MdagM for CG based predictor
	Handle< DiffLinearOperator<Phi,P,Q> > MdagM(FA.lMdagM(state));
	Phi M_dag_prec_phi;

	// M_dag_prec phi = M^{dag}_prec \phi - the RHS
	(*M_prec)(M_dag_prec_phi, getPhi(), MINUS);

	// Now done in the solver
	//(getMDSolutionPredictor())(X, *MdagM, M_dag_prec_phi);

	// Solve MdagM X = eta
	res = (*invMdagM)(X, M_dag_prec_phi, getMDSolutionPredictor());
	// Now done in the solver
	// (getMDSolutionPredictor()).newVector(X);
      }

      END_CODE();

      return res.n_count;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi() = 0;    

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct() const
      {return getNumerFermAct();}

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action for preconditioner
    virtual const WilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get parameters for the inverter
    virtual const GroupXML_t& getNumerInvParams() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getDenomActionInvParams() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getDenomForceInvParams() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getDenomFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getDenomSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getDenomSIPFE() const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned RatioConvRat type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor unpreconditioned RatioConvRat tpye fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial : public TwoFlavorExactRatioConvRatWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial");

      Phi X;
      
      // Energy calc doesnt use Chrono Predictor
      // Now X here is (M^{dag}M)^{-1} M^{dag}_prec \phi
      //
      // We now need to multiply by M_prec afterwords
      X = zero;
      QDPIO::cout << "TwoFlavRatioConvRatWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();

      int n_count = this->getX(X,s);


      // Get the fermion action for the preconditioner
      const WilsonTypeFermAct<Phi,P,Q>& S_prec = getDenomFermAct();
      Handle< FermState<Phi,P,Q> > f_state(S_prec.createState(s.getQ()));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(S_prec.linOp(f_state));      

      Phi phi_tmp=zero;
      (*M_prec)(phi_tmp, X, PLUS);

      // Action on the entire lattice
      Double action = innerProductReal(getPhi(), phi_tmp);

      
      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi() = 0;    

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const UnprecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned RatioConvRat type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned RatioConvRat type fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial : public TwoFlavorExactRatioConvRatWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& FA = getNumerFermAct();

      Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      // Get the X fields
      Phi X;

      // Action calc doesnt use chrono predictor use zero guess
      X[ lin->subset() ] = zero;

      // getX noe always uses chrono predictor. Best to Nuke it therefore
      QDPIO::cout << "TwoFlavRatioConvRatWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();

      int n_count = this->getX(X, s);

      const WilsonTypeFermAct<Phi,P,Q>& S_prec = getDenomFermAct();
      Handle< FermState<Phi,P,Q> > f_state(S_prec.createState(s.getQ()));
      Handle< DiffLinearOperator<Phi,P,Q> > M_prec(S_prec.linOp(f_state));      

      Phi phi_tmp=zero;
      (*M_prec)(phi_tmp, X, PLUS);


      Double action = innerProductReal(getPhi(), phi_tmp, lin->subset());
      
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
      push(xml_out, "TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }

  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const Phi& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual Phi& getPhi() = 0;    

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor4D<Phi>& getMDSolutionPredictor() = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned RatioConvRat type fermact monomial
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned RatioConvRat type fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecConstDetRatioConvRatWilsonTypeFermMonomial : public TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecConstDetRatioConvRatWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    };

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at the preconditioned fermion actions
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;
  };

}


#endif
