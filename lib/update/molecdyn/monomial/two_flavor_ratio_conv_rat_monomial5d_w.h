// -*- C++ -*-
// $Id: two_flavor_ratio_conv_rat_monomial5d_w.h,v 3.1 2008-05-23 21:31:34 edwards Exp $

/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __two_flavor_ratio_conv_rat_monomial5d_w_h__
#define __two_flavor_ratio_conv_rat_monomial5d_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/monomial/remez_coeff.h"
#include "update/molecdyn/predictor/chrono_predictor.h"

#include <typeinfo> // For bad_cast

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor RatioConvRat like fermact monomial in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor RatioConvRat like fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   * 
   * Note most of the code is the same here as in the usual 2 flavour
   * case. The only real change is that instead of the PV field, I now
   * use the preconditioning matrix. -- This leads to evil code duplication.
   * Think about abstracting this better.
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D : public ExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*M_2*(M^dag*M)^(-1)*M_2^dag*chi     
      // Here, M is some 5D operator and M_2 is the preconditioning
      // linop, which plays a role identical to the Pauli-Villars matrix
      // in the normal case. The Pauli Villars matrices cancel between 
      // numerator and denominator.
      //
      // Need
      // dS_f/dU =  chi^dag * dM_2 * (M^dag*M)^(-1) * M_2^dag * chi 
      //         -  chi^dag * M_2 * (M^dag*M)^(-1) * [d(M^dag)*M + M^dag*dM] * M_2^dag * (M^dag*M)^(-1) * chi
      //         +  chi^dag * M_2 * (M^dag*M)^(-1) * d(M_2^dag) * chi 
      //
      //         =  chi^dag * dM_2 * psi
      //         -  psi^dag * [d(M^dag)*M + M^dag*dM] * psi
      //         +  psi^dag * d(M_2^dag) * chi 
      //
      // where  psi = (M^dag*M)^(-1) * M_2^dag * chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = chi^dag*\dot(M_2)*X - X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X + X*\dot{M_2}^dag*chi
      // where
      //    X = (M^dag*M)^(-1)*M_2^dag*chi   Y = M*X = (M^dag)^(-1)*M_2^dag*chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      const WilsonTypeFermAct5D<Phi,P,Q>& precFA = getDenomFermAct();

      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Get linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
	
      // Get the mass-style preconditioning linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M_2(precFA.linOp(state));
	
      // Get/construct the pseudofermion solution
      multi1d<Phi> X(FA.size()), Y(FA.size());

      // Move these to get X
      int n_count = this->getX(X,s);

      (*M)(Y, X, PLUS);

      // First M_2 contribution
      M_2->deriv(F, getPhi(), X, PLUS);

      // Last M_2 contribution
      P F_tmp;
      M_2->deriv(F_tmp, X, getPhi(), MINUS);
      F += F_tmp;   // NOTE SIGN

      P FM;
      M->deriv(FM, X, Y, MINUS);
      
      // fold M^dag into X^dag ->  Y  !!
      M->deriv(F_tmp, Y, X, PLUS);
      FM += F_tmp;   // NOTE SIGN

      // Combine forces
      F -= FM;  

      // Recurse only once
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
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      const WilsonTypeFermAct5D<Phi,P,Q>& precFA = getDenomFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(f_state));
      
      // Get the mass-style preconditioning linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M_2(precFA.linOp(f_state));
	
      const int N5 = FA.size();
      multi1d<Phi> eta(N5);
      eta = zero;

      getPhi() = zero;


      // Fill the eta field with gaussian noise
      for(int s=0; s < N5; ++s) {
	gaussian(eta[s], M->subset());
      }

      // Account for fermion BC by modifying the proposed field
      FA.getFermBC().modifyF(eta);
      
      // Temporary: Move to correct normalisation
      for(int s=0; s < N5; ++s)
	eta[s][M->subset()] *= sqrt(0.5);
      
      // Build  phi = M_2 * (M_2^dag*M_2)^(-1) * M^dag * eta
      multi1d<Phi> tmp(N5);

      tmp=zero;

      (*M)(tmp, eta, MINUS);

      // Solve  (V^dag*V)*eta = tmp
      // Get system solver
      Handle< MdagMSystemSolverArray<Phi> > invMdagM(precFA.invMdagM(f_state, getNumerInvParams()));

      // Do the inversion
      SystemSolverResults_t res = (*invMdagM)(eta, tmp);

      // Finally, get phi
      (*M_2)(getPhi(), eta, PLUS);

      // Reset the chronological predictor
      QDPIO::cout << "TwoFlavRatioConvRatWilson5DMonomial: resetting Predictor at end of field refresh" << endl;
      getMDSolutionPredictor().reset();
    
      END_CODE();
    }				    

    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try {
	const TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D<P,Q,Phi>& fm = dynamic_cast< const TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D<P,Q,Phi>& >(m);

	// Do a resize here -- otherwise if the fields have not yet
	// been refreshed there may be trouble
	getPhi().resize(fm.getPhi().size());

	for(int i=0 ; i < fm.getPhi().size(); i++) { 
	  (getPhi())[i] = (fm.getPhi())[i];
	}
      }
      catch(std::bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D" << endl;
	QDP_abort(1);
      }

    
      END_CODE();
    }
  
    //! Reset predictors
    virtual void resetPredictors() {
      getMDSolutionPredictor().reset();

    }

  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi() = 0;    

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const
      {return getNumerFermAct();}

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const =0;

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

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor() = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      const WilsonTypeFermAct5D<Phi,P,Q>& FA_prec = getDenomFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get linop
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
      // Get PV
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M_prec(FA_prec.linOp(state));

      multi1d<Phi> MPrecDagPhi(FA.size());
    
      (*M_prec)(MPrecDagPhi, getPhi(), MINUS);

      // Get system solver
      Handle< MdagMSystemSolverArray<Phi> > invMdagM(FA.invMdagM(state, getNumerInvParams()));

      // CG Chrono predictor needs MdagM
      Handle< DiffLinearOperatorArray<Phi,P,Q> > MdagM(FA.lMdagM(state));
      (getMDSolutionPredictor())(X, *MdagM, MPrecDagPhi);

      // Do the inversion
      SystemSolverResults_t res = (*invMdagM)(X, MPrecDagPhi);

      // Register the new vector
      (getMDSolutionPredictor()).newVector(X);
 
      return res.n_count;
    
      END_CODE();
    }

   };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor unpreconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor unpreconditioned fermact monomial.
   * 
   * CAVEAT: I assume there is only 1 pseudofermion field in the following
   * so called TwoFlavorExact actions.
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial5D : public TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      // SelfEncapsulation/Identification Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactUnprecRatioConvRatWilsonTypeFermMonomial5D");

      // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P,Q>& FA_prec = getDenomFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA_prec.createState(s.getQ()));
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M_prec(FA_prec.linOp(f_state));
 
      multi1d<Phi> X(FA_prec.size());
      multi1d<Phi> tmp(FA_prec.size());

      // Paranoia -- to deal with subsets.
      tmp = zero; 

      // Energy calc does not use chrono predictor
      X = zero;

      // X is now (M^dagM)^{-1} V^{dag} phi

      // getX() now always uses Chrono predictor. Best to Nuke it for
      // energy calcs
      QDPIO::cout << "TwoFlavRatioConvRatWilson5DMonomial: resetting Predictor before energy calc solve" << endl;
      getMDSolutionPredictor().reset();
      int n_count = this->getX(X,s);

      // tmp is now V (M^dag M)^{-1} V^{dag} phi
      (*M_prec)(tmp, X, PLUS);

      // Action on the entire lattice
      Double action = zero;
      for(int s=0; s < FA_prec.size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s]);

      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi() = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor() = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial5D : public TwoFlavorExactRatioConvRatWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& FA_prec = getDenomFermAct();
      
      Handle< FermState<Phi,P,Q> > bc_g_state(FA.createState(s.getQ()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperatorArray<Phi,P,Q> > lin(FA.linOp(bc_g_state));

      Handle< EvenOddPrecLinearOperatorArray<Phi,P,Q> > M_prec(FA_prec.linOp(bc_g_state));
      // Get the X fields
      multi1d<Phi> X(FA.size());

      // X is now (M^dag M)^{-1} V^dag phi

      // Chrono predictor not used in energy calculation
      X = zero;

      // Get X now always uses predictor. Best to nuke it therefore
      QDPIO::cout << "TwoFlavRatioConvRatWilson5DMonomial: resetting Predictor before energy calc solve" << endl;
      getMDSolutionPredictor().reset();
      int n_count = this->getX(X, s);

      multi1d<Phi> tmp(FA.size());
      (*M_prec)(tmp, X, PLUS);

      Double action = zero;
      // Total odd-subset action. NOTE: QDP has norm2(multi1d) but not innerProd
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s], lin->subset());


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
      push(xml_out, "TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial5D");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor() = 0;

    //! Accessor for pseudofermion
    virtual const multi1d<Phi>& getPhi() const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi() = 0;    
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 2 degen flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 2 degen flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class TwoFlavorExactEvenOddPrecConstDetRatioConvRatWilsonTypeFermMonomial5D : public TwoFlavorExactEvenOddPrecRatioConvRatWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecConstDetRatioConvRatWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecConstDetWilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;
    virtual const EvenOddPrecConstDetWilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;
  };



}


#endif
