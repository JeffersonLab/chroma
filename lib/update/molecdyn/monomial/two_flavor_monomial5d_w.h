// -*- C++ -*-
// $Id: two_flavor_monomial5d_w.h,v 3.8 2008-05-23 21:31:34 edwards Exp $

/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __two_flavor_monomial5d_w_h__
#define __two_flavor_monomial5d_w_h__

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
  //! Exact 2 degen flavor fermact monomial in extra dimensions
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s) = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactWilsonTypeFermMonomial5D");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*(M^dag*M)^(-1)*chi     
      // Here, M is some 5D operator
      //
      // Need
      // dS_f/dU = -chi^dag * (M^dag*M)^(-1) * [d(M^dag)*M + M^dag*dM] * (M^dag*M)^(-1) * chi
      //
      // where  psi = (M^dag*M)^(-1) * chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = chi^dag*X - X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X + X*chi
      // where
      //    X = (M^dag*M)^(-1)*chi   Y = M*X = (M^dag)^(-1)*chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Get linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
	
      // Get/construct the pseudofermion solution
      multi1d<Phi> X(FA.size()), Y(FA.size());

      // Move these to get X
      //(getMDSolutionPredictor())(X);
      int n_count = getX(X,s);
      //(getMDSolutionPredictor()).newVector(X);

      (*M)(Y, X, PLUS);

      M->deriv(F, X, Y, MINUS);
      
      // fold M^dag into X^dag ->  Y  !!
      P F_tmp;
      M->deriv(F_tmp, Y, X, PLUS);
      F += F_tmp;
 
      for(int mu=0; mu < F.size(); ++mu) {
	F[mu] *= Real(-1);
      }

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
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(f_state));
      
      const int N5 = FA.size();
      multi1d<Phi> eta(N5);
      eta = zero;
      
      // Fill the eta field with gaussian noise
      for(int s=0; s < N5; ++s)
	gaussian(eta[s], M->subset());

      // Account for fermion BC by modifying the proposed field
      FA.getFermBC().modifyF(eta);
      
      // Temporary: Move to correct normalisation
      for(int s=0; s < N5; ++s)
	eta[s][M->subset()] *= sqrt(0.5);
      
      // Build  phi = V * (V^dag*V)^(-1) * M^dag * eta
      multi1d<Phi> tmp(N5);
      (*M)(tmp, eta, MINUS);

      // Reset the chronological predictor
      QDPIO::cout << "TwoFlavWilson5DMonomial: resetting Predictor at end of field refresh" << endl;
      getMDSolutionPredictor().reset();
    
      END_CODE();
    }				    

    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try 
      {
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
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams() const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor() = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get system solver
      Handle< MdagMSystemSolverArray<Phi> > invMdagM(FA.invMdagM(state, getInvParams()));

      // CG Chrono predictor needs MdagM
      Handle< DiffLinearOperatorArray<Phi,P,Q> > MdagM(FA.lMdagM(state));
      (getMDSolutionPredictor())(X, *MdagM, getPhi());

      // Do the inversion
      SystemSolverResults_t res = (*invMdagM)(X, getPhi());

      // Register the new vector
      (getMDSolutionPredictor()).newVector(X);
 
      END_CODE();

      return res.n_count;
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
  class TwoFlavorExactUnprecWilsonTypeFermMonomial5D : public TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactUnprecWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      // SelfEncapsulation/Identification Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "TwoFlavorExactUnprecWilsonTypeFermMonomial5D");

      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      multi1d<Phi> X(FA.size());

      // Energy calc does not use chrono predictor
      X = zero;

      // getX() now always uses Chrono predictor. Best to Nuke it for
      // energy calcs
      QDPIO::cout << "TwoFlavWilson5DMonomial: resetting Predictor before energy calc solve" << endl;
      getMDSolutionPredictor().reset();
      int n_count = getX(X,s);

      // Action on the entire lattice
      Double action = zero;
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], X[s]);

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
    virtual const UnprecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams() const = 0;

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
  class TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D : public TwoFlavorExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      Handle< FermState<Phi,P,Q> > bc_g_state(FA.createState(s.getQ()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperatorArray<Phi,P,Q> > M(FA.linOp(bc_g_state));

      // Get the X fields
      multi1d<Phi> X(FA.size());

      // Chrono predictor not used in energy calculation
      X = zero;

      // Get X now always uses predictor. Best to nuke it therefore
      QDPIO::cout << "TwoFlavWilson5DMonomial: resetting Predictor before energy calc solve" << endl;
      getMDSolutionPredictor().reset();
      int n_count = getX(X, s);

      Double action = zero;
      // Total odd-subset action. NOTE: QDP has norm2(multi1d) but not innerProd
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], X[s], M->subset());

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
      push(xml_out, "TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
    
      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getInvParams() const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor() = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi() const = 0;

    //! mutator for pseudofermion with Pf index i 
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
  class TwoFlavorExactEvenOddPrecConstDetWilsonTypeFermMonomial5D : public TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~TwoFlavorExactEvenOddPrecConstDetWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecConstDetWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;
  };



}


#endif
