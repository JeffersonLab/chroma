// -*- C++ -*-
// $Id: two_flavor_monomial5d_w.h,v 3.1 2006-06-13 20:18:18 bjoo Exp $

/*! @file
 * @brief Two flavor Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __two_flavor_monomial5d_w_h__
#define __two_flavor_monomial5d_w_h__

#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "actions/ferm/invert/invert.h"
#include "invtype.h"

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
      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactWilsonTypeFermMonomial5D");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*V*(M^dag*M)^(-1)*V^dag*chi     
      // Here, M is some 5D operator and V is the Pauli-Villars field
      //
      // Need
      // dS_f/dU =  chi^dag * dV * (M^dag*M)^(-1) * V^dag * chi 
      //         -  chi^dag * V * (M^dag*M)^(-1) * [d(M^dag)*M + M^dag*dM] * V^dag * (M^dag*M)^(-1) * chi
      //         +  chi^dag * V * (M^dag*M)^(-1) * d(V^dag) * chi 
      //
      //         =  chi^dag * dV * psi
      //         -  psi^dag * [d(M^dag)*M + M^dag*dM] * psi
      //         +  psi^dag * d(V^dag) * chi 
      //
      // where  psi = (M^dag*M)^(-1) * V^dag * chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = chi^dag*\dot(V)*X - X^dag*\dot{M}^\dag*Y - Y^dag\dot{M}*X + X*\dot{V}^dag*chi
      // where
      //    X = (M^dag*M)^(-1)*V^dag*chi   Y = M*X = (M^dag)^(-1)*V^dag*chi
      // In Robert's notation,  X -> psi .
      //
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Get linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
	
      // Get Pauli-Villars linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(state));
	
      // Get/construct the pseudofermion solution
      multi1d<Phi> X(FA.size()), Y(FA.size());

      // Move these to get X
      //(getMDSolutionPredictor())(X);
      int n_count = getX(X,s);
      //(getMDSolutionPredictor()).newVector(X);

      (*M)(Y, X, PLUS);

      // First PV contribution
      PV->deriv(F, getPhi(), X, PLUS);

      // Last PV contribution
      P F_tmp;
      PV->deriv(F_tmp, X, getPhi(), MINUS);
      F += F_tmp;   // NOTE SIGN

      // F now holds PV contribution with respect to fat links
      // Now derive it with respect to the thin links
      // I could do this at the end but this way I can continue
      // monitoring
      state->deriv(F);
      Double F_pv_sq = norm2(F);  // monitor force


      // First interior term
      P FM;
      M->deriv(FM, X, Y, MINUS);
      
      // fold M^dag into X^dag ->  Y  !!
      M->deriv(F_tmp, Y, X, PLUS);
      FM += F_tmp;   // NOTE SIGN

      // Now get derivative with respect to thin links if action is 
      // fat linked
      state->deriv(FM);
      
      Double F_m_sq = norm2(FM);  // monitor force

      F -= FM;  // NOTE SIGN
      Double F_sq = norm2(F);  // monitor force

      write(xml_out, "n_count", n_count);
      write(xml_out, "F_m_sq", F_m_sq);
      write(xml_out, "F_pv_sq", F_pv_sq);
      write(xml_out, "F_sq", F_sq);
      pop(xml_out);
    }
  
    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) 
    {
      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(field_state.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(f_state));
      
      // Get Pauli-Villars linear operator
      Handle< DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(f_state));
	
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

      // Solve  (V^dag*V)*eta = tmp
      int n_pv_count = getXPV(eta, tmp, field_state);

      // Finally, get phi
      (*PV)(getPhi(), eta, PLUS);

      // Reset the chronological predictor
      getMDSolutionPredictor().reset();
    }				    

    virtual void setInternalFields(const Monomial<P,Q>& m) {
      try {
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

      // Reset the chronological predictor
      getMDSolutionPredictor().reset();
    }
  

  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor(void) = 0;

    //! Get (M^dagM)^{-1} phi
    virtual int getX(multi1d<Phi>& X, const AbsFieldState<P,Q>& s)
    {
      const InvertParam_t& inv_param = getInvParams();

      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get linop
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
      // Get PV
      Handle< DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(state));

      multi1d<Phi> VdagPhi(FA.size());
    
      (*PV)(VdagPhi, getPhi(), MINUS);

      int n_count;

      // Get the chrono prediction
      switch( inv_param.invType) {
      case CG_INVERTER:
      {
	// CG Chrono predictor needs MdagM
	Handle< DiffLinearOperatorArray<Phi,P,Q> > MdagM(FA.lMdagM(state));
	(getMDSolutionPredictor())(X, *MdagM, VdagPhi);

	// Do the inversion
	n_count = invert(X, *M, VdagPhi);

	// Register the new vector
	(getMDSolutionPredictor()).newVector(X);
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

  
    //! Get X = (PV^dag*PV)^{-1} eta
    virtual int getXPV(multi1d<Phi>& X, const multi1d<Phi>& eta, const AbsFieldState<P,Q>& s) const
    {
      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
   
      // Get linop
      Handle< DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOpPV(state));

      // Do the inversion...
      int n_count = invert(X, *M, eta);
      return n_count;
    }


    // Get X = (A^dag*A)^{-1} eta
    virtual int invert(multi1d<Phi>& X, const DiffLinearOperatorArray<Phi,P,Q>& M,
		       const multi1d<Phi>& eta) const
    {
      const InvertParam_t& inv_param = getInvParams();

      int n_count =0;

      // Do the inversion...
      switch( inv_param.invType) {
      case CG_INVERTER:
      {
	// Solve MdagM X = eta
	InvCG2(M, eta, X, inv_param.RsdCG, inv_param.MaxCG, n_count);
	QDPIO::cout << "2Flav5D::invert,  n_count = " << n_count << endl;
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
      // SelfEncapsulation/Identification Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactUnprecWilsonTypeFermMonomial5D");

           // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(s.getQ()));
      Handle< DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(f_state));
 
      multi1d<Phi> X(FA.size());
      multi1d<Phi> tmp(FA.size());

      // Paranoia -- to deal with subsets.
      tmp = zero; 

      // Energy calc does not use chrono predictor
      X = zero;

      // X is now (M^dagM)^{-1} V^{dag} phi

      // getX() now always uses Chrono predictor. Best to Nuke it for
      // energy calcs
      getMDSolutionPredictor().reset();
      int n_count = getX(X,s);

      // tmp is now V (M^dag M)^{-1} V^{dag} phi
      (*PV)(tmp, X, PLUS);

      // Action on the entire lattice
      Double action = zero;
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s]);

      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }


  protected:
    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor(void) = 0;
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
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      Handle< FermState<Phi,P,Q> > bc_g_state(FA.createState(s.getQ()));

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperatorArray<Phi,P,Q> > lin(FA.linOp(bc_g_state));

      Handle< EvenOddPrecLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(bc_g_state));
      // Get the X fields
      multi1d<Phi> X(FA.size());

      // X is now (M^dag M)^{-1} V^dag phi

      // Chrono predictor not used in energy calculation
      X = zero;

      // Get X now always uses predictor. Best to nuke it therefore
      getMDSolutionPredictor().reset();
      int n_count = getX(X, s);

      multi1d<Phi> tmp(FA.size());
      (*PV)(tmp, X, PLUS);

      Double action = zero;
      // Total odd-subset action. NOTE: QDP has norm2(multi1d) but not innerProd
      for(int s=0; s < FA.size(); ++s)
	action += innerProductReal(getPhi()[s], tmp[s], lin->subset());


      write(xml_out, "n_count", n_count);
      write(xml_out, "S_oo", action);
      pop(xml_out);

      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)  {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial5D");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
      return action;

    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Get the initial guess predictor
    virtual AbsChronologicalPredictor5D<Phi>& getMDSolutionPredictor(void) = 0;

    //! Accessor for pseudofermion with Pf index i (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion with Pf index i 
    virtual multi1d<Phi>& getPhi(void) = 0;    
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
