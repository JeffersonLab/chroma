// -*- C++ -*-
// $Id: one_flavor_ratio_rat_rat_monomial5d_w.h,v 3.1 2008-05-23 21:31:33 edwards Exp $

/*! @file
 * @brief One flavor ratio of rational monomials using RHMC
 */

#ifndef __one_flavor_ratio_rat_rat_monomial5d_w_h__
#define __one_flavor_ratio_rat_rat_monomial5d_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/monomial/remez_coeff.h"

#include <typeinfo>

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor fermact monomial in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 1 flavor fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D : public ExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	// SelfIdentification/Encapsultaion Rule
	XMLWriter& xml_out = TheXMLLogWriter::Instance();
	push(xml_out, "OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D");

	/**** Identical code for unprec and even-odd prec case *****/
      
	// S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi
	//
	// Here, M is some 5D operator and V is the Pauli-Villars field
	// N and D makeup the rat. poly of the M term and P and & makeup the rat.poly of the denom term
	//
	// Need
	// dS_f/dU = - \sum_i psi_i^dag * n_i * [d(M^dag)*M + M^dag*dM] * psi
	//           - \sum_i psiPV_i^dag * p_i * [d(V^dag)*V + V^dag*dV] * psiPV
	//
	//    where    psi_i   = (M^dag*M + n_i)^(-1)*chi
	//    and      psiPV_i = (V^dag*V + p_i)^(-1)*chiPV
	//
	// In Balint's notation, the result is  
	// \dot{S} = -\sum_i p_i [ X_i^dag*\dot{M}^\dag*Y_i - Y_i^dag\dot{M}*X_i]
	// where
	//    X_i = (M^dag*M + q_i)^(-1)*chi   Y_i = M*X_i
	// In Robert's notation,  X_i -> psi_i .
	//
	const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      
	// Create a state for linop
	Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
	// Start the force
	F.resize(Nd);
	F = zero;

	// Force term for the linop
	{
	  // Get multi-shift system solver
	  Handle< MdagMMultiSystemSolverArray<Phi> > invMdagM(FA.mInvMdagM(state, getNumerForceInvParams()));

	  // Get linear operator
	  Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
	
	  // Partial Fraction Expansion coeffs for force
	  const RemezCoeff_t& fpfe = getNumerFPFE();

	  multi1d< multi1d<Phi> > X;
	  multi1d<Phi> Y;
	  P  F_1, F_2, F_tmp(Nd);
	  multi1d<int> n_m_count(getNPF());

	  // Loop over the pseudoferms
	  for(int n=0; n < getNPF(); ++n)
	  {
	    // The multi-shift inversion
	    SystemSolverResults_t res = (*invMdagM)(X, fpfe.pole, getPhi()[n]);
	    n_m_count[n] = res.n_count;

	    // Loop over solns and accumulate force contributions
	    F_tmp = zero;
	    for(int i=0; i < X.size(); ++i)
	    {
	      (*M)(Y, X[i], PLUS);

	      // The  d(M^dag)*M  term
	      M->deriv(F_1, X[i], Y, MINUS);
      
	      // The  M^dag*d(M)  term
	      M->deriv(F_2, Y, X[i], PLUS);
	      F_1 += F_2;

	      // Reweight each contribution in partial fraction
	      for(int mu=0; mu < F.size(); mu++)
		F_tmp[mu] -= fpfe.res[i] * F_1[mu];
	    }

	    // We are not monitoring by pole anymore, so I will
	    // Just accumulate this and take one derivative at the end
	    F += F_tmp;                  // add on base force terms
	  }

	  // Take the derivative with respect to possibly fat links
	  state->deriv(F);

	  // Write out the n_count array and the Forces for the Base Operator
	  write(xml_out, "n_m_count", n_m_count);    
	  monitorForces(xml_out, "ForcesOperator", F);
	}

	monitorForces(xml_out, "Forces", F);

	pop(xml_out);
    
	END_CODE();
      }

  
    //! Refresh pseudofermions
    /*!
     * This routine calculates the pseudofermion field (chi) for the case
     * of rational evolution
     *
     *           chi =  n(Q)*[d(Q)]^(-1) * eta
     * Where:    Q   = M^dag*M
     *           d(Q) = (Q+q_1)*(Q+q_2)*...*(Q+q_m)
     *           n(Q) = (Q+p_1)*(Q+p_2)*...  *(Q+p_m)
     *           m  = SIRatDeg
     *	     							
     * The rational function n(x)/d(x) is the optimal rational
     * approximation to the inverse square root of N(x)/D(x) which in
     * turn is the optimal rational approximation to x^(-alpha).
     * Here, alpha = 1/2
     *
     * To solve {n(Q)*[d(Q)]^(-1) * eta} the partial fraction expansion is
     * used in combination with a multishift solver.
     */
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	// SelfIdentification/Encapsultaion Rule
	XMLWriter& xml_out = TheXMLLogWriter::Instance();
	push(xml_out, "OneFlavorRatioRatRatExactWilsonTypeFermMonomial5DRefresh");

	// Heatbath all the fields
      
	// Get at the ferion action for piece i
	const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();
      
	// Create a Connect State, apply fermionic boundaries
	Handle< FermState<Phi,P,Q> > f_state(FA.createState(s.getQ()));
      
	// Force terms
	const int N5 = FA.size();

	// Pseudofermions for M term
	multi1d<int> n_m_count(getNPF());
	getPhi().resize(getNPF()); // Will hold nth-root pseudoferms
	{ 
	  // Get multi-shift system solver
	  Handle< MdagMMultiSystemSolverArray<Phi> > invMdagM(FA.mInvMdagM(f_state, getNumerActionInvParams()));

	  // Get linear operator
	  Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(f_state));
      
	  // Partial Fraction Expansion coeffs for heat-bath
	  const RemezCoeff_t& sipfe = getNumerSIPFE();

	  multi1d<Phi> eta(N5);
	
	  // Loop over the pseudoferms
	  for(int n=0; n < getNPF(); ++n)
	  {
	    // Fill the eta field with gaussian noise
	    eta = zero;
	    for(int i=0; i < N5; ++i)
	      gaussian(eta[i], M->subset());

	    // Account for fermion BC by modifying the proposed field
	    FA.getFermBC().modifyF(eta);
      
	    // Temporary: Move to correct normalisation
	    for(int i=0; i < N5; ++i)
	      eta[i][M->subset()] *= sqrt(0.5);
      
	    // The multi-shift inversion
	    multi1d< multi1d<Phi> > X;
	    SystemSolverResults_t res = (*invMdagM)(X, sipfe.pole, eta);
	    n_m_count[n] = res.n_count;

	    // Sanity checks
	    if (X.size() != sipfe.pole.size())
	      QDP_error_exit("%s : sanity failure, internal size not getSIPartFracRoot size", __func__);
	  
	    if (X[0].size() != N5)
	      QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	    // Weight solns to make final PF field
	    getPhi()[n].resize(N5);
	    getPhi()[n] = zero;
	    // Loop over each 5d slice
	    for(int j=0; j < N5; ++j)
	    {
	      getPhi()[n][j][M->subset()] = sipfe.norm * eta[j];
	      for(int i=0; i < X.size(); ++i)
		getPhi()[n][j][M->subset()] += sipfe.res[i] * X[i][j];
	    }
	  }
	}

	write(xml_out, "n_m_count", n_m_count);
	pop(xml_out);
    
	END_CODE();
      }


    //! Copy internal fields
    virtual void setInternalFields(const Monomial<P,Q>& m) 
      {
	START_CODE();

	try {
	  const OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D<P,Q,Phi>& fm = dynamic_cast< const OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D<P,Q,Phi>& >(m);

	  // Do a resize here -- otherwise if the fields have not yet
	  // been refreshed there may be trouble
	  getPhi().resize(fm.getPhi().size());
	  for(int i=0 ; i < fm.getPhi().size(); i++) { 
	    (getPhi())[i] = (fm.getPhi())[i];
	  }
	}
	catch(bad_cast) { 
	  QDPIO::cerr << "Failed to cast input Monomial to OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D" << endl;
	  QDP_abort(1);
	}
    
	END_CODE();
      }
  

    //! Compute action on the appropriate subset
    /*! 
     * This may only be a Piece Of The Action
     *
     * This function measures the pseudofermion contribution to the Hamiltonian
     * for the case of rational evolution (with polynomials n(x) and d(x),
     * of degree SRatDeg
     *
     * S_f = chi_dag * (n(A)*d(A)^(-1))^2* chi
     *
     * where A is M^dag*M
     *
     * The rational function n(x)/d(x) is the optimal rational
     * approximation to the square root of N(x)/D(x) which in
     * turn is the optimal rational approximation to x^(-alpha).
     * Here, alpha = 1/2
     */
    virtual Double S_subset(const AbsFieldState<P,Q>& s) const
      {
	START_CODE();

	XMLWriter& xml_out = TheXMLLogWriter::Instance();
	push(xml_out, "S_subset");

	const WilsonTypeFermAct5D<Phi,P,Q>& FA = getNumerFermAct();

	// Create a Connect State, apply fermionic boundaries
	Handle< FermState<Phi,P,Q> > bc_g_state(FA.createState(s.getQ()));

	// Force terms
	const int N5 = FA.size();

	// Action for M term
	Double action_m = zero;
	multi1d<int> n_m_count(getNPF());
	{
	  // Get multi-shift system solver
	  Handle< MdagMMultiSystemSolverArray<Phi> > invMdagM(FA.mInvMdagM(bc_g_state, getNumerActionInvParams()));

	  // Get linear operator
	  Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(bc_g_state));

	  // Partial Fraction Expansion coeffs for action
	  const RemezCoeff_t& spfe = getNumerSPFE();

	  // Get X out here via multisolver
	  multi1d< multi1d<Phi> > X;
	  multi1d<Phi> tmp(N5);

	  // Loop over the pseudoferms
	  for(int n=0; n < getNPF(); ++n)
	  {
	    // The multi-shift inversion
	    SystemSolverResults_t res = (*invMdagM)(X, spfe.pole, getPhi()[n]);
	    n_m_count[n] = res.n_count;

	    // Sanity checks
	    if (X.size() != spfe.pole.size())
	      QDP_error_exit("%s : sanity failure, internal size not getSPartFracRoot size", __func__);
	
	    if (X[0].size() != N5)
	      QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	    // Weight solns
	    // Loop over each 5d slice
	    for(int j=0; j < N5; ++j)
	    {
	      tmp[j][M->subset()] = spfe.norm * getPhi()[n][j];
	      for(int i=0; i < X.size(); ++i)
		tmp[j][M->subset()] += spfe.res[i] * X[i][j];
	    }

	    // Action on the subset
	    action_m += norm2(tmp, M->subset());
	  }
	}

	write(xml_out, "n_m_count", n_m_count);
	write(xml_out, "S_m", action_m);
	Double action = action_m;
	write(xml_out, "S", action);
	pop(xml_out);
    
	END_CODE();

	return action;
      }


  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const
      {return getNumerFermAct();}

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getNumerActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getNumerForceInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getDenomActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getDenomForceInvParams(void) const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getNumerFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getNumerSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getNumerSIPFE() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getDenomFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getDenomSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getDenomSIPFE() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d< multi1d<Phi> >& getPhi(void) = 0;    

    //! Return number of pseudofermions
    virtual int getNPF() const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor unpreconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 1 flavor unpreconditioned fermact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatRatExactUnprecWilsonTypeFermMonomial5D : public OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~OneFlavorRatioRatRatExactUnprecWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	XMLWriter& xml_out=TheXMLLogWriter::Instance();
	push(xml_out, "OneFlavorRatioRatRatExactUnprecWilsonTypeFermMonomial5D");

	Double action = this->S_subset(s);

	write(xml_out, "S", action);
	pop(xml_out);
    
	END_CODE();

	return action;
      }


  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatRatExactEvenOddPrecWilsonTypeFermMonomial5D : public OneFlavorRatioRatRatExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~OneFlavorRatioRatRatExactEvenOddPrecWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) 
      {
	return this->S_subset(s);
      }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	XMLWriter& xml_out=TheXMLLogWriter::Instance();
	push(xml_out, "OneFlavorRatioRatRatExactEvenOddPrecWilsonTypeFermMonomial5D");

	Double action_e = S_even_even(s);
	Double action_o = S_odd_odd(s);
	Double action   = action_e + action_o;

	write(xml_out, "S_even_even", action_e);
	write(xml_out, "S_odd_odd", action_o);
	write(xml_out, "S", action);
	pop(xml_out);

	END_CODE();

	return action;
      }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatRatExactEvenOddPrecConstDetWilsonTypeFermMonomial5D : 
    public OneFlavorRatioRatRatExactEvenOddPrecWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~OneFlavorRatioRatRatExactEvenOddPrecConstDetWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }


  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getDenomFermAct() const = 0;
  };

}

#endif
