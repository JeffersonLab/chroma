// -*- C++ -*-
// $Id: one_flavor_ratio_rat_conv_monomial_w.h,v 3.1 2008-05-23 21:31:33 edwards Exp $

/*! @file
 * @brief One flavor ratio of rational monomials using RHMC
 */

#ifndef __one_flavor_ratio_rat_conv_monomial_w_h__
#define __one_flavor_ratio_rat_conv_monomial_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "eoprec_constdet_wilstype_fermact_w.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "update/molecdyn/monomial/remez_coeff.h"
#include <typeinfo>

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor fermact monomial using rational polynomials
  /*! @ingroup monomial
   *
   * Exact 1 flavor fermact monomial using Rational Polynomial. 
   * Preconditioning is not specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatConvExactWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatioRatConvExactWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*f_2(A_2)*f_1(A_1)*f_2(A_2)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatioRatConvExactWilsonTypeFermMonomial");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*M_2*f_1(A_1)*M_2^dag*chi     
      // Here, A_1=M_1^dag*M_1 is some 4D operator and M_2
      // is the preconditioning linop. The function 
      //    f(A) = norm + \sum_i p_i/(A + q_i)
      //
      // We apply f(A)*chi and define the multi-shift solutions
      //    (A+q_i)\hat{chi}_i = chi
      //    (A+q_i)\hat{phi}_i = phi
      //
      // A core piece of math is that
      // S_sample = chi^dag*f(A)*phi
      //
      // dS_sample/dU 
      //       = chi^dag*df(A)*phi
      //       = - \sum_i p_i * chi^dag* (A+q_i)^{-1} * dA * (A+q_i)^{-1}*phi
      //       = - \sum_i p_i * \hat{chi}_i^dag * dA * \hat{phi}_i
      // so, two solutions are needed
      // Also,   dA = dM^dag*M + M^dag*dM
      //
      // Need
      // dS_f/dU = chi^dag*df_2(A_2)*f_1(A_1)*f_2(A_2)*chi
      //         + chi^dag*f_2(A_2)*f_1(A_1)*df_2(A_2)*chi
      //         + chi^dag*f_2(A_2)*df_1(A_1)*f_2(A_2)*chi
      //
      // Fermion action
      const WilsonTypeFermAct<Phi,P,Q>& FA_num = getNumerFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA_num.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > M_num(FA_num.linOp(state));

      // Get multi-shift system solver
      Handle< MdagMMultiSystemSolver<Phi> > invMdagM_num(FA_num.mInvMdagM(state, getNumerForceInvParams()));

      // Fermion action
      const WilsonTypeFermAct<Phi,P,Q>& FA_den = getDenomFermAct();
      
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > M_den(FA_den.linOp(state));

      // Partial Fraction Expansion coeffs for force
      const RemezCoeff_t& fpfe_num = getNumerFPFE();

      multi1d<Phi> X;
      Phi Y;

      P  F_1, F_2, F_tmp(Nd);
      F.resize(Nd);
      F = zero;

      // Loop over all the pseudoferms
      multi1d<int> n_count(getNPF());

      for(int n=0; n < getNPF(); ++n)
      {
	// M_dag_den phi = M^{dag}_den \phi - the RHS
	Phi M_dag_den_phi;
	(*M_den)(M_dag_den_phi, getPhi()[n], MINUS);

	// The multi-shift inversion
	SystemSolverResults_t res = (*invMdagM_num)(X, fpfe_num.pole, M_dag_den_phi);
	n_count[n] = res.n_count;

	// Loop over solns and accumulate force contributions
	F_tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  (*M_num)(Y, X[i], PLUS);

	  // The  d(M^dag)*M  term
	  M_num->deriv(F_1, X[i], Y, MINUS);
      
	  // The  M^dag*d(M)  term
	  M_num->deriv(F_2, Y, X[i], PLUS);
	  F_1 += F_2;

	  // Reweight each contribution in partial fraction
	  for(int mu=0; mu < F.size(); mu++)
	    F_tmp[mu] -= fpfe_num.res[i] * F_1[mu];
	}

	// Concious choice. Don't monitor forces by pole 
	// Just accumulate it
	F += F_tmp;
      }

      state->deriv(F);
      write(xml_out, "n_count", n_count);
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
     *           m  = HBRatDeg
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

      // Heatbath all the fields
      
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatioRatConvExactWilsonTypeFermMonomial");

      // Get at the fermion action for piece i
      const WilsonTypeFermAct<Phi,P,Q>& FA_num = getNumerFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > state(FA_num.createState(s.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > M_num(FA_num.linOp(state));
      
      // Get multi-shift system solver
      Handle< MdagMMultiSystemSolver<Phi> > invMdagM_num(FA_num.mInvMdagM(state, getNumerActionInvParams()));

      // Get at the fermion action for piece i
      const WilsonTypeFermAct<Phi,P,Q>& FA_den = getDenomFermAct();
      
      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > M_den(FA_den.linOp(state));
      
      // Partial Fraction Expansion coeffs for heat-bath
      const RemezCoeff_t& sipfe = getNumerSIPFE();

      // Loop over pseudoferms
      getPhi().resize(getNPF());
      multi1d<int> n_count(getNPF());
      Phi eta;

      for(int n=0; n < getNPF(); ++n)
      {
	// Fill the eta field with gaussian noise
	eta = zero;
	gaussian(eta, M_num->subset());
      
	// Account for fermion BC by modifying the proposed field
	FA_num.getFermBC().modifyF(eta);

	// Temporary: Move to correct normalisation
	eta *= sqrt(0.5);
      
	// The multi-shift inversion
	multi1d<Phi> X;
	SystemSolverResults_t res = (*invMdagM_num)(X, sipfe.pole, eta);
	n_count[n] = res.n_count;

	// Weight solns to make final PF field
	getPhi()[n][M_num->subset()] = sipfe.norm * eta;
	for(int i=0; i < X.size(); ++i)
	  getPhi()[n][M_num->subset()] += sipfe.res[i] * X[i];
      }

      write(xml_out, "n_count", n_count);
      pop(xml_out);

      END_CODE();
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      START_CODE();

      try {
	const OneFlavorRatioRatConvExactWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const OneFlavorRatioRatConvExactWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to OneFlavorRatioRatConvExactWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }

      END_CODE();
    }


    //! Compute the action on the appropriate subset
    /*!
     * This measures the pseudofermion contribution to the Hamiltonian
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

      const WilsonTypeFermAct<Phi,P,Q>& FA_num = getNumerFermAct();

      Handle< FermState<Phi,P,Q> > state = FA_num.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > M_num(FA_num.linOp(state));
 
      // Get multi-shift system solver
      Handle< MdagMMultiSystemSolver<Phi> > invMdagM_num(FA_num.mInvMdagM(state, getNumerActionInvParams()));

      // Fermion action
      const WilsonTypeFermAct<Phi,P,Q>& FA_den = getDenomFermAct();
      
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > M_den(FA_den.linOp(state));

      // Partial Fraction Expansion coeffs for action
      const RemezCoeff_t& spfe = getNumerSPFE();

      // Compute energy
      // Get X out here via multisolver
      multi1d<Phi> X;

      // Loop over all the pseudoferms
      multi1d<int> n_count(getNPF());
      Double action = zero;
      Phi psi;

      for(int n=0; n < getNPF(); ++n)
      {
	// The multi-shift inversion
	SystemSolverResults_t res = (*invMdagM_num)(X, spfe.pole, getPhi()[n]);
	n_count[n] = res.n_count;

	// Weight solns to make final PF field
	psi[M_num->subset()] = spfe.norm * getPhi()[n];
	for(int i=0; i < X.size(); ++i)
	  psi[M_num->subset()] += spfe.res[i] * X[i];

	action += norm2(psi, M_num->subset());
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }


  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const
      {return getNumerFermAct();}

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getNumerFermAct(void) const = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getDenomFermAct(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getNumerActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getNumerForceInvParams(void) const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getNumerFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getNumerSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getNumerSIPFE() const = 0;

    //! Return number of roots in used
    virtual int getNPF() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    

  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor unpreconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 1 flavor unpreconditioned fermact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatConvExactUnprecWilsonTypeFermMonomial : public OneFlavorRatioRatConvExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatioRatConvExactUnprecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatioRatConvExactUnprecWilsonTypeFermMonomial");

      Double action = this->S_subset(s);

      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getNumerFermAct(void) const = 0;

    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getDenomFermAct(void) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatConvExactEvenOddPrecWilsonTypeFermMonomial : public OneFlavorRatioRatConvExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatioRatConvExactEvenOddPrecWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s)
    {
      return this->S_subset(s);
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatioRatConvExactEvenOddPrecWilsonTypeFermMonomial");

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
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial constant determinant
  //  Can fill out the S_odd_odd piece

  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatioRatConvExactEvenOddPrecConstDetWilsonTypeFermMonomial : 
    public OneFlavorRatioRatConvExactEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatioRatConvExactEvenOddPrecConstDetWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getNumerFermAct() const = 0;

    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getDenomFermAct() const = 0;
  };

}


#endif
