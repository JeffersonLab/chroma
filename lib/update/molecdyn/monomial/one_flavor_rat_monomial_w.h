// -*- C++ -*-
// $Id: one_flavor_rat_monomial_w.h,v 3.12 2009-02-06 15:25:17 bjoo Exp $

/*! @file
 * @brief One flavor monomials using RHMC
 */

#ifndef __one_flavor_rat_monomial_w_h__
#define __one_flavor_rat_monomial_w_h__

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
  class OneFlavorRatExactWilsonTypeFermMonomial : public ExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial");

      /**** Identical code for unprec and even-odd prec case *****/
      
      /* The pseudofermionic action is of the form  S_pf = chi^dag N*D^(-1) chi
       * where N(Q) and D(Q) are polynomials in the rescaled kernel:
       *
       *     Q   = M^dag*M
       *     M   = is some operator
       *
       * The function x^(-alpha) is approximated by c*N(x)/D(x)
       * with c=FRatNorm. This has to be taken into account when comparing
       * the force with the standard HMC routine dsduf for alpha=1.  To solve
       * ( N(Q)/D(Q) )psi the rational function is reformed as a partial fraction
       * expansion and a multishift solver used to find the solution.
       *
       * Need
       * dS_f/dU = - \sum_i psi_i^dag * p_i * [d(M^dag)*M + M^dag*dM] * psi
       *    where    psi_i = (M^dag*M + q_i)^(-1)*chi
       *
       * In Balint's notation, the result is  
       * \dot{S} = -\sum_i p_i [ X_i^dag*\dot{M}^\dag*Y_i - Y_i^dag\dot{M}*X_i]
       * where
       *    X_i = (M^dag*M + q_i)^(-1)*chi   Y_i = M*X_i
       * In Robert's notation,  X_i -> psi_i .
       */
      const WilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.linOp(state));

      // Get multi-shift system solver
      Handle< MdagMMultiSystemSolver<Phi> > invMdagM(FA.mInvMdagM(state, getForceInvParams()));

      // Partial Fraction Expansion coeffs for force
      const RemezCoeff_t& fpfe = getFPFE();

      multi1d<Phi> X;
      Phi Y;

      P  F_1;
      F.resize(Nd);
      F = zero;

      // Loop over all the pseudoferms
      multi1d<int> n_count(getNPF());
      QDPIO::cout << "num_pf = " << getNPF() << endl;

      for(int n=0; n < getNPF(); ++n)
      {
	// The multi-shift inversion
	SystemSolverResults_t res = (*invMdagM)(X, fpfe.pole, getPhi()[n]);
	n_count[n] = res.n_count;

	// Loop over solns and accumulate force contributions


#if 0 

	P F_2;
	P F_tmp(Nd);

	F_tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  (*lin)(Y, X[i], PLUS);

	  // The  d(M^dag)*M  term
	  lin->deriv(F_1, X[i], Y, MINUS);
      
	  // The  M^dag*d(M)  term
	  lin->deriv(F_2, Y, X[i], PLUS);
	  F_1 += F_2;

	  // Reweight each contribution in partial fraction
	  for(int mu=0; mu < F.size(); mu++) {
	    F_tmp[mu] -= fpfe.res[i] * F_1[mu];
	  }
	}
	F += F_tmp;
#else
	// New code with new force term
	multi1d<Phi> Y(X.size());
	for(int i=0; i < X.size(); i++) {
	  (*lin)(Y[i], X[i], PLUS);
	  Y[i]*= -fpfe.res[i];
	}

	lin->derivMultipole(F_1, X, Y, MINUS);
	F += F_1;
	lin->derivMultipole(F_1, Y, X, PLUS);
	F += F_1;
#endif
	// Concious choice. Don't monitor forces by pole 
	// Just accumulate it
	
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
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial");

      // Get at the ferion action for piece i
      const WilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(s.getQ()));
      
      // Create a linear operator
      Handle< DiffLinearOperator<Phi,P,Q> > M(FA.linOp(f_state));
      
      // Get multi-shift system solver
#if 1
      Handle< MdagMMultiSystemSolver<Phi> > invMdagM(FA.mInvMdagM(f_state, getActionInvParams()));
#else
      Handle< MdagMMultiSystemSolverAccumulate<Phi> > invMdagM(FA.mInvMdagMAcc(f_state, getActionInvParams()));
#endif
      // Partial Fraction Expansion coeffs for heat-bath
      const RemezCoeff_t& sipfe = getSIPFE();

      // Loop over pseudoferms
      getPhi().resize(getNPF());
      multi1d<int> n_count(getNPF());
      Phi eta;

      for(int n=0; n < getNPF(); ++n)
      {
	// Fill the eta field with gaussian noise
	eta = zero;
	gaussian(eta, M->subset());
      
	// Account for fermion BC by modifying the proposed field
	FA.getFermBC().modifyF(eta);

	// Temporary: Move to correct normalisation
	eta *= sqrt(0.5);
      
	// The multi-shift inversion
#if 1
	multi1d<Phi> X;
	SystemSolverResults_t res = (*invMdagM)(X, sipfe.pole, eta);
#else
	SystemSolverResults_t res = (*invMdagM)(getPhi()[n], sipfe.norm, sipfe.res,sipfe.pole, eta);
#endif
	n_count[n] = res.n_count;

	// Weight solns to make final PF field
#if 1
	getPhi()[n][M->subset()] = sipfe.norm * eta;
	for(int i=0; i < X.size(); ++i)
	  getPhi()[n][M->subset()] += sipfe.res[i] * X[i];
#endif
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
	const OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to OneFlavorRatExactWilsonTypeFermMonomial " << endl;
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

      const WilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();

      Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< DiffLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
 
      // Get multi-shift system solver
#if 0
      Handle< MdagMMultiSystemSolverAccumulate<Phi> > invMdagM(FA.mInvMdagMAcc(bc_g_state, getActionInvParams()));
#else
     Handle< MdagMMultiSystemSolver<Phi> > invMdagM(FA.mInvMdagM(bc_g_state, getActionInvParams()));
#endif

      // Partial Fraction Expansion coeffs for action
      const RemezCoeff_t& spfe = getSPFE();

      // Compute energy
      // Get X out here via multisolver
#if 1
      multi1d<Phi> X;
#endif
      // Loop over all the pseudoferms
      multi1d<int> n_count(getNPF());
      Double action = zero;
      Phi psi;

      for(int n=0; n < getNPF(); ++n)
      {
#if 0
	// The multi-shift inversion
	SystemSolverResults_t res = (*invMdagM)(psi, spfe.norm, spfe.res,spfe.pole, getPhi()[n]);
#else
	// The multi-shift inversion
	SystemSolverResults_t res = (*invMdagM)(X, spfe.pole, getPhi()[n]);
#endif
	n_count[n] = res.n_count;
	LatticeDouble site_S=zero;

	// Take a volume factor out - redefine zero point energy
	// this constant should have no physical effect, but by
	// making S fluctuate around 0, it should remove a volume
	// factor from the energies
	site_S[ lin->subset() ] = -Double(12);
#if 1	
	psi[lin->subset()] = spfe.norm * getPhi()[n];
	for(int i=0; i < X.size(); ++i) {
	  psi[lin->subset()] += spfe.res[i] * X[i];
	}
#endif
	// Accumulate locally 
	site_S[ lin->subset()] += localNorm2(psi);

	// action += norm2(psi, lin->subset());

	// Sum
	action += sum(site_S, lin->subset());
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }


  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getForceInvParams(void) const = 0;

    //! Return number of roots in used
    virtual int getNPF() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

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
  class OneFlavorRatExactUnprecWilsonTypeFermMonomial : public OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactUnprecWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)
    {
      START_CODE();

      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "OneFlavorRatExactUnprecWilsonTypeFermMonomial");

      Double action = this->S_subset(s);

      write(xml_out, "S", action);
      pop(xml_out);

      END_CODE();

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getForceInvParams(void) const = 0;

    //! Return number of roots in used
    virtual int getNPF() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial
  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial : public OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial() {}

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
      push(xml_out, "OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial");

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
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getActionInvParams(void) const = 0;

    //! Get inverter params
    virtual const GroupXML_t& getForceInvParams(void) const = 0;

    //! Return number of roots in used
    virtual int getNPF() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    
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
  class OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial : 
    public OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }

  protected:
    //!  Replace thiw with PrecConstDet
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;
  };

}


#endif
