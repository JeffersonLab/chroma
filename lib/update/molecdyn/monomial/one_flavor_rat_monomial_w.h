// -*- C++ -*-
// $Id: one_flavor_rat_monomial_w.h,v 1.1 2005-01-28 02:15:32 edwards Exp $

/*! @file
 * @brief One flavor monomials using RHMC
 */

#ifndef __one_flavor_monomial_w_h__
#define __one_flavor_monomial_w_h__

#include "update/molecdyn/monomial/abs_monomial.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor fermact monomial using rational polynomials
  /*! @ingroup actions
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
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)
    {
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // The pseudofermionic action is of the form  S_pf = chi^dag N*D^(-1) chi
      // where N(Q) and D(Q) are polynomials in the rescaled kernel:
      
      //     Q   = M^dag*M
      //     M   = is some operator

      // The function x^(-alpha) is approximated by c*N(x)/D(x)
      // with c=FRatNorm. This has to be taken into account when comparing
      // the force with the standard HMC routine dsduf for alpha=1.  To solve
      // ( N(Q)/D(Q) )psi the rational function is reformed as a partial fraction
      // expansion and a multishift solver used to find the solution.

      // Need
      // dS_f/dU = - \sum_i psi_i^dag * p_i * [d(M^dag)*M + M^dag*dM] * psi
      //    where    psi_i = (M^dag*M + q_i)^(-1)*chi
      //
      // In Balint's notation, the result is  
      // \dot{S} = -\sum_i p_i [ X_i^dag*\dot{M}^\dag*Y_i - Y_i^dag\dot{M}*X_i]
      // where
      //    X_i = (M^dag*M + q_i)^(-1)*chi   Y_i = M*X_i
      // In Robert's notation,  X_i -> psi_i .
      //
      const WilsonTypeFermAct<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const DiffLinearOperator<Phi,P> > lin(FA.linOp(state));
	
      multi1d<Phi> X;
      Phi Y;

      // Get X out here via multisolver
      int n_count = getX(X,getFPartFracRoot(),getPhi(),s);

      // Loop over solns and accumulate force contributions
      P  F_1, F_2;
      F = zero;

      for(int i=0; i < X.size(); ++i)
      {
	(*lin)(Y, X[i], PLUS);

	// The  d(M^dag)*M  term
	lin->deriv(F_1, X[i], Y, MINUS);
      
	// The  M^dag*d(M)  term
	lin->deriv(F_2, Y, X[i], PLUS);
	F_1 += F_2;

	// Reweight each contribution in partial fraction
	for(int mu=0; mu < F.size(); mu++)
	  F[mu] -= getFPartFracCoeff()[i] * F_1[mu];
      }

      write(xml_out, "n_count", n_count);
      pop(xml_out);
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
      // Heatbath all the fields
      
      // Self Description/Encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial");

      // Get at the ferion action for piece i
      const WilsonTypeFermAct<Phi,P>& S_f = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(S_f.createState(s.getQ()));
      
      // Create a linear operator
      Handle< const LinearOperator<Phi> > M(S_f.linOp(f_state));
      
      Phi eta = zero;
      
      // Fill the eta field with gaussian noise
      gaussian(eta, M->subset());
      
      // Temporary: Move to correct normalisation
      eta *= sqrt(0.5);
      
      // Get X out here via multisolver
      multi1d<Phi> X;
      int n_count = getX(X,getHBPartFracRoot(),eta,s);

      // Weight solns to make final PF field
      getPhi() = zero;
      for(int i=0; i < X.size(); ++i)
	getPhi() += getHBPartFracCoeff()[i] * X[i];

      write(xml_out, "n_count", n_count);
      pop(xml_out);
    }				    
  
    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) {
      try {
	const OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>& fm = dynamic_cast<  const OneFlavorRatExactWilsonTypeFermMonomial<P,Q,Phi>& >(m);

	getPhi() = fm.getPhi();
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to OneFlavorRatExactWilsonTypeFermMonomial " << endl;
	QDP_abort(1);
      }
    }

  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual Phi& getPhi(void) = 0;    

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d<Phi>& X, 
		     const multi1d<Real>& shifts, 
		     const Phi& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor unpreconditioned fermact monomial
  /*! @ingroup actions
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
    /*!
     * This measures the pseudofermion contribution to the Hamiltonian
     * for the case of rational evolution (with polynomials N(x) and D(x),
     * of degree FRatDeg)
     *
     * S_f = chi_dag * N(A)*D(A)^(-1)* chi
     *
     * where A is M^dag*M
     */
    virtual Double S(const AbsFieldState<P,Q>& s) const
    {
      // Self identification/encapsulation Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactUnprecWilsonTypeFermMonomial");

      // Compute energy
      // Get X out here via multisolver
      multi1d<Phi> X;
      int n_count = getX(X,getFPartFracRoot(),getPhi(),s);

      // Weight solns to make final PF field
      Phi psi = zero;
      for(int i=0; i < X.size(); ++i)
	psi += getFPartFracCoeff()[i] * X[i];

      // Action on the entire lattice
      Double action = innerProductReal(getPhi(), psi);
      
      write(xml_out, "n_count", n_count);
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct<Phi,P>& getFermAct(void) const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual Phi& getPhi(void) = 0;    

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d<Phi>& X, 
		     const multi1d<Real>& shifts, 
		     const Phi& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial
  /*! @ingroup actions
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
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) const
    {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct<Phi,P>& FA = getFermAct();

      Handle<const ConnectState> bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< const EvenOddPrecLinearOperator<Phi,P> > lin(FA.linOp(bc_g_state));

      // Compute energy
      // Get X out here via multisolver
      multi1d<Phi> X;
      int n_count = getX(X,getFPartFracRoot(),getPhi(),s);

      // Weight solns to make final PF field
      Phi psi;
      psi[lin->subset()] = zero;
      for(int i=0; i < X.size(); ++i)
	psi[lin->subset()] += getFPartFracCoeff()[i] * X[i];

      Double action = innerProductReal(getPhi(), psi, lin->subset());

      write(xml_out, "n_count", n_count);
      write(xml_out, "S_oo", action);
      pop(xml_out);

      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s) const 
    {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct<Phi,P>& getFermAct() const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const Phi& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual Phi& getPhi(void) = 0;    

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d<Phi>& X, 
		     const multi1d<Real>& shifts, 
		     const Phi& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;
  };

}


#endif
