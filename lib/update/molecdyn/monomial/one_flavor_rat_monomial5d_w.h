// -*- C++ -*-
// $Id: one_flavor_rat_monomial5d_w.h,v 1.1 2005-01-28 02:15:32 edwards Exp $

/*! @file
 * @brief One flavor monomials using RHMC
 */

#ifndef __one_flavor_monomial5d_w_h__
#define __one_flavor_monomial5d_w_h__

#include "update/molecdyn/monomial/abs_monomial.h"

namespace Chroma
{
  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor fermact monomial in extra dimensions
  /*! @ingroup actions
   *
   * Exact 1 flavor fermact monomial. Preconditioning is not
   * specified yet.
   * Can supply a default dsdq and pseudoferm refresh algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatExactWilsonTypeFermMonomial5D : public ExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute dsdq for the system... 
    /*! Actions of the form  chi^dag*(M^dag*M)*chi */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial5D");

      /**** Identical code for unprec and even-odd prec case *****/
      
      // S_f = chi^dag*N(M^dag*M)/D(M^dag*M)*chi  +  chiPV^dag*P(V^dag*V)/Q(V^dag*V)*chiPV
      //
      // Here, M is some 5D operator and V is the Pauli-Villars field
      // N and D makeup the rat. poly of the M term and P and & makeup the rat.poly of the PV term
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
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();
      
      // Create a state for linop
      Handle< const ConnectState> state(FA.createState(s.getQ()));
	
      // Start the force
      int n_count;
      int n_pv_count;
      F = zero;

      // Force term for the linop
      {
	// Get linear operator
	Handle< const DiffLinearOperator<multi1d<Phi>, P> > M(FA.linOp(state));
	
	multi1d< multi1d<Phi> > X;
	multi1d<Phi> Y;

	// Get X out here via multisolver
	n_count = getX(X,getFPartFracRoot(),getPhi(),s);

	// Loop over solns and accumulate force contributions
	P  F_1, F_2;

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
	    F[mu] -= getFPartFracCoeff()[i] * F_1[mu];
	}
      }

      // Force term for the PV
      {
	// Get Pauli-Villars linear operator
	Handle< const DiffLinearOperator<multi1d<Phi>, P> > PV(FA.linOpPV(state));
	
	multi1d< multi1d<Phi> > X;
	multi1d<Phi> Y;

	// Get X out here via multisolver
	n_count = getXPV(X,getFPVPartFracRoot(),getPhiPV(),s);

	// Loop over solns and accumulate force contributions
	P  F_1, F_2;
	
	for(int i=0; i < X.size(); ++i)
	{
	  (*PV)(Y, X[i], PLUS);

	  // The  d(M^dag)*M  term
	  PV->deriv(F_1, X[i], Y, MINUS);
      
	  // The  M^dag*d(M)  term
	  PV->deriv(F_2, Y, X[i], PLUS);
	  F_1 += F_2;

	  // Reweight each contribution in partial fraction
	  for(int mu=0; mu < F.size(); mu++)
	    F[mu] -= getFPVPartFracCoeff()[i] * F_1[mu];
	}
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "n_pv_count", n_pv_count);
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
      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial5DRefresh");

      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< const ConnectState > f_state(FA.createState(s.getQ()));
      
      // Force terms
      const int N5 = FA.size();
      int n_count;
      int n_pv_count;

      // Pseudofermions for M term
      { 
	Handle< const LinearOperator< multi1d<Phi> > > M(FA.linOp(f_state));
      
	multi1d<Phi> eta(N5);
	eta = zero;
      
	// Fill the eta field with gaussian noise
	for(int i=0; i < N5; ++i)
	  gaussian(eta[i], M->subset());
      
	// Temporary: Move to correct normalisation
	for(int i=0; i < N5; ++i)
	  eta[i][M->subset()] *= sqrt(0.5);
      
	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getX(X,getHBPartFracRoot(),eta,s);

	// Weight solns to make final PF field
	if (X.size() != getHBPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	getPhi() = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    getPhi()[j] += getHBPartFracCoeff()[i] * X[i][j];
	}
      }

      // Pseudofermions for PV term
      { 
	Handle< const LinearOperator< multi1d<Phi> > > PV(FA.linOpPV(f_state));
	
	multi1d<Phi> eta(N5);
	eta = zero;
      
	// Fill the eta field with gaussian noise
	for(int i=0; i < N5; ++i)
	  gaussian(eta[i], PV->subset());
      
	// Temporary: Move to correct normalisation
	for(int i=0; i < N5; ++i)
	  eta[i][PV->subset()] *= sqrt(0.5);
      
	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getXPV(X,getHBPVPartFracRoot(),eta,s);

	// Weight solns to make final PF field
	if (X.size() != getHBPVPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	getPhi() = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    getPhi()[j] += getHBPVPartFracCoeff()[i] * X[i][j];
	}
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "n_pv_count", n_pv_count);
      pop(xml_out);
    }

    virtual void setInternalFields(const Monomial<P,Q>& m) 
    {
      try {
	const OneFlavorRatExactWilsonTypeFermMonomial5D<P,Q,Phi>& fm = dynamic_cast< const OneFlavorRatExactWilsonTypeFermMonomial5D<P,Q,Phi>& >(m);

	// Do a resize here -- otherwise if the fields have not yet
	// been refreshed there may be trouble
	getPhi().resize(fm.getPhi().size());
	for(int i=0 ; i < fm.getPhi().size(); i++) { 
	  (getPhi())[i] = (fm.getPhi())[i];
	}

	getPhiPV().resize(fm.getPhiPV().size());
	for(int i=0 ; i < fm.getPhiPV().size(); i++) { 
	  (getPhiPV())[i] = (fm.getPhiPV())[i];
	}
      }
      catch(bad_cast) { 
	QDPIO::cerr << "Failed to cast input Monomial to OneFlavorRatExactWilsonTypeFermMonomial5D" << endl;
	QDP_abort(1);
      }
    }
  

  protected:
    //! Get at fermion action
    virtual const WilsonTypeFermAct5D<Phi,P>& getFermAct(void) const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d<Phi>& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d<Phi>& getPhiPV(void) = 0;    

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracRoot() const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d< multi1d<Phi> >& X, 
		     const multi1d<Real>& shifts, 
		     const multi1d<Phi>& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi for PV using partfrac
    virtual int getXPV(multi1d< multi1d<Phi> >& X, 
		       const multi1d<Real>& shifts, 
		       const multi1d<Phi>& chi, 
		       const AbsFieldState<P,Q>& s) const = 0;
   };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor unpreconditioned fermact monomial living in extra dimensions
  /*! @ingroup actions
   *
   * Exact 1 flavor unpreconditioned fermact monomial.
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatExactUnprecWilsonTypeFermMonomial5D : public OneFlavorRatExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactUnprecWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s) const
    {
      // SelfEncapsulation/Identification Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactUnprecWilsonTypeFermMonomial5D");

      // Get at the fermion action
      const WilsonTypeFermAct5D<Phi,P>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle<const ConnectState> bc_g_state(FA.createState(s.getQ()));

      // Force terms
      const int N5 = FA.size();

      // Action for M term
      Double action_m = zero;
      int n_count;
      {
	Handle< const LinearOperator< multi1d<Phi> > > M(FA.linOp(bc_g_state));

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getX(X,getFPartFracRoot(),getPhi(),s);

	// Weight solns
	if (X.size() != getFPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	multi1d<Phi> tmp(N5);
	tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    tmp[j] += getFPartFracCoeff()[i] * X[i][j];
	}

	// Action on the entire lattice
	for(int i=0; i < N5; ++i)
	  action_m += innerProductReal(getPhi()[i], tmp[i]);
      }

      // Action for PV term
      Double action_pv = zero;
      int n_pv_count;
      {
	Handle< const LinearOperator< multi1d<Phi> > > PV(FA.linOpPV(bc_g_state));

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getX(X,getFPVPartFracRoot(),getPhi(),s);

	// Weight solns
	if (X.size() != getFPVPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	multi1d<Phi> tmp(N5);
	tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    tmp[j] += getFPVPartFracCoeff()[i] * X[i][j];
	}

	// Action on the entire lattice
	for(int i=0; i < N5; ++i)
	  action_pv += innerProductReal(getPhi()[i], tmp[i]);
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "n_pv_count", n_pv_count);
      write(xml_out, "S_m", action_m);
      write(xml_out, "S_PV", action_pv);
      Double action = action_m + action_pv;
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }


  protected:
    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P>& getFermAct(void) const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d<Phi>& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d<Phi>& getPhiPV(void) = 0;    

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracRoot() const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d< multi1d<Phi> >& X, 
		     const multi1d<Real>& shifts, 
		     const multi1d<Phi>& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi for PV using partfrac
    virtual int getXPV(multi1d< multi1d<Phi> >& X, 
		       const multi1d<Real>& shifts, 
		       const multi1d<Phi>& chi, 
		       const AbsFieldState<P,Q>& s) const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup actions
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D : public OneFlavorRatExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) const = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) const
    {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct5D<Phi,P>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle<const ConnectState> bc_g_state(FA.createState(s.getQ()));

      // Force terms
      const int N5 = FA.size();

      // Action for M term
      Double action_m = zero;
      int n_count;
      {
	Handle< const LinearOperator< multi1d<Phi> > > M(FA.linOp(bc_g_state));

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getX(X,getFPartFracRoot(),getPhi(),s);

	// Weight solns
	if (X.size() != getFPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	multi1d<Phi> tmp(N5);
	tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    tmp[j][M->subset()] += getFPartFracCoeff()[i] * X[i][j];
	}

	// Action on the entire lattice
	for(int i=0; i < N5; ++i)
	  action_m += innerProductReal(getPhi()[i], tmp[i], M->subset());
      }

      // Action for PV term
      Double action_pv = zero;
      int n_pv_count;
      {
	Handle< const LinearOperator< multi1d<Phi> > > PV(FA.linOpPV(bc_g_state));

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	n_count = getX(X,getFPVPartFracRoot(),getPhi(),s);

	// Weight solns
	if (X.size() != getFPVPartFracRoot().size())
	  QDP_error_exit("%s : sanity failure, internal size not getHNPartFracRoot size", __func__);

	multi1d<Phi> tmp(N5);
	tmp = zero;
	for(int i=0; i < X.size(); ++i)
	{
	  if (X[i].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	    tmp[j][PV->subset()] += getFPVPartFracCoeff()[i] * X[i][j];
	}

	// Action on the entire lattice
	for(int i=0; i < N5; ++i)
	  action_pv += innerProductReal(getPhi()[i], tmp[i], PV->subset());
      }

      write(xml_out, "n_count", n_count);
      write(xml_out, "n_pv_count", n_pv_count);
      write(xml_out, "S_m", action_m);
      write(xml_out, "S_PV", action_pv);
      Double action = action_m + action_pv;
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s)  const {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D");

      Double action = S_even_even(s) + S_odd_odd(s);

      write(xml_out, "S", action);
      pop(xml_out);
      return action;

    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P>& getFermAct() const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d<Phi>& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d<Phi>& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d<Phi>& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d<Phi>& getPhiPV(void) = 0;    

    //! Return the numerator coefficient in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion
    virtual const multi1d<Real>& getFPartFracRoot() const = 0;

    //! Return the numerator coefficient in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracCoeff() const = 0;

    //! Return the denominator roots in force calc. partial fraction expansion for PV
    virtual const multi1d<Real>& getFPVPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion
    virtual const multi1d<Real>& getHBPartFracRoot() const = 0;

    //! Return the numerator coefficient in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracCoeff() const = 0;

    //! Return the denominator roots in heat-bath partial fraction expansion for PV
    virtual const multi1d<Real>& getHBPVPartFracRoot() const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d< multi1d<Phi> >& X, 
		     const multi1d<Real>& shifts, 
		     const multi1d<Phi>& chi, 
		     const AbsFieldState<P,Q>& s) const = 0;

    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi for PV using partfrac
    virtual int getXPV(multi1d< multi1d<Phi> >& X, 
		       const multi1d<Real>& shifts, 
		       const multi1d<Phi>& chi, 
		       const AbsFieldState<P,Q>& s) const = 0;
  };

}


#endif
