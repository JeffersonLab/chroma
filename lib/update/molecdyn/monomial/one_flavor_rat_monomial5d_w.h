// -*- C++ -*-
// $Id: one_flavor_rat_monomial5d_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $

/*! @file
 * @brief One flavor monomials using RHMC
 */

#ifndef __one_flavor_monomial5d_w_h__
#define __one_flavor_monomial5d_w_h__

#include "update/molecdyn/monomial/abs_monomial.h"
#include "actions/ferm/fermacts/remez_coeff.h"
#include "actions/ferm/invert/invert.h"
#include "invtype.h"

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
  class OneFlavorRatExactWilsonTypeFermMonomial5D : public ExactWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... 
    /*! Monomial of the form  chi^dag*(M^dag*M)*chi */
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
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      // Start the force
      F.resize(Nd);
      F = zero;

      // Force term for the linop
      multi1d<int> n_m_count(getNthRoot());
      multi1d<Real> F_m_sq(getNthRoot());
      {
	// Get linear operator
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(state));
	
	// Partial Fraction Expansion coeffs for force
	const RemezCoeff_t& fpfe = getFPFE();

	multi1d< multi1d<Phi> > X;
	multi1d<Phi> Y;
	P  F_1, F_2, F_tmp(Nd);

	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRoot(); ++n)
	{
	  // Get X out here via multisolver
	  n_m_count[n] = getX(X,fpfe.pole,getPhi()[n],s);

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

	  // Take compute the force with respect to thin links
	  // This is wasteful to do here, 
	  // and could be done at the very end, but that would 
	  // screw up the monitoring
	  state->deriv(F_tmp);

	  F_m_sq[n] = norm2(F_tmp);    // monitor force form each nth-root
	  F += F_tmp;                  // add on base force terms
	}
      }

      // Force term for the PV
      multi1d<int> n_pv_count(getNthRootPV());
      multi1d<Real> F_pv_sq(getNthRootPV());
      {
	// Get Pauli-Villars linear operator
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(state));
	
	// Partial Fraction Expansion coeffs for force in PV
	const RemezCoeff_t& fpvpfe = getFPVPFE();

	multi1d< multi1d<Phi> > X;
	multi1d<Phi> Y;
	P  F_1, F_2, F_tmp(Nd);

	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRootPV(); ++n)
	{
	  // Get X out here via multisolver
	  n_pv_count[n] = getXPV(X,fpvpfe.pole,getPhiPV()[n],s);

	  // Loop over solns and accumulate force contributions
	  F_tmp = zero;
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

	      F_tmp[mu] -= fpvpfe.res[i] * F_1[mu];
	  }

	  // Take compute the force with respect to thin links
	  // This is wasteful to do here, 
	  // and could be done at the very end, but that would 
	  // screw up the monitoring
	  state->deriv(F_tmp);


	  F_pv_sq[n] = norm2(F_tmp);    // monitor force form each nth-root
	  F += F_tmp;   // add on PV force term
	}
      }



      Real F_sq = norm2(F);

      write(xml_out, "n_m_count", n_m_count);
      write(xml_out, "n_pv_count", n_pv_count);
      write(xml_out, "F_m_sq", F_m_sq);
      write(xml_out, "F_pv_sq", F_pv_sq);
      write(xml_out, "F_sq", F_sq);
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
      // SelfIdentification/Encapsultaion Rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactWilsonTypeFermMonomial5DRefresh");

      // Heatbath all the fields
      
      // Get at the ferion action for piece i
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();
      
      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > f_state(FA.createState(s.getQ()));
      
      // Force terms
      const int N5 = FA.size();

      // Pseudofermions for M term
      multi1d<int> n_m_count(getNthRoot());
      getPhi().resize(getNthRoot()); // Will hold nth-root pseudoferms
      { 
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(f_state));
      
	// Partial Fraction Expansion coeffs for heat-bath
	const RemezCoeff_t& sipfe = getSIPFE();

	multi1d<Phi> eta(N5);
	
	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRoot(); ++n)
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
      
	  // Get X out here via multisolver
	  multi1d< multi1d<Phi> > X;
	  n_m_count[n] = getX(X,sipfe.pole,eta,s);

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


      // Pseudofermions for PV term
      multi1d<int> n_pv_count(getNthRootPV());
      getPhiPV().resize(getNthRootPV()); // Will hold nth-root pseudoferms
      { 
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(f_state));
	
	// Partial Fraction Expansion coeffs for heat-bath in pv
	const RemezCoeff_t& sipvpfe = getSIPVPFE();

	multi1d<Phi> eta(N5);

	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRootPV(); ++n)
	{
	  // Fill the eta field with gaussian noise
	  eta = zero;
	  for(int i=0; i < N5; ++i)
	    gaussian(eta[i], PV->subset());
      
	  // Temporary: Move to correct normalisation
	  for(int i=0; i < N5; ++i)
	    eta[i][PV->subset()] *= sqrt(0.5);
      
	  // Get X out here via multisolver
	  multi1d< multi1d<Phi> > X;
	  n_pv_count[n] = getXPV(X,sipvpfe.pole,eta,s);

	  // Sanity checks
	  if (X.size() != sipvpfe.pole.size())
	    QDP_error_exit("%s : sanity failure, internal size not getSIPartFracRoot size", __func__);

	  if (X[0].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);
	  
	  // Weight solns to make final PF field
	  getPhiPV()[n].resize(N5);

	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	  {
	    getPhiPV()[n][j][PV->subset()] = sipvpfe.norm * eta[j];
	    for(int i=0; i < X.size(); ++i)
	      getPhiPV()[n][j][PV->subset()] += sipvpfe.res[i] * X[i][j];
	  }
	}
      }

      write(xml_out, "n_m_count", n_m_count);
      write(xml_out, "n_pv_count", n_pv_count);
      pop(xml_out);
    }


    //! Copy internal fields
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
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "S_subset");

      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Create a Connect State, apply fermionic boundaries
      Handle< FermState<Phi,P,Q> > bc_g_state(FA.createState(s.getQ()));

      // Force terms
      const int N5 = FA.size();

      // Action for M term
      Double action_m = zero;
      multi1d<int> n_m_count(getNthRoot());
      {
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > M(FA.linOp(bc_g_state));

	// Partial Fraction Expansion coeffs for action
	const RemezCoeff_t& spfe = getSPFE();

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	multi1d<Phi> tmp(N5);

	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRoot(); ++n)
	{
	  n_m_count[n] = getX(X,spfe.pole,getPhi()[n],s);

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

      // Action for PV term
      Double action_pv = zero;
      multi1d<int> n_pv_count(getNthRootPV());
      {
	Handle< const DiffLinearOperatorArray<Phi,P,Q> > PV(FA.linOpPV(bc_g_state));

	// Partial Fraction Expansion coeffs for action
	const RemezCoeff_t& spvpfe = getSPVPFE();

	// Get X out here via multisolver
	multi1d< multi1d<Phi> > X;
	multi1d<Phi> tmp(N5);

	// Loop over nth-roots, so the pseudoferms
	for(int n=0; n < getNthRootPV(); ++n)
	{
	  n_pv_count[n] = getXPV(X,spvpfe.pole,getPhiPV()[n],s);

	  // Sanity checks
	  if (X.size() != spvpfe.pole.size())
	    QDP_error_exit("%s : sanity failure, internal size not getSPVPartFracRoot size", __func__);

	  if (X[0].size() != N5)
	    QDP_error_exit("%s : sanity failure, internal size not N5", __func__);

	  // Weight solns
	  // Loop over each 5d slice
	  for(int j=0; j < N5; ++j)
	  {
	    tmp[j][PV->subset()] = spvpfe.norm * getPhiPV()[n][j];
	    for(int i=0; i < X.size(); ++i)
	      tmp[j][PV->subset()] += spvpfe.res[i] * X[i][j];
	  }

	  // Action on the subset
	  action_pv += norm2(tmp, PV->subset());
	}
      }

      write(xml_out, "n_m_count", n_m_count);
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
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d< multi1d<Phi> >& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d< multi1d<Phi> >& getPhiPV(void) = 0;    

    //! Return number of roots in used
    virtual int getNthRoot() const = 0;

    //! Return number of roots used in PV
    virtual int getNthRootPV() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

    //! Return the partial fraction expansion for the force calc in PV
    virtual const RemezCoeff_t& getFPVPFE() const = 0;

    //! Return the partial fraction expansion for the action calc in PV
    virtual const RemezCoeff_t& getSPVPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath in PV
    virtual const RemezCoeff_t& getSIPVPFE() const = 0;

    //! Get X = (A^dag*A + q_i)^{-1} eta
    virtual int invert(multi1d< multi1d<Phi> >& X, 
		       const multi1d<Real>& shifts, 
		       const DiffLinearOperatorArray<Phi,P,Q>& A,
		       const multi1d<Phi>& eta) const
    {
      const InvertParam_t& inv_param = getInvParams();

      int n_count = 0;
      multi1d<Real> RsdCG(shifts.size());
      RsdCG = inv_param.RsdCG;
      
      // X allocated and passed in
//    X=zero;

      // Do the inversion...
      switch( inv_param.invType) {
      case CG_INVERTER:
      {
	// Solve A^dag*M X = eta
	MInvCG(A, eta, X, shifts, RsdCG, inv_param.MaxCG, n_count);
	QDPIO::cout << "1Flav5D::invert, n_count = " << n_count << endl;
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


    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi  using partfrac
    virtual int getX(multi1d< multi1d<Phi> >& X, 
		     const multi1d<Real>& shifts, 
		     const multi1d<Phi>& chi, 
		     const AbsFieldState<P,Q>& s) const
    {
      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get linop
      Handle< const DiffLinearOperatorArray<Phi,P,Q> > MdagM(FA.lMdagM(state));

      int n_count = invert(X, shifts, *MdagM, chi);
      return n_count;
    }

  
    //! Multi-mass solver  (M^dagM + q_i)^{-1} chi for PV using partfrac
    virtual int getXPV(multi1d< multi1d<Phi> >& X, 
		       const multi1d<Real>& shifts, 
		       const multi1d<Phi>& chi, 
		       const AbsFieldState<P,Q>& s) const
    {
      // Grab the fermact
      const WilsonTypeFermAct5D<Phi,P,Q>& FA = getFermAct();

      // Make the state
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));

      // Get linop
      Handle< const DiffLinearOperatorArray<Phi,P,Q> > 
	MdagM(new DiffMdagMLinOpArray<Phi,P,Q>(FA.linOpPV(state)));
    
      // Do the inversion...
      int n_count = invert(X, shifts, *MdagM, chi);
      return n_count;
    }

  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor unpreconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
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
    virtual Double S(const AbsFieldState<P,Q>& s) 
    {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactUnprecWilsonTypeFermMonomial5D");

      Double action = S_subset(s);

      write(xml_out, "S", action);
      pop(xml_out);
      return action;
    }


  protected:
    //! Get at fermion action
    virtual const UnprecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d< multi1d<Phi> >& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d< multi1d<Phi> >& getPhiPV(void) = 0;    

    //! Return number of roots in used
    virtual int getNthRoot() const = 0;

    //! Return number of roots used in PV
    virtual int getNthRootPV() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

    //! Return the partial fraction expansion for the force calc in PV
    virtual const RemezCoeff_t& getFPVPFE() const = 0;

    //! Return the partial fraction expansion for the action calc in PV
    virtual const RemezCoeff_t& getSPVPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath in PV
    virtual const RemezCoeff_t& getSIPVPFE() const = 0;
  };


  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
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
    virtual Double S_even_even(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute the odd odd contribution (eg
    virtual Double S_odd_odd(const AbsFieldState<P,Q>& s) 
    {
      return S_subset(s);
    }

    //! Compute the total action
    Double S(const AbsFieldState<P,Q>& s) 
    {
      XMLWriter& xml_out=TheXMLOutputWriter::Instance();
      push(xml_out, "OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D");

      Double action_e = S_even_even(s);
      Double action_o = S_odd_odd(s);
      Double action   = action_e + action_o;

      write(xml_out, "S_even_even", action_e);
      write(xml_out, "S_odd_odd", action_o);
      write(xml_out, "S", action);
      pop(xml_out);

      return action;
    }

  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

    //! Get inverter params
    virtual const InvertParam_t getInvParams(void) const = 0;

    //! Accessor for pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhi(void) const = 0;

    //! mutator for pseudofermion
    virtual multi1d< multi1d<Phi> >& getPhi(void) = 0;    

    //! Accessor for PV pseudofermion (read only)
    virtual const multi1d< multi1d<Phi> >& getPhiPV(void) const = 0;

    //! mutator for PV pseudofermion 
    virtual multi1d< multi1d<Phi> >& getPhiPV(void) = 0;    

    //! Return number of roots in used
    virtual int getNthRoot() const = 0;

    //! Return number of roots used in PV
    virtual int getNthRootPV() const = 0;

    //! Return the partial fraction expansion for the force calc
    virtual const RemezCoeff_t& getFPFE() const = 0;

    //! Return the partial fraction expansion for the action calc
    virtual const RemezCoeff_t& getSPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath
    virtual const RemezCoeff_t& getSIPFE() const = 0;

    //! Return the partial fraction expansion for the force calc in PV
    virtual const RemezCoeff_t& getFPVPFE() const = 0;

    //! Return the partial fraction expansion for the action calc in PV
    virtual const RemezCoeff_t& getSPVPFE() const = 0;

    //! Return the partial fraction expansion for the heat-bath in PV
    virtual const RemezCoeff_t& getSIPVPFE() const = 0;
  };

  //-------------------------------------------------------------------------------------------
  //! Exact 1 flavor even-odd preconditioned fermact monomial living in extra dimensions
  /*! @ingroup monomial
   *
   * Exact 1 flavor even-odd preconditioned fermact monomial.
   * Can supply a default dsdq algorithm
   */
  template<typename P, typename Q, typename Phi>
  class OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial5D : 
    public OneFlavorRatExactEvenOddPrecWilsonTypeFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~OneFlavorRatExactEvenOddPrecConstDetWilsonTypeFermMonomial5D() {}

    //! Even even contribution (eg ln det Clover)
    virtual Double S_even_even(const AbsFieldState<P,Q>& s) {
      return Double(0);
    }


  protected:
    //! Get at fermion action
    virtual const EvenOddPrecWilsonTypeFermAct5D<Phi,P,Q>& getFermAct() const = 0;

  };

}


#endif
