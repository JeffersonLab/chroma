// $Id: overlap_fermact_base_w.cc,v 3.4 2007-02-22 21:11:45 bjoo Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#include "chromabase.h"
//#include "actions/ferm/fermstates/overlap_state.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/invert/invcg1.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/invert/inv_rel_cg1.h"
#include "actions/ferm/invert/inv_rel_cg2.h"

#include "actions/ferm/invert/invsumr.h"
#include "actions/ferm/invert/inv_rel_sumr.h"
#include "actions/ferm/invert/minvsumr.h"
#include "actions/ferm/invert/minv_rel_sumr.h"
#include "actions/ferm/invert/minvcg.h"
#include "actions/ferm/invert/minv_rel_cg.h"
#include "actions/ferm/invert/inv_rel_gmresr_sumr.h"
#include "actions/ferm/invert/inv_rel_gmresr_cg.h"
#include "actions/ferm/linop/lopscl.h"
#include "meas/eig/ischiral_w.h"

#include "actions/ferm/invert/syssolver_cg_params.h"
#include "actions/ferm/invert/multi_syssolver_cg_params.h"

namespace Chroma 
{ 
  //! Propagator for unpreconditioned overlap-like fermion actions
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for DWF-like fermions
   */
  class Ovlap4DQprop : public SystemSolver<LatticeFermion>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Constructor
    /*!
     * \param S_f_       Fermion action ( Read )
     * \param state_     state ( Read )
     */
    Ovlap4DQprop(Handle< OverlapFermActBase> S_f_,
		 Handle< FermState<T,P,Q> > state_,
		 const SysSolverCGParams& invParam_) : 
      S_f(S_f_), state(state_), invParam(invParam_) {}
    
    //! Destructor is automatic
    ~Ovlap4DQprop() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (LatticeFermion& psi, const LatticeFermion& chi) const
      {
	START_CODE();

	SystemSolverResults_t res;
	Real mass = S_f->getQuarkMass();

//	if( invType == "CG_INVERTER") 
	{
	  // We do our solve into tmp, so make sure this is the supplied initial guess
	  LatticeFermion tmp = psi;

	  Handle< LinearOperator<T> > M(S_f->linOp(state));
    
	  // Check whether the source is chiral.
	  Chirality ichiral = isChiralVector(chi);
	  if( ichiral == CH_NONE || ( S_f->isChiral() == false )) 
	  { 
	    Handle<LinearOperator<T> > MM(S_f->lMdagM(state));

	    // Do this at the end as otherwise it may mix chiralities?
	    (*M)(tmp, chi, MINUS);
    

	    // Source is not chiral. In this case we should use,
	    // InvCG2 with M
	    res = InvCG2(*M, tmp, psi, invParam.RsdCG, invParam.MaxCG);

	  }
	  else 
	  {
	    // If we have a chiral source we have M^{dag} M = M^{2}
	    // but applying M at the front might mix chiralities hurting
	    // convergence so we do it last
	    Handle<LinearOperator<T> > MM(S_f->lMdagM(state,ichiral));

	    // Source is chiral. In this case we should use InvCG1
	    // with the special MdagM
	    res = InvCG1(*MM, chi,tmp, invParam.RsdCG, invParam.MaxCG);
	    (*M)(psi,tmp, MINUS);
	  }
  
	  LatticeFermion Mpsi;
	  (*M)(Mpsi, psi, PLUS);
	  Mpsi -= chi;
    
	  QDPIO::cout << "OvQprop || chi - D psi ||/||chi|| = " 
		      << sqrt(norm2(Mpsi))/sqrt(norm2(chi)) 
		      << "  n_count = " << res.n_count << " iters" << endl;
	}
#if 0
	else if (invType == "REL_CG_INVERTER")
	{
	  LatticeFermion tmp = psi;
	  Handle<LinearOperator<T> > M(S_f->linOp(state));


	  // Check whether the source is chiral.
	  Chirality ichiral = isChiralVector(chi);
	  if( ichiral == CH_NONE || ( S_f->isChiral() == false )) 
	  { 
	    
	    // Source is not chiral. In this case we should use,
	    // InvCG2 with M
	    (*M)(tmp, chi, MINUS); 
	    InvRelCG2(*M, tmp, psi, invParam.RsdCG, invParam.MaxCG, res.n_count);
	  }
	  else 
	  {
	    // Source is chiral. In this case we should use InvCG1
	    // with the special MdagM
	    Handle<LinearOperator<T> > MM( dynamic_cast< LinearOperator<LatticeFermion>* >( S_f->lMdagM(state, ichiral) ) );

	    (*M)(tmp, chi, MINUS);	
	    InvRelCG1(*MM, tmp, psi, invParam.RsdCG, invParam.MaxCG, res.n_count);
	  }

	  LatticeFermion Mpsi;
	  (*M)(Mpsi, psi, PLUS);
	  Mpsi -= chi;
	  QDPIO::cout << "OvQprop || chi - D psi ||/||chi|| = " 
		      << sqrt(norm2(Mpsi))/sqrt(norm2(chi))
		      << "  n_count = " << res.n_count << " iters" << endl;

	}
#if 0
	else if (invType == "MR_INVERTER")
	{
	  // psi = D^(-1)* chi
	  InvMR (*M, chi, psi, MRover, invParam.RsdCG, invParam.MaxCG, res.n_count);
	}
	else if (invType == "BICG_INVERTER")
	{
	  // psi = D^(-1) chi
	  InvBiCG (*M, chi, psi, invParam.RsdCG, invParam.MaxCG, res.n_count);
	}
#endif
	else if (invType == "SUMR_INVERTER")
	{
	  // Solve by SUMR solver -- for shifted unitary matrices
	  //
	  // Solve zeta I + rho gamma_5 eps(H)
	  // where gamma_5 eps(H) is unitary
	  //
	  // zeta = (1 + mu)/(1-mu)
	  // rho  = 1
	  Real rho = Real(1);
	  Real mu = S_f->getQuarkMass();
	  Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
	  {
	    Handle< LinearOperator<T> > U(S_f->lgamma5epsH(state));
	
	    // Now solve:
	    InvSUMR(*U, chi, psi, zeta, rho, invParam.RsdCG, invParam.MaxCG, res.n_count);
	
	    // Restore to normal scaling
	    Real fact = Real(2)/(Real(1) - mu);
	    psi *= fact;
	  }

	  // Check back:
	  // Get a proper operator and compute chi- Dpsi
	  Handle< LinearOperator<T> > D(S_f->linOp(state));
	  LatticeFermion Dpsi;
	  (*D)(Dpsi, psi, PLUS);
	  Dpsi -= chi;
	  QDPIO::cout << "OvQprop || chi - D psi || = " 
		      << sqrt(norm2(Dpsi))/sqrt(norm2(chi))
		      << "  n_count = " << res.n_count << " iters" << endl;
	}
	else if (invType == "REL_SUMR_INVERTER")
	{
	  // Solve by Relaxed SUMR solver -- for shifted unitary matrices
	  //
	  // Solve zeta I + rho gamma_5 eps(H)
	  // where gamma_5 eps(H) is unitary
	  //
	  // zeta = (1 + mu)/(1-mu)
	  // rho  = 1
	  Real rho = Real(1);
	  Real mu = S_f->getQuarkMass();
	  Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
	  {
	    Handle< LinearOperator<T> > U( dynamic_cast< LinearOperator<LatticeFermion>* >(S_f->lgamma5epsH(state)) );

	    Real fact = Real(2)/(Real(1) - mu);
	
	    // Now solve:
	    InvRelSUMR(*U, chi, psi, zeta, rho, invParam.RsdCG, invParam.MaxCG, res.n_count);
	
	    // Restore to normal scaling
	    psi *= fact;
	
	
	  }

	  // Check back:
	  // Get a proper operator and compute chi- Dpsi
	  Handle< LinearOperator<T> > D(S_f->linOp(state));
	  LatticeFermion Dpsi;
	  (*D)(Dpsi, psi, PLUS);
	  Dpsi -= chi;
	  QDPIO::cout << "OvQprop || chi - D psi || = " 
		      << sqrt(norm2(Dpsi)) / sqrt(norm2(chi)) 
		      << "  n_count = " << res.n_count << " iters" << endl;
	}
	else if (invType == "REL_GMRESR_SUMR_INVERTER")
	{
	  // Solve by Relaxed SUMR solver -- for shifted unitary matrices
	  //
	  // Solve zeta I + rho gamma_5 eps(H)
	  // where gamma_5 eps(H) is unitary
	  //
	  // zeta = (1 + mu)/(1-mu)
	  // rho  = 1
	  Real rho = Real(1);
	  Real mu = S_f->getQuarkMass();
	  Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
	  {
	    Handle< LinearOperator<T> > UnprecU( dynamic_cast< LinearOperator<T>* >(S_f->lgamma5epsH(state)) );

	    Handle< LinearOperator<T> > PrecU( dynamic_cast< LinearOperator<T>* >(S_f->lgamma5epsHPrecondition(state)) );

	
	    Real fact = Real(2)/(Real(1) - mu);
	
	    // Now solve:
	    InvRelGMRESR_SUMR(*PrecU, zeta, rho, *UnprecU, chi, psi, invParam.RsdCG, invParam.RsdCGPrec, invParam.MaxCG, invParam.MaxCGPrec, res.n_count);
	
	    // Restore to normal scaling
	    psi *= fact;
	
	
	  }

	  // Check back:
	  // Get a proper operator and compute chi- Dpsi
	  Handle< LinearOperator<T> > D(S_f->linOp(state));
	  LatticeFermion Dpsi;
	  (*D)(Dpsi, psi, PLUS);
	  Dpsi -= chi;

	  QDPIO::cout << "OvQprop || chi - D psi || = " 
		      << sqrt(norm2(Dpsi))/sqrt(norm2(chi))
		      << "  n_count = " << res.n_count << " iters" << endl;
	}
	else if (invType == "REL_GMRESR_GG_INVERTER")
	{
      
	  LatticeFermion tmp= psi;
	  Chirality ichiral = isChiralVector(chi);
	  LinearOperator<T> *MM_ptr;
	  LinearOperator<T> *MM_prec_ptr;

	  if( ichiral == CH_NONE || ( S_f->isChiral() == false )) { 
	
	    MM_ptr =  dynamic_cast< LinearOperator<T>* >( S_f->lMdagM(state));
	    MM_prec_ptr =  dynamic_cast< LinearOperator<T>* >( S_f->lMdagMPrecondition(state));
	  
	  }
	  else {
	
	    // Source is chiral. In this case we should use InvCG1
	    // with the special MdagM
	    MM_ptr = dynamic_cast< LinearOperator<T>* >( S_f->lMdagM(state, ichiral) );
	    MM_prec_ptr = dynamic_cast< LinearOperator<T>* >( S_f->lMdagMPrecondition(state, ichiral) );

	  }

	  Handle< LinearOperator<T> > MM(MM_ptr);
	  Handle< LinearOperator<T> > MM_prec(MM_prec_ptr);
	  Handle< LinearOperator<T> > M(S_f->linOp(state));

	  (*M)(tmp,chi,MINUS);
	  InvRelGMRESR_CG(*MM_prec,*MM, tmp, psi, invParam.RsdCG, invParam.RsdCGPrec, invParam.MaxCG, invParam.MaxCGPrec, res.n_count);
	
   

	  T Mpsi;
	  (*M)(Mpsi, psi, PLUS);
	  Mpsi -= chi;
	  QDPIO::cout << "OvQprop || chi - D psi ||/||chi|| = "
		      << sqrt(norm2(Mpsi)) / sqrt(norm2(chi))
		      << "  n_count = " << res.n_count << " iters" << endl;
	}
	else
	{
	  QDPIO::cerr << "Zolotarev4DFermActBj::qprop Solver Type not implemented" << endl;
	  QDP_abort(1);
	}
#endif

	if ( res.n_count == invParam.MaxCG ) { 
	  QDP_error_exit("Zolotarev4DFermAct::qprop: No convergence in solver: n_count = %d\n", res.n_count);
	}

	// Compute residual
	{
	  Handle< LinearOperator<T> > M(S_f->linOp(state));
    
	  T  r;
	  (*M)(r, psi, PLUS);
	  r -= chi;
	  res.resid = sqrt(norm2(r));
	}

	// Normalize and remove contact term 
	Real ftmp = Real(1) / ( Real(1) - mass );
  
	psi -= chi;
	psi *= ftmp;

	END_CODE();

	return res;
      }

  private:
    // Hide default constructor
    Ovlap4DQprop() {}

    Handle< OverlapFermActBase > S_f;
    Handle< FermState<T,P,Q> > state;
    const SysSolverCGParams invParam;
  };

 
// Propagator for unpreconditioned overlap-like fermion actions
/* Yuk, just make a clone of the current action and pass it around */
  SystemSolver<LatticeFermion>* 
  OverlapFermActBase::qprop(Handle< FermState<T,P,Q> > state,
			    const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return new Ovlap4DQprop(Handle<OverlapFermActBase>(clone()), state, 
			    SysSolverCGParams(paramtop,invParam.path));
  }



/* This routine is Wilson type Overlap fermions */

/* Compute multiple quark mass propagators for an unpreconditioned overla fermion */
/* using the single mass source in "chi" - so, the source can */
/* be of any desired form. The results will appear in "psi", which is */
/* ignored on input (and set to zero initially) */

/* This routine can return the solution to the following systems,
   when using CG */
/*   D * psi           = chi   (n_soln = 1) */
/*   D^dag * psi       = chi   (n_soln = 2) */
/*  and all of the above for n_soln = 3 */
/*   (D * D^dag) * psi = chi   (for n_soln = 4) */
 
/* u        -- gauge field ( Read ) */
/* chi      -- source ( Modify ) */
/* Mass     -- quark masses in lattice units ( Read ) */
/* psi      -- quark propagators ( Write ) */
/* ncg_had  -- number of CG iterations ( Modify ) */
  void 
  OverlapFermActBase::multiQprop(multi1d<T>& psi,
				 const multi1d<Real>& masses,
				 Handle< FermState<T,P,Q> > state, 
				 const T& chi, 
				 const GroupXML_t& invParam_,
				 const int n_soln,
				 int& n_count) const
  {

    if ( toBool( getQuarkMass() != 0 ) ) {
      QDP_error_exit("Multi Mass Qprop only works if action has mass = 0 (strictly)");
    }

    std::istringstream  is(invParam_.xml);
    XMLReader  paramtop(is);
	
    MultiSysSolverCGParams invParam(paramtop, invParam_.path);

    int nsets ;
    int msets ;
    multi1d<enum PlusMinus> isign_set(2);


    switch(n_soln) {
    case 1:             //  D psi = chi
      nsets = 1;
      isign_set[0] = MINUS;
      break;

    case 2:
      nsets = 1;         // D^dag psi = chi
      isign_set[0] = PLUS;
      break;

    case 3:              //  both  D psi = chi and D dag psi = chi
      nsets = 2;
      isign_set[0] = MINUS;
      isign_set[1] = PLUS;
      break;

    case 4:
      nsets = 0;           // D^dag D psi = chi
      isign_set[0] = PLUS;
      break;

    default:
      QDP_error_exit("solution type incorrect", n_soln);
    }
  

    const int n_mass = masses.size();

    // Generic feature of multi mass solvers, initial guesses psi have to be zero
    for(int n=0; n < n_mass; n++) { 
      psi[n] = zero;
    }
  
  
    Real ftmp;
//    switch( invParam.invType ) {

//    case CG_INVERTER: 
    {
      // This is M_scaled = 2 D(0). 
      // the fact that D has zero masss is enforced earlier 
      //
      Handle< LinearOperator<T> > 
	M_scaled(new lopscl<T, Real>( linOp(state), Real(2) ) );

      Chirality ischiral;
      LinearOperator<T>* MdagMPtr;
      multi1d<Real> shifted_masses(n_mass);
    
      for(int i = 0; i < n_mass; i++) { 
	shifted_masses[i] = ( Real(1) + masses[i] )/( Real(1) - masses[i] );
	shifted_masses[i] += ( Real(1) / shifted_masses[i] );
	shifted_masses[i] -= Real(2);
      }

      
      // This is icky. I have to new and delete, 
      // and I can't drop in a handle as I
      //    a) Can't declare variables in a case statement
      //       in case it jumps.
      // 
      // Hence I have to deal directly with the f***ing pointer
      // using new and delete.
      //
      //
      // MdagMPtr = 4 D^{dag}(0) D(0)
      //
      // The Zero is enforced elsewhere
      ischiral = isChiralVector(chi);
      if( ischiral == CH_NONE || ( isChiral() == false ) ) {

	// I have one scaled by 2. The MdagM fixes up the scale
	// factor of 4. I don't need to recompute coeffs etc here.
	//
      
	MdagMPtr = new MdagMLinOp<T>(M_scaled);
      }
      else { 
	// Special lovddag version
	MdagMPtr = new lopscl<T, Real> ( lMdagM( state, ischiral ), Real(4) );
      }
      
      // Do the solve
      //
      // This really solves with Kerel:
      //
      // (1/(1-m^2)) 4 D(m)^{dag} D(m)
      //  
      //  =  4 D(0)^{dag} D(0) + shift
      //
      // with shift = (1+m)/(1-m) + (1-m)/(1+m) - 2
      // 
      //
      // result of  solution is:
      //
      // psi = [ (1 - m^2 ) / 4 ] [ D^{dag}(m) D(m) ]^{-1} chi
      //
      //     = [ ( 1 - m^2 ) / 2 ] (2 D^{dag}(m))^{-1} D(m)^{-1} chi
      //
      //     = [ ( 1 - m^2 ) / 2 ] ( 2 D(m))^{-1} D^{-dag}(m) chi
      MInvCG(*MdagMPtr, chi, psi, shifted_masses, invParam.RsdCG, invParam.MaxCG, n_count);

      // Delete MdagMPtr. Do it right away before I forget.
      delete MdagMPtr;

      if ( n_count == invParam.MaxCG ) {
	QDP_error_exit("no convergence in the inverter", n_count);    
      }
      
      // Now make the solution(s)
      switch(n_soln) { 
      case 4:

	// Compensate for the  4/(1-m^2) factor

	for(int i = 0; i < n_mass; i++) { 
	  psi[i] *= Real(4) / ( Real(1) - masses[i]*masses[i] );
	}
	break;
	
      case 3:
	// Copy ofver the solutions of D^{-1} for solutions to D_dag^{-1}
	for( int i = 0 ; i < n_mass; i++) {
	  psi[n_mass + i] = psi[i];
	}
	// Fall through to the case below:

      case 1:
      case 2:
	for(int iset=0; iset < nsets; iset++) {
	  for(int i=0; i < n_mass; i++) {
	    int j = i+iset*n_mass;
	    T tmp1;
	    	    
	    // Need to multiply:
	    //
	    //  (2/(1-m^2)) (2D(m))  or  (2/(1-m^2)) (2D^{dag}(m))
	    //
	    //  as appropriate. D(m) = m + (1 - m) D(0)
	    //  
	    // so  2D(m) = 2m + (1-m) 2D(0)
	    //           = [1+m-(1-m)] + (1-m) 2D(0)
	    //
	    // and 2/(1-m^2) = 2/((1+m)(1-m)
	    //
	    // finally: 
	    //
	    //  (2/(1-m^2)) (2D(m))  
	    //  = 2/(1+m) * [ (1+m)/(1-m) - 1 + 2D(0) ]
	    //  = 2/(1+m) * [ {(1+m)/(1-m) - 1} + M_scaled ]
	    //
	    //  
	    (*M_scaled)(tmp1, psi[j], isign_set[iset]);
	    
	    ftmp = (Real(1) + masses[i])/(Real(1) - masses[i]) - Real(1);
	
	    tmp1 += ftmp*psi[j];

	    tmp1 *= Real(2)/(Real(1) + masses[i]);
	    
	    
	    // Subtract off contact term 
	    psi[j] = tmp1 - chi;

	    // overall noramalisation 
	    ftmp = Real(1) / ( Real(1) - masses[i] );
	    psi[j] *= ftmp;
	    
	  }
	}
	break;
      default:
	QDP_error_exit("This value of n_soln is not known about %d\n", n_soln);
	break;
      }  // End switch over n_soln;
    }    
#if 0
    break; // End of CG case

    case REL_CG_INVERTER: 
    {
      // This is M_scaled = 2 D(0). 
      // the fact that D has zero masss is enforced earlier 
      //
      Handle< LinearOperator<T> > 
	// This is a really ugly construction but should work
	M_scaled(new approx_lopscl<T, Real>(dynamic_cast<LinearOperator<T>* >(linOp(state)),Real(2) ) );

      Chirality ischiral;
      LinearOperator<T>* MdagMPtr;
      multi1d<Real> shifted_masses(n_mass);

      for(int i = 0; i < n_mass; i++) { 
	shifted_masses[i] = ( Real(1) + masses[i] )/( Real(1) - masses[i] );
	shifted_masses[i] += ( Real(1) / shifted_masses[i] );
	shifted_masses[i] -= Real(2);
      }

      
      // This is icky. I have to new and delete, 
      // and I can't drop in a handle as I
      //    a) Can't declare variables in a case statement
      //       in case it jumps.
      // 
      // Hence I have to deal directly with the f***ing pointer
      // using new and delete.
      //
      //
      // MdagMPtr = 4 D^{dag}(0) D(0)
      //
      // The Zero is enforced elsewhere
      ischiral = isChiralVector(chi);
      if( ischiral == CH_NONE || ( isChiral() == false ) ) {

	// I have one scaled by 2. The MdagM fixes up the scale
	// factor of 4. I don't need to recompute coeffs etc here.
	//      
	MdagMPtr = new approx_lmdagm<T>(M_scaled);
      }
      else { 
	// Special lovddag version
	// This is yet another ugly cast.
	MdagMPtr = new approx_lopscl<T, Real> ( dynamic_cast<LinearOperator<T>* >(lMdagM( state, ischiral )), Real(4) );
      }
      
      // Do the solve
      //
      // This really solves with Kerel:
      //
      // (1/(1-m^2)) 4 D(m)^{dag} D(m)
      //  
      //  =  4 D(0)^{dag} D(0) + shift
      //
      // with shift = (1+m)/(1-m) + (1-m)/(1+m) - 2
      // 
      //
      // result of  solution is:
      //
      // psi = [ (1 - m^2 ) / 4 ] [ D^{dag}(m) D(m) ]^{-1} chi
      //
      //     = [ ( 1 - m^2 ) / 2 ] (2 D^{dag}(m))^{-1} D(m)^{-1} chi
      //
      //     = [ ( 1 - m^2 ) / 2 ] ( 2 D(m))^{-1} D^{-dag}(m) chi

      
      MInvRelCG(*MdagMPtr, chi, psi, shifted_masses, invParam.RsdCG, invParam.MaxCG, n_count);

      // Delete MdagMPtr. Do it right away before I forget.
      delete MdagMPtr;

      if ( n_count == invParam.MaxCG ) {
	QDP_error_exit("no convergence in the inverter", n_count);    
      }
      
      // Now make the solution(s)
      switch(n_soln) { 
      case 4:


	// Compensate for the  4/(1-m^2) factor

	for(int i = 0; i < n_mass; i++) { 
	  psi[i] *= Real(4) / ( Real(1) - masses[i]*masses[i] );
	}
	break;
	
      case 3:
	// Copy ofver the solutions of D^{-1} for solutions to D_dag^{-1}
	for( int i = 0 ; i < n_mass; i++) {
	  psi[n_mass + i] = psi[i];
	}
	// Fall through to the case below:

      case 1:
      case 2:
	for(int iset=0; iset < nsets; iset++) {
	  for(int i=0; i < n_mass; i++) {
	    int j = i+iset*n_mass;
	    T tmp1;
	    	    
	    // Need to multiply:
	    //
	    //  (2/(1-m^2)) (2D(m))  or  (2/(1-m^2)) (2D^{dag}(m))
	    //
	    //  as appropriate. D(m) = m + (1 - m) D(0)
	    //  
	    // so  2D(m) = 2m + (1-m) 2D(0)
	    //           = [1+m-(1-m)] + (1-m) 2D(0)
	    //
	    // and 2/(1-m^2) = 2/((1+m)(1-m)
	    //
	    // finally: 
	    //
	    //  (2/(1-m^2)) (2D(m))  
	    //  = 2/(1+m) * [ (1+m)/(1-m) - 1 + 2D(0) ]
	    //  = 2/(1+m) * [ {(1+m)/(1-m) - 1} + M_scaled ]
	    //
	    //  
	    (*M_scaled)(tmp1, psi[j], isign_set[iset]);
	    
	    ftmp = (Real(1) + masses[i])/(Real(1) - masses[i]) - Real(1);
	
	    tmp1 += ftmp*psi[j];

	    tmp1 *= Real(2)/(Real(1) + masses[i]);
	    
	    
	    // Subtract off contact term 
	    psi[j] = tmp1 - chi;

	    // overall noramalisation 
	    ftmp = Real(1) / ( Real(1) - masses[i] );
	    psi[j] *= ftmp;
	    
	  }
	}
	break;
      default:
	QDP_error_exit("This value of n_soln is not known about %d\n", n_soln);
	break;
      }  // End switch over n_soln;
    }    
    break; // End of CG case

    case SUMR_INVERTER:
    {

      /* Only do psi = D^(-1) chi */
      if (n_soln != 1) {
	QDP_error_exit("SUMMR inverter does not solve for D_dag");
      }

      Handle< LinearOperator<T> > U = lgamma5epsH(state);

      multi1d<Complex> shifted_masses(n_mass);

      // This is the code from SZIN. I checked the shifts . Should work ok
      for(int i=0; i < n_mass; ++i) {
	shifted_masses[i] = (Real(1) + masses[i]) / ( Real(1) - masses[i] );
      }
      
      // For the case of the overlap rho, is always 1
      multi1d<Real> rho(n_mass);
      rho = Real(1);

      multi1d<Real> scaledRsdCG(n_mass);
      for(int i=0; i < n_mass; i++) {
	scaledRsdCG[i] = RsdCG[i]*(Real(1) - masses[i])/Real(2);
      }

      // Do the solve
      MInvSUMR(*U, chi, psi, shifted_masses, rho, scaledRsdCG, invParam.MaxCG,n_count);

#if 1 
      for(int s=0; s < n_mass; s++) { 

	// Check back solutions
	T r;
	(*U)(r, psi[s], PLUS);
	r *= rho[s];
	r += shifted_masses[s]*psi[s];


	r -= chi;
	Double r_norm = sqrt(norm2(r))/sqrt(norm2(chi));
	QDPIO::cout << "Check: shift="<<s<<" || r ||/||b|| = " << r_norm << " RsdCG = " << RsdCG[s] << endl;
	
      }
#endif
      
      /* psi <- (D^(-1) - 1) chi */
      for(int i=0; i < n_mass; ++i) {
	
	// Go back to (1/2)( 1 + mu + (1 - mu) normalisation
	ftmp = Real(2) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
	// Remove conact term
	psi[i] -= chi;
	
	// overall noramalisation 
	ftmp = Real(1) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
      }
    }
    break; // End of MR inverter case
    case REL_SUMR_INVERTER:
    {

      /* Only do psi = D^(-1) chi */
      if (n_soln != 1) {
	QDP_error_exit("SUMMR inverter does not solve for D_dag");
      }

      Handle< LinearOperator<T> > U(lgamma5epsH(state));

      multi1d<Complex> shifted_masses(n_mass);

      // This is the code from SZIN. I checked the shifts . Should work ok
      for(int i=0; i < n_mass; ++i) {
	shifted_masses[i] = (Real(1) + masses[i]) / ( Real(1) - masses[i] );
      }
      
      // For the case of the overlap rho, is always 1
      multi1d<Real> rho(n_mass);
      rho = Real(1);

      // We are solving with a scaled operator 
      //  shifted_masses[i] + gamma_5 eps
      //
      // which is the original scaled by 2/(1 - m)
      // 
      // So I should rescale each desired RsdCG by (1-m)/2
      // where m is the smallest mass.
      multi1d<Real> scaledRsdCG(n_mass);
      for(int i=0; i < n_mass; i++) {
	scaledRsdCG[i] = RsdCG[i]*(1-masses[i])/Real(2);
      }

      // Do the solve
      MInvRelSUMR(*U, chi, psi, shifted_masses, rho, scaledRsdCG, invParam.MaxCG,n_count);
 
#if 1 
      for(int s=0; s < n_mass; s++) { 

	// Check back solutions
	T r;
	(*U)(r, psi[s], PLUS);
	r *= rho[s];
	r += shifted_masses[s]*psi[s];

	r -= chi;
	Double r_norm = sqrt(norm2(r))/sqrt(norm2(chi));
	QDPIO::cout << "Check: shift="<<s<<" || r ||/||b|| = " << r_norm << " RsdCG = " << RsdCG[s] << endl;
	
      }
#endif


      /* psi <- (D^(-1) - 1) chi */
      for(int i=0; i < n_mass; ++i) {
	
	// Go back to (1/2)( 1 + mu + (1 - mu) normalisation
	ftmp = Real(2) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
	// Remove conact term
	psi[i] -= chi;
	
	// overall noramalisation 
	ftmp = Real(1) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
      }
    }
    break; // End of MR inverter case
#if 0
    case MR_INVERTER:
    {

      /* psi = D^(-1) chi */
      if (n_soln != 1) {
	QDP_error_exit("MR inverter does not solve for D_dag");
      }

      // This is the code from SZIN. I checked the shifts . Should work ok
      multi1d<Real> shifted_masses(n_mass);
      for(int i=0; i < n_mass; ++i) {
	shifted_masses[i] = (Real(1) + masses[i]) / ( Real(1) - masses[i] ) - 1;
      }
      
      // Do the solve
      MInvMR (*M_scaled, chi, psi, shifted_masses, invParam.RsdCG, ncg_had);
      if ( n_count == invParam.MaxCG ) {
	QDP_error_exit("no convergence in the inverter", n_count);    
      }
      
      /* psi <- (D^(-1) - 1) chi */
      for(i=0; i < n_mass; ++i) {
	
	// Go back to (1/2)( 1 + mu + (1 - mu) normalisation
	ftmp = Real(2) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
	// Remove conact term
	psi[i] -= chi;
	
	// overall noramalisation 
	ftmp = Real(1) / ( Real(1) - masses[i] );
	psi[i] *= ftmp;
	
      }
    }
    break; // End of MR inverter case
#endif


    default:
      QDP_error_exit("Unknown inverter type %d\n", invParam.invType);
    }
#endif


    END_CODE();
  }
  
} // End Namespace Chroma

  
