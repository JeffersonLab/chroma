// $Id: overlap_fermact_base_w.cc,v 1.13 2004-05-21 12:03:13 bjoo Exp $
/*! \file
 *  \brief Base class for unpreconditioned overlap-like fermion actions
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "actions/ferm/invert/invcg1.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/invert/inv_rel_cg1.h"
#include "actions/ferm/invert/inv_rel_cg2.h"

#include "actions/ferm/invert/invsumr.h"
#include "actions/ferm/invert/inv_rel_sumr.h"
#include "actions/ferm/invert/minvsumr.h"
#include "actions/ferm/invert/minvcg.h"
#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/linop/lopscl.h"
#include "meas/eig/ischiral_w.h"


using namespace QDP;

//! Propagator for unpreconditioned overlap-like fermion actions
/*!
 * \param psi      quark propagator ( Modify )
 * \param state_   gauge field ( Read )
 * \param chi      source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void 
OverlapFermActBase::qprop(LatticeFermion& psi, 
			  Handle<const ConnectState> state, 
			  const LatticeFermion& chi, 
			  enum InvType invType,
			  const Real& RsdCG, 
			  int MaxCG, int& n_count) const
{
  START_CODE("OverlapFermActBase::qprop");

  Real mass = quark_mass();

  switch( invType ) {
  case CG_INVERTER:
    {

      LatticeFermion tmp;
      Handle< const LinearOperator<LatticeFermion> > M(linOp(state));
	
      // Check whether the source is chiral.
      Chirality ichiral = isChiralVector(chi);
      if( ichiral == CH_NONE || ( isChiral() == false )) { 
	
	

	(*M)(tmp, chi, MINUS);
      
	// Source is not chiral. In this case we should use,
	// InvCG2 with M
	InvCG2(*M, tmp, psi, RsdCG, MaxCG, n_count);
      }
      else {
	
	// Source is chiral. In this case we should use InvCG1
	// with the special MdagM
	Handle< const LinearOperator<LatticeFermion> > MM(lMdagM(state, ichiral));

	InvCG1(*MM, chi, tmp, RsdCG, MaxCG, n_count);
	(*M)(psi, tmp, MINUS);	
      }
  
      LatticeFermion Mpsi;
      (*M)(Mpsi, psi, PLUS);
      Mpsi = chi - Mpsi;
      Mpsi /= sqrt(norm2(chi));
      QDPIO::cout << "OvQprop || chi - D psi ||/||chi|| = " << sqrt(norm2(Mpsi))
		  << "  n_count = " << n_count << " iters" << endl;

    }
    break;
  case REL_CG_INVERTER:
    {

      LatticeFermion tmp;
      Handle< const ApproxLinearOperator<LatticeFermion> > M( dynamic_cast<const ApproxLinearOperator<LatticeFermion>* > ( linOp(state) ) );


      // Check whether the source is chiral.
      Chirality ichiral = isChiralVector(chi);
      if( ichiral == CH_NONE || ( isChiral() == false )) { 
	
	(*M)(tmp, chi, MINUS);
      
	// Source is not chiral. In this case we should use,
	// InvCG2 with M
	InvRelCG2(*M, tmp, psi, RsdCG, MaxCG, n_count);
      }
      else {
	
	// Source is chiral. In this case we should use InvCG1
	// with the special MdagM
	Handle< const ApproxLinearOperator<LatticeFermion> > MM( dynamic_cast<const ApproxLinearOperator<LatticeFermion>* >( lMdagM(state, ichiral) ) );

	InvRelCG1(*MM, chi, tmp, RsdCG, MaxCG, n_count);
	(*M)(psi, tmp, MINUS);	
      }

      LatticeFermion Mpsi;
      (*M)(Mpsi, psi, PLUS);
      Mpsi = chi - Mpsi;
      Mpsi /= sqrt(norm2(chi));
      QDPIO::cout << "OvQprop || chi - D psi ||/||chi|| = " << sqrt(norm2(Mpsi))
		  << "  n_count = " << n_count << " iters" << endl;

    }
    break;

#if 0
  case MR_INVERTER:
    // psi = D^(-1)* chi
    InvMR (*M, chi, psi, MRover, RsdCG, MaxCG, n_count);
    break;

  case BICG_INVERTER:
    // psi = D^(-1) chi
    InvBiCG (*M, chi, psi, RsdCG, MaxCG, n_count);
    break;
#endif

  case SUMR_INVERTER:
    {
      // Solve by SUMR solver -- for shifted unitary matrices
      //
      // Solve zeta I + rho gamma_5 eps(H)
      // where gamma_5 eps(H) is unitary
      //
      // zeta = (1 + mu)/(1-mu)
      // rho  = 1
      Real rho = Real(1);
      Real mu = quark_mass();
      Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
      {
	Handle<const LinearOperator<LatticeFermion> > U(lgamma5epsH(state));
	
	// Now solve:
	InvSUMR(*U, chi, psi, zeta, rho, RsdCG, MaxCG, n_count);
	
	// Restore to normal scaling
	Real fact = Real(2)/(Real(1) - mu);
	psi *= fact;
	
	
      }

      // Check back:
      // Get a proper operator and compute chi- Dpsi
      Handle<const LinearOperator<LatticeFermion> > D(linOp(state));
      LatticeFermion Dpsi;
      (*D)(Dpsi, psi, PLUS);
      Dpsi = chi - Dpsi;
      Dpsi /= sqrt(norm2(chi));
      QDPIO::cout << "OvQprop || chi - D psi || = " << sqrt(norm2(Dpsi))
		  << "  n_count = " << n_count << " iters" << endl;
    }
    break;

  case REL_SUMR_INVERTER:
    {
      // Solve by Relaxed SUMR solver -- for shifted unitary matrices
      //
      // Solve zeta I + rho gamma_5 eps(H)
      // where gamma_5 eps(H) is unitary
      //
      // zeta = (1 + mu)/(1-mu)
      // rho  = 1
      Real rho = Real(1);
      Real mu = quark_mass();
      Complex zeta = cmplx(( Real(1) + mu ) / (Real(1) - mu),0);
      {
	Handle<const ApproxLinearOperator<LatticeFermion> > U( dynamic_cast<const ApproxLinearOperator<LatticeFermion>* >(lgamma5epsH(state)) );

	Real fact = Real(2)/(Real(1) - mu);
	
	// Now solve:
	InvRelSUMR(*U, chi, psi, zeta, rho, RsdCG, MaxCG, n_count);
	
	// Restore to normal scaling
	psi *= fact;
	
	
      }

      // Check back:
      // Get a proper operator and compute chi- Dpsi
      Handle<const LinearOperator<LatticeFermion> > D(linOp(state));
      LatticeFermion Dpsi;
      (*D)(Dpsi, psi, PLUS);
      Dpsi = chi - Dpsi;
      Dpsi /= sqrt(norm2(chi));
      QDPIO::cout << "OvQprop || chi - D psi || = " << sqrt(norm2(Dpsi))
		  << "  n_count = " << n_count << " iters" << endl;
    }
    break;

  default:
    QDP_error_exit("Zolotarev4DFermActBj::qprop Solver Type not implemented\n");
    break;
  };

  if ( n_count == MaxCG ) { 
    QDP_error_exit("Zolotarev4DFermAct::qprop: No convergence in solver: n_count = %d\n", n_count);
  }

  // Normalize and remove contact term 
  Real ftmp = Real(1) / ( Real(1) - mass );
  
  psi -= chi;
  psi *= ftmp;

  END_CODE("OverlapFermActBase::qprop");
}
  

/* This routine is Wilson type Overlap fermions */

/* Compute multiple quark mass propagators for an unpreconditioned overla fermion */
/* using the single mass source in "chi" - so, the source can */
/* be of any desired form. The results will appear in "psi", which is */
/* ignored on input (and set to zero initially) */

/* This routine can return the solution to the following systems,
   when using CG */
/*   D * psi           = chi   (nsoln = 1) */
/*   D^dag * psi       = chi   (nsoln = 2) */
/*  and all of the above for nsoln = 3 */
/*   (D * D^dag) * psi = chi   (for nsoln = 4) */

/* u        -- gauge field ( Read ) */
/* chi      -- source ( Modify ) */
/* Mass     -- quark masses in lattice units ( Read ) */
/* numMass  -- number of quark masses ( Read ) */
/* RsdCG    -- CG or MR residual ( Read ) */
/* psi      -- quark propagators ( Write ) */
/* ncg_had  -- number of CG iterations ( Modify ) */
void 
OverlapFermActBase::multiQprop(multi1d<LatticeFermion>& psi,
			       const multi1d<Real>& masses,
			       Handle<const ConnectState> state, 
			       const LatticeFermion& chi, 
			       enum InvType invType,
			       const multi1d<Real>& RsdCG, 
			       int nsoln,
			       int MaxCG, 
			       int& n_count) const

{

  if ( toBool( quark_mass() != 0 ) ) {
    QDP_error_exit("Multi Mass Qprop only works if action has mass = 0 (strictly)");
  }


  int nsets ;
  int msets ;
  multi1d<enum PlusMinus> isign_set(2);


  switch(nsoln) {
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
    QDP_error_exit("solution type incorrect", nsoln);
  }
  

  const int n_mass = masses.size();



  
  
  Real ftmp;
  switch( invType ) {

  case CG_INVERTER: 
    {
      // This is M_scaled = 2 D(0). 
      // the fact that D has zero masss is enforced earlier 
      //
      const Handle< const LinearOperator<LatticeFermion> > 
	M_scaled( new lopscl<LatticeFermion, Real>( linOp(state), Real(2) ) );

      Chirality ischiral;
      LinearOperator<LatticeFermion>* MdagMPtr;
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
      
	MdagMPtr = new lmdagm<LatticeFermion>(M_scaled);
      }
      else { 
	// Special lovddag version
	MdagMPtr = new lopscl<LatticeFermion, Real> ( lMdagM( state, ischiral ), Real(4) );
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
      MInvCG(*MdagMPtr, chi, psi, shifted_masses, RsdCG, MaxCG, n_count);

      // Delete MdagMPtr. Do it right away before I forget.
      delete MdagMPtr;

      if ( n_count == MaxCG ) {
	QDP_error_exit("no convergence in the inverter", n_count);    
      }
      
      // Now make the solution(s)
      switch(nsoln) { 
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
	    LatticeFermion tmp1;
	    	    
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
	QDP_error_exit("This value of nsoln is not known about %d\n", nsoln);
	break;
      }  // End switch over nsoln;
    }    
    break; // End of CG case

  case SUMR_INVERTER:
    {

      /* Only do psi = D^(-1) chi */
      if (nsoln != 1) {
	QDP_error_exit("SUMMR inverter does not solve for D_dag");
      }

      Handle<const LinearOperator<LatticeFermion> > U = lgamma5epsH(state);

      multi1d<Complex> shifted_masses(n_mass);

      // This is the code from SZIN. I checked the shifts . Should work ok
      for(int i=0; i < n_mass; ++i) {
	shifted_masses[i] = (Real(1) + masses[i]) / ( Real(1) - masses[i] );
      }
      
      // For the case of the overlap rho, is always 1
      multi1d<Real> rho(n_mass);
      rho = Real(1);

      // Do the solve
      MInvSUMR(*U, chi, psi, shifted_masses, rho, RsdCG, MaxCG,n_count);

      if ( n_count == MaxCG ) {
	QDP_error_exit("no convergence in the inverter", n_count);    
      }
      
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
      if (nsoln != 1) {
	QDP_error_exit("MR inverter does not solve for D_dag");
      }

      // This is the code from SZIN. I checked the shifts . Should work ok
      multi1d<Real> shifted_masses(n_mass);
      for(int i=0; i < n_mass; ++i) {
	shifted_masses[i] = (Real(1) + masses[i]) / ( Real(1) - masses[i] ) - 1;
      }
      
      // Do the solve
      MInvMR (*M_scaled, chi, psi, shifted_masses, RsdCG, ncg_had);
      if ( n_count == MaxCG ) {
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
    QDP_error_exit("Unknown inverter type %d\n", invType);
  }


  END_CODE("qprop");
}
  
  
