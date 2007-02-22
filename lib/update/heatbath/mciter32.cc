// $Id: mciter32.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#error "Not tested (or even compiled). However, reasonably well converted."

#error "THIS VERSION SHOULD BE USED IN PLACE OF MCITER.CC - IT IS MORE GENERAL. REMOVE THIS FILE"

#warning "THERE ARE INSTANCES OF shift(...,mu)*shift(...,mu)  THAT SHOULD BE OPTIMIZED"


#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "update/heatbath/su3over.h"
#include "update/heatbath/su3hb.h"

namespace Chroma {


//! One heatbath interation of updating the gauge field configuration
/*!
 * \ingroup heatbath
 *
 * Make one interation of updating the gauge field configuration:
 *      for Wilson or Symmanzik improved pure gauge action:
 *      this consists of n_over overrelaxation sweeps followed
 *      by one heatbath sweep with nheat trials.
 * In the case of SU(3), for each link we loop over the 3 SU(2) subgroups.
 *
 * Because of the option for the Symanzik improved action we need 2^(d+1)
 *      sublattices, deviding first into 2^d hypercubes and then
 *       checkerboarding those. 
 *
 *
 * Warning: this works only for Nc = 2 and 3 !
 *
 * \param u        gauge field ( Modify )
 * \param n_over   number of overrelaxation sweeps ( Read )
 * \param nheat    number of heatbath trials ( Read )
 * \param ntrials  total number of individual heatbath trials ( Modify )
 * \param nfails   total number of individual heatbath failures ( Modify ) 
 */

void mciter32(multi1d<LatticeColorMatrix>& u, 
	      int n_over, int nheat,
	      int& ntrials, int& nfails,
	      const Set& ss,
	      const multi3d<int>& neighsubl)
{
  START_CODE();

  // Expect neighsubl to be declared  multi3d<int> neighsubl(Nd,2,ss.numSubsets());

  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix u_staple;
  LatticeColorMatrix u_staple2;


#if 0
  igluetmp = fabs(GlueImp);
  if (igluetmp > 2)
    QDP_error_exit("Illegal value of GlueImp", GlueImp);
#endif

  
  for(int iter = 0; iter <= n_over; ++iter)
  {
    for(int subl = 0; subl < ss.numSubsets(); ++subl)
      for(int mu = 0; mu < Nd; ++mu)
      {
	u_staple[ss[subl]] = 0;

	int fsublm = neighsubl[mu][1][subl];
	int bsublm = neighsubl[mu][0][subl];
	
	for(int nu = 0; nu < Nd; ++nu)
	{
	  if(nu == mu)  continue;

	  int fsubln = neighsubl[nu][1][subl];
	  int bsubln = neighsubl[nu][0][subl];

	  /* Forward staple */
	  /* tmp_1(x) = u(x+mu,nu)*u_dag(x+nu,mu) */
	  tmp_1[ss[subl]] = shift(u[nu], FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

	  /* u_staple(x) +=  tmp_1 * u_dag(x,nu)
	     += u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  u_staple[ss[subl]] += tmp_1 * adj(u[nu]);

	  /* Backward staple */
	  int bfsubl = neighsubl[mu][1][bsubln];

	  /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
	  tmp_1[ss[bsubln]] = u[mu] * shift(u[nu], FORWARD, mu);

	  /* u_staple(x) += tmp_1_dag(x-nu) * u(x-nu,nu)
	     += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu) */
	  u_staple[ss[subl]] += shift(adj(tmp_1), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);
	}

	/* Symmanzik improvement terms */
	if (igluetmp >= 1)
	{
	  u_staple2 = 0;
	  
	  for(int nu = 0; nu < Nd; ++nu)
	  {
	    if(nu == mu)  continue;

	    int fsubln = neighsubl[nu][1][subl];
	    int bsubln = neighsubl[nu][0][subl];

	    /* forward 1x2 rectangle in mu,nu directions, nu!=mu */
	    int fsublmn = neighsubl[mu][1][fsubln];
	    int fsubl2 = neighsubl[nu][1][fsubln];

	    /* tmp_1(x) = u(x+mu,nu)*u_dag(x+nu,mu) */
	    tmp_1[ss[fsubln]] = shift(u[nu], FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

	    /* tmp_2(x) = tmp_1(x+nu)*u_dag(x,nu)
	       = u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	    tmp_2[ss[fsubln]] = tmp_1 * adj(u[nu]);

	    /* tmp_1(x) = tmp_2(x+nu)*u_dag(x,nu)
	       = u(x+nu+mu,nu)*u_dag(x+2nu,mu)*
	       u_dag(x+nu,nu)*u_dag(x,nu) */
	    tmp_1[ss[subl]] = shift(tmp_2, FORWARD, nu) * adj(u[nu]);

	    /* u_staple2 += u(x+mu,nu)*tmp_1(x) */
	    u_staple2[ss[subl]] += shift(u[nu], FORWARD, mu) * tmp_1;

	    /* backward 1x2 rectangle in mu,nu directions, nu!=mu */
	    int bsubl2 = neighsubl[nu][0][bsubln];
	    int bfsubl = neighsubl[mu][1][bsubl2];

	    /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
	    tmp_1[ss[bsubl2]] = u[mu] * shift(u[nu], FORWARD, mu);

	    /* tmp_2(x) = tmp_1_dag(x-nu) * u(x-nu,nu)
	       = u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu) */
	    tmp_2[ss[bsubln]] = shift(adj(tmp_1), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	    /* tmp_1(x) = u_dag(x+mu,nu)*tmp_2(x)
	       = u_dag(x+mu,nu)*u_dag(x+mu-nu,nu)*
	       u_dag(x-nu,mu)*u(x-nu,nu) */
	    int bfsubl = neighsubl[mu][1][bsubln];
	    tmp_1[ss[bsubln]] = shift(adj(u[nu]), FORWARD, mu) * tmp_2;

	    /* u_staple2(x) += tmp_1(x-nu) * u(x-nu,nu) */
	    u_staple2[ss[subl]] += shift(tmp_1, BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	    /* forward 2x1 rectangle in mu,nu directions open at bottom left */
	    int fsubl2 = neighsubl[mu][1][fsublm];
	    int fsublmn = neighsubl[nu][1][fsublm];

	    /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
	    tmp_1[ss[fsublm]] = u[mu] * shift(u[nu], FORWARD, mu);

	    /* tmp_2(x) = tmp_1(x)*u_dag(x+nu,mu)
	       = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu) */
	    tmp_2[ss[fsublm]] = tmp_1 * shift(adj(u[mu]), FORWARD, nu);

	    /* tmp_1(x) = tmp_2(x+mu)*u_dag(x+nu,mu)
	       = u(x+mu,mu)*u(x+2mu,nu)*
	       u_dag(x+mu+nu,mu)*u_dag(x+nu,mu) */
	    tmp_1[ss[subl]] = shift(tmp_2, FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

	    /* u_staple2(x) +=  tmp_1 * u_dag(x,nu) */
	    u_staple2[ss[subl]] += tmp_1 * adj(u[nu]);

	    /* backward 2x1 rectangle in mu,nu directions open at top left */
	    int bfsubl = neighsubl[mu][1][bsubln];
	    int fsubl2 = neighsubl[mu][1][bfsubl];

	    /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
	    tmp_1[ss[bfsubl]] = u[mu] * shift(u[nu], FORWARD, mu);

	    /* tmp_2(x) = u(x,mu)*tmp_1(x+mu)
	       = u(x,mu)*u(x+mu,mu)*u(x+2mu,nu) */
	    tmp_2[ss[bsubln]] = u[mu] * shift(tmp_1, FORWARD, mu);

	    /* tmp_1(x) = tmp_2_dag(x-nu)*u(x-nu,nu)
	       = u_dag(x-nu+2mu,nu)*u_dag(x-nu+mu,mu)*
	       u_dag(x-nu,mu)*u(x-nu,nu) */
	    tmp_1[ss[subl]] = shift(adj(tmp_2), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	    /* u_staple2(x) +=  u(x+mu,mu) * tmp_1 */
	    u_staple2[ss[subl]] += shift(u[mu], FORWARD, mu) * tmp_1;

	    /* forward 2x1 rectangle in mu,nu directions open at bottom right */
	    /* tmp_1(x) = u(x,nu)*u(x+nu,mu) */
	    int bfsubl = neighsubl[nu][1][bsublm];
	    tmp_1[ss[bsublm]] = u[nu] * shift(u[mu], FORWARD, nu);

	    /* tmp_2(x) = tmp_1_dag(x-mu)*u(x-mu,mu)
	       = u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu) */
	    tmp_2[ss[subl]] = shift(adj(tmp_1), BACKWARD, mu) * shift(u[mu], BACKWARD, mu);

	    /* tmp_1(x) = u_dag(x+nu,mu)*tmp_2(x)
	       = u_dag(x+nu,mu)*u_dag(x-mu+nu,mu)*
	       u_dag(x-mu,nu)*u(x-mu,mu) */
	    tmp_1[ss[subl]] = shift(adj(u[mu]), FORWARD, nu) * tmp_2;

	    /* u_staple2(x) +=  u(x+mu,nu) * tmp_1 */
	    u_staple2[ss[subl]] += shift(u[nu], FORWARD, mu) * tmp_1;

	    /* backward 2x1 rectangle in mu,nu directions open at top right */
	    int bfsubl = neighsubl[mu][1][bsubln];
	    int bsublmn = neighsubl[mu][0][bsubln];

	    /* tmp_1(x) = u_dag(x,mu)*u(x,nu) */
	    tmp_1[ss[bsublmn]] = adj(u[mu]) * u[nu];

	    /* tmp_2(x) = u_dag(x,mu)*tmp_1(x-mu)
	       = u_dag(x,mu)*u_dag(x-mu,mu)*u(x-mu,nu) */
	    tmp_2[ss[bsubln]] = adj(u[mu]) * shift(tmp_1, BACKWARD, mu);

	    /* tmp_1(x) = u_dag(x+mu,nu)*tmp_2(x)
	       = u_dag(x+mu,nu)*u_dag(x,mu)*
	       u_dag(x-mu,mu)*u(x-mu,nu) */
	    tmp_1[ss[bsubln]] = shift(adj(u[nu]), FORWARD, mu) * tmp_2;

	    /* u_staple2(x) += tmp_1(x-nu) * u(x-mu,mu) */
	    u_staple2[ss[subl]] += shift(tmp_1, BACKWARD, nu) * shift(u[mu], BACKWARD, mu);
	  }

	  u_staple[ss[subl]] += u_staple2 * GlueCoeffRT;
	}

	/* Symmanzik 6-link parallelograms */
	if (igluetmp == 2)
	{
	  u_staple2 = 0;
	  
	  /* In each cube, labeled by eta<nu<rho there are 4 ways
	     to get a 6-link parallelogram: (eta,nu,rho,-eta,-nu,-rho),
	     (nu,eta,rho,-nu,-eta,-rho), (eta,rho,nu,-eta,-rho,-nu)
	     and (eta,-nu,rho,-eta,nu,-rho). The first three can be
	     obtained by keeping only (eta,nu,rho,-eta,-nu,-rho) 
	     while restricting only to eta<rho, nu!=eta, nu!=rho. */
	  for(int eta = 0; eta < Nd-1; ++eta)
	    for(int nu = 0; nu < Nd; ++nu)
	    {
	      if(nu == eta)  continue;

	      int fsuble = neighsubl[eta][1][subl];
	      int bsuble = neighsubl[eta][0][subl];
	      int fsubln = neighsubl[nu][1][subl];
	      int bsubln = neighsubl[nu][0][subl];

	      for(int rho = eta+1; rho < Nd; ++rho)
	      {
	        if(rho == nu)  continue;
		
	        int fsublr = neighsubl[rho][1][subl];
	        int bsublr = neighsubl[rho][0][subl];

	        if(mu == eta)
	        {
	          /* "update" link is first link */
	          /* tmp_1(x) = u(x,nu)*u(x+nu,mu) */
	          int fsublmn = neighsubl[nu][1][fsublr];
	          tmp_1[ss[fsublr]] = u[nu] * shift(u[mu], FORWARD, nu);

	          /* tmp_2(x) = u(x,rho)*tmp_1(x+rho)
		     = u(x,rho)*u(x+rho,nu)*u(x+rho+nu,mu) */
	          tmp_2[ss[subl]] = u[rho] * shift(tmp_1, FORWARD, rho);

	          /* tmp_1(x) = u(x,nu)*u(x+nu,rho) */
	          int fsublmn = neighsubl[nu][1][fsublm];
	          tmp_1[ss[fsublm]] = u[nu] * shift(u[rho], FORWARD, nu);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2_dag(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * adj(tmp_2);

	          /* "update" link is fourth link */
	          int bsubl2 = neighsubl[rho][0][bsubln];
	          int bfsubl = neighsubl[mu][1][bsubl2];
	          int fsubl2 = neighsubl[nu][1][bfsubl];

	          /* tmp_1(x) = u(x,nu)*u(x+nu,rho) */
	          tmp_1[ss[bfsubl]] = u[nu] * shift(u[rho], FORWARD, nu);

	          /* tmp_2(x) = u(x,mu)*tmp_1(x+mu)
		     = u(x,mu)*u(x+mu,nu)*u(x+mu+nu,rho) */
	          tmp_2[ss[bsubl2]] = u[mu] * shift(tmp_1, FORWARD, mu);

	          /* tmp_1(x) = tmp_2_dag(x-rho)*u(x-rho,rho)
		     = u_dag(x-rho+mu+nu,rho)*u_dag(x-rho+mu,nu)*
		     u_dag(x-rho,mu)*u(x-rho,rho) */
	          tmp_1[ss[bsubln]] = shift(adj(tmp_2), BACKWARD, rho) * shift(u[rho], BACKWARD, rho);

	          /* u_staple2(x) += tmp_1(x-nu) * u(x-nu,nu) */
	          u_staple2[ss[subl]] += shift(tmp_1, BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	        }
	        else if(mu == nu)
	        {
	          /* "update" link is second link */
	          int bfsubl = neighsubl[rho][1][bsuble];
	          int fsubl2 = neighsubl[mu][1][bfsubl];

	          /* tmp_1(x) = u(x,mu)*u(x+mu,eta) */
	          tmp_1[ss[bfsubl]] = u[mu] * shift(u[eta], FORWARD, mu);

	          /* tmp_2(x) = u(x,rho)*tmp_1(x+rho)
		     = u(x,rho)*u(x+rho,mu)*u(x+rho+mu,eta) */
	          tmp_2[ss[bsuble]] = u[rho] * shift(tmp_1, FORWARD, rho);

	          /* tmp_1(x) = tmp_2_dag(x)*u(x,eta)
		     = u_dag(x+rho+mu,eta)*u_dag(x+rho,mu)*
		     u_dag(x,rho)*u(x,eta) */
	          tmp_1[ss[bsuble]] = adj(tmp_2) * u[eta];

	          /* u_staple2(x) += u(x+mu,rho) * tmp_1(x-eta) */
	          u_staple2[ss[subl]] += shift(u[rho], FORWARD, mu) * shift(tmp_1, BACKWARD, eta);

	          /* "update" link is fifth link */
	          int bfsubl = neighsubl[eta][1][bsublr];
	          int fsubl2 = neighsubl[mu][1][bfsubl];

	          /* tmp_1(x) = u(x,mu)*u(x+mu,rho) */
	          tmp_1[ss[bfsubl]] = u[mu] * shift(u[rho], FORWARD, mu);

	          /* tmp_2(x) = u(x,eta)*tmp_1(x+eta)
		     = u(x,eta)*u(x+eta,mu)*u(x+eta+mu,rho) */
	          tmp_2[ss[bsublr]] = u[eta] * shift(tmp_1, FORWARD, eta);

	          /* tmp_1(x) = tmp_2_dag(x)*u(x,rho)
		     = u_dag(x+eta+mu,eta)*u_dag(x+eta,mu)*
		     u_dag(x,eta)*u(x,rho) */
	          tmp_1[ss[bsublr]] = adj(tmp_2) * u[rho];

	          /* u_staple2(x) += u(x+mu,eta) * tmp_1(x-rho) */
	          u_staple2[ss[subl]] += shift(u[eta], FORWARD, mu) * shift(tmp_1, BACKWARD, rho);

	        }
	        else if(mu == rho)
	        {
	          /* "update" link is third link */
	          int bsubl2 = neighsubl[eta][0][bsubln];
	          int bfsubl = neighsubl[mu][1][bsubl2];
	          int fsubl2 = neighsubl[nu][1][bfsubl];

	          /* tmp_1(x) = u(x,nu)*u(x+nu,eta) */
	          tmp_1[ss[bfsubl]] = u[nu] * shift(u[eta], FORWARD, nu);

	          /* tmp_2(x) = u(x,mu)*tmp_1(x+mu)
		     = u(x,mu)*u(x+mu,nu)*u(x+mu+nu,eta) */
	          tmp_2[ss[bsubl2]] = u[mu] * shift(tmp_1, FORWARD, mu);

	          /* tmp_1(x) = tmp_2_dag(x-eta)*u(x-eta,eta)
		     = u_dag(x-eta+mu+nu,eta)*u_dag(x-eta+mu,nu)*
		     u_dag(x-eta,mu)*u(x-eta,eta) */
	          tmp_1[ss[bsubln]] = shift(adj(tmp_2), BACKWARD, eta) * shift(u[eta], BACKWARD, eta);

	          /* u_staple2(x) += tmp_1(x-nu) * u(x-nu,nu) */
	          u_staple2[ss[subl]] += shift(tmp_1, BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	          /* "update" link is sixth link */
	          /* tmp_1(x) = u(x,nu)*u(x+nu,mu) */
	          int fsublmn = neighsubl[nu][1][fsuble];
	          tmp_1[ss[fsuble]] = u[nu] * shift(u[mu], FORWARD, nu);

	          /* tmp_2(x) = u(x,eta)*tmp_1(x+eta)
		     = u(x,eta)*u(x+eta,nu)*u(x+eta+nu,mu) */
	          tmp_2[ss[subl]] = u[eta] * shift(tmp_1, FORWARD, eta);

	          /* tmp_1(x) = u(x,nu)*u(x+nu,eta) */
	          int fsublmn = neighsubl[nu][1][fsublm];
	          tmp_1[ss[fsublm]] = u[nu] * shift(u[eta], FORWARD, nu);
		  
	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2_dag(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * adj(tmp_2);
	        }
	      }
	    }

	  /* Finally the (eta,-nu,rho,-eta,nu,-rho) parallelogram */
	  for(int eta = 0; eta < Nd-2; ++eta)
	    for(int nu = eta+1; nu < Nd-1; ++nu)
	    {
	      int fsuble = neighsubl[eta][1][subl];
	      int bsuble = neighsubl[eta][0][subl];
	      int fsubln = neighsubl[nu][1][subl];
	      int bsubln = neighsubl[nu][0][subl];

	      for(int rho = nu+1; rho < Nd; ++rho)
	      {
	        int fsublr = neighsubl[rho][1][subl];
	        int bsublr = neighsubl[rho][0][subl];

	        if(mu == eta)
	        {
	          /* "update" link is first link */
	          /* tmp_1(x) = u_dag(x-nu,mu)*u(x-nu,nu) */
	          int bfsubl = neighsubl[nu][0][fsublr];
	          tmp_1[ss[fsublr]] = shift(adj(u[mu]), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	          /* tmp_2(x) = tmp_1(x+rho)*u_dag(x,rho)
		     = u_dag(x+rho-nu,mu)*u(x+rho-nu,nu)*u_dag(x,rho) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, rho) * adj(u[rho]);

	          /* tmp_1(x) = u_dag(x-nu,nu)*u(x-nu,rho) */
	          int bfsubl = neighsubl[nu][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[nu]), BACKWARD, nu) * shift(u[rho], BACKWARD, nu);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;

	          /* "update" link is forth link */
	          /* tmp_1(x) = u_dag(x-rho,mu)*u(x-rho,rho) */
	          int bfsubl = neighsubl[rho][0][fsubln];
	          tmp_1[ss[fsubln]] = shift(adj(u[mu]), BACKWARD, rho) * shift(u[rho], BACKWARD, rho);

	          /* tmp_2(x) = tmp_1(x+nu)*u_dag(x,nu)
		     = u_dag(x+nu-rho,mu)*u(x+nu-rho,rho)*u_dag(x,nu) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, nu) * adj(u[nu]);

	          /* tmp_1(x) = u_dag(x-rho,rho)*u(x-rho,nu) */
	          int bfsubl = neighsubl[rho][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[rho]), BACKWARD, rho) * shift(u[nu], BACKWARD, rho);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;

	        }
	        else if(mu == nu)
	        {
	          /* "update" link is second link */
	          /* tmp_1(x) = u_dag(x-eta,mu)*u(x-eta,eta) */
	          int bfsubl = neighsubl[eta][0][fsublr];
	          tmp_1[ss[fsublr]] = shift(adj(u[mu]), BACKWARD, eta) * shift(u[eta], BACKWARD, eta);

	          /* tmp_2(x) = tmp_1(x+rho)*u_dag(x,rho)
		     = u_dag(x+rho-eta,mu)*u(x+rho-eta,eta)*
		     u_dag(x,rho) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, rho) * adj(u[rho]);

	          /* tmp_1(x) = u_dag(x-eta,eta)*u(x-eta,rho) */
	          int bfsubl = neighsubl[eta][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[eta]), BACKWARD, eta) * shift(u[rho], BACKWARD, eta);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;

	          /* "update" link is fifth link */
	          /* tmp_1(x) = u_dag(x-rho,mu)*u(x-rho,rho) */
	          int bfsubl = neighsubl[rho][0][fsuble];
	          tmp_1[ss[fsuble]] = shift(adj(u[mu]), BACKWARD, rho) * shift(u[rho], BACKWARD, rho);

	          /* tmp_2(x) = tmp_1(x+eta)*u_dag(x,eta)
		     = u_dag(x+eta-rho,mu)*u(x+eta-rho,rho)*
		     u_dag(x,eta) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, eta) * adj(u[eta]);

	          /* tmp_1(x) = u_dag(x-rho,rho)*u(x-rho,eta) */
	          int bfsubl = neighsubl[rho][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[rho]), BACKWARD, rho) * shift(u[eta], BACKWARD, rho);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;

	        }
	        else if(mu == rho)
	        {
	          /* "update" link is third link */
	          /* tmp_1(x) = u_dag(x-eta,mu)*u(x-eta,eta) */
	          int bfsubl = neighsubl[eta][0][fsubln];
	          tmp_1[ss[fsubln]] = shift(adj(u[mu]), BACKWARD, eta) * shift(u[eta], BACKWARD, eta);

	          /* tmp_2(x) = tmp_1(x+nu)*u_dag(x,nu)
		     = u_dag(x+nu-eta,mu)*u(x+nu-eta,eta)*
		     u_dag(x,nu) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, nu) * adj(u[nu]);

	          /* tmp_1(x) = u_dag(x-eta,eta)*u(x-eta,nu) */
	          bfsubl = neighsubl[eta][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[eta]), BACKWARD, eta) * shift(u[nu], BACKWARD, eta);

	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;

	          /* "update" link is sixth link */
	          /* tmp_1(x) = u_dag(x-nu,mu)*u(x-nu,nu) */
	          bfsubl = neighsubl[nu][0][fsuble];
	          tmp_1[ss[fsuble]] = shift(adj(u[mu]), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);

	          /* tmp_2(x) = tmp_1(x+eta)*u_dag(x,eta)
		     = u_dag(x+eta-nu,mu)*u(x+eta-nu,nu)*
		     u_dag(x,eta) */
	          tmp_2[ss[subl]] = shift(tmp_1, FORWARD, eta) * adj(u[eta]);

	          /* tmp_1(x) = u_dag(x-nu,nu)*u(x-nu,eta) */
	          bfsubl = neighsubl[nu][0][fsublm];
	          tmp_1[ss[fsublm]] = shift(adj(u[nu]), BACKWARD, nu) * shift(u[eta], BACKWARD, nu);
		  
	          /* u_staple2(x) += tmp_1(x+mu) * tmp_2(x) */
	          u_staple2[ss[subl]] += shift(tmp_1, FORWARD, mu) * tmp_2;
	        }
	      }
	    }

	  u_staple[ss[subl]] += u_staple2 * GlueCoeffPG;
	}

	
	if ( iter < n_over )
	{
	  /* Do an overrelaxation step */
          /*# Loop over SU(2) subgroup index */
          for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
            su3over(u[mu], u_staple, su2_index, ss[subl]);
	}
	else
	{
	  /* Do a heatbath step */
	  /*# Loop over SU(2) subgroup index */
          for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	  {
	    int ntry;
	    int nfail;

	    su3hb(u[mu], u_staple, su2_index, nheat, ntry, nfail, ss[subl]);
	    ntrials = ntrials + ntry;
	    nfails = nfails + nfail;
	  }

          
	  /* Reunitarize */
/*	  reunit (u[mu][subl], lbad, OPTION[REUNITARIZE_ERROR], numbad);  */
	  reunit (u[mu][subl]);

	}

#if 0
	/* If using Schroedinger functional, reset the boundaries */
	if ( SchrFun > 0 )
	{
	  copymask(u[mu][subl], lSFmask[mu][subl], SFBndFld[mu][subl], REPLACE);
	}
#endif
      }
  }
  
  END_CODE();
}

}  // end namespace Chroma
