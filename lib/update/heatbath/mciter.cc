// $Id: mciter.cc,v 1.2 2004-01-22 22:51:37 edwards Exp $
/*! \file
 *  \brief One heatbath interation of updating the gauge field configuration
 */

#error "Not tested (or even compiled). However, reasonably well converted."

#warning "THERE ARE INSTANCES OF shift(...,mu)*shift(...,mu)  THAT SHOULD BE OPTIMIZED"

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "update/heatbath/su3over.h"
#include "update/heatbath/su3hb.h"

using namespace QDP;

//! One heatbath interation of updating the gauge field configuration
/*!
 * \ingroup heatbath
 *
 * Make one interation of updating the gauge field configuration:
 *      this consists of n_over overrelaxation sweeps followed
 *      by one heatbath sweep with nheat trials.
 * In the case of SU(3), for each link we loop over the 3 SU(2) subgroups.

 * Warning: this works only for Nc = 2 and 3 !

 * \param u        gauge field ( Modify )
 * \param n_over   number of overrelaxation sweeps ( Read )
 * \param nheat    number of heatbath trials ( Read )
 * \param ntrials  total number of individual heatbath trials ( Modify )
 * \param nfails   total number of individual heatbath failures ( Modify ) 
 */

void mciter(multi1d<LatticeColorMatrix>& u, 
	    int n_over, int nheat,
	    int& ntrials, int& nfails)
{
  START_CODE("mciter");

  LatticeColorMatrix tmp_1;
  LatticeColorMatrix u_staple;

#if 0
  Real xi02;
  if( AnisoP == YES )
    xi02 = xi_0 * xi_0;  
#else
  bool AnisoP = false;
  Real xi02 = 1.0;
#endif

  for(int iter = 0; iter <= n_over; ++iter)
  {
    for(int cb = 0; cb < 2; ++cb)
      for(int mu = 0; mu < Nd; ++mu)
      {
	u_staple = 0;

	for(int nu = 0; nu < Nd; ++nu)
	{
	  if (nu == mu) continue;

	  /* Forward staple */
	  /* tmp_1(x) = u(x+mu,nu)*u_dag(x+nu,mu) */
	  tmp_1[rb[cb]] = shift(u[nu], FORWARD, mu) * shift(adj(u[mu]), FORWARD, nu);

          if( AnisoP )  	
	    if( mu == t_dir || nu == t_dir )
	      tmp_1[rb[cb]] *= xi02;

	  /* u_staple(x) +=  tmp_1 * u_dag(x,nu)
	     += u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  u_staple[rb[cb]] += tmp_1 * adj(u[nu]);
          

	  /* Backward staple */
	  /* tmp_1(x) = u(x,mu)*u(x+mu,nu) */
	  tmp_1[rb[1-cb]] = u[mu] * shift(u[nu], FORWARD, mu);

          if( AnisoP )  	
	    if( mu == t_dir || nu == t_dir )
	      tmp_1[rb[1-cb]] = xi02 * tmp_1;

	  /* u_staple(x) += tmp_1_dag(x-nu) * u(x-nu,nu)
	     += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu) */
          u_staple[rb[cb]] += shift(adj(tmp_1), BACKWARD, nu) * shift(u[nu], BACKWARD, nu);
	}  /* closes nu loop */

	if ( iter < n_over )
	{
	  /* Do an overrelaxation step */
          /*# Loop over SU(2) subgroup index */
          for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
            su3over(u[mu], u_staple, su2_index, rb[cb]);
	}
	else
	{
	  /* Do a heatbath step */
	  /*# Loop over SU(2) subgroup index */
          for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	  {
	    int ntry;
	    int nfail;

	    su3hb(u[mu], u_staple, su2_index, nheat, ntry, nfail, rb[cb]);
	    ntrials = ntrials + ntry;
	    nfails = nfails + nfail;
	  }

          
	  /* Reunitarize */
/*	  reunit (u[mu][cb], lbad, OPTION[REUNITARIZE_ERROR], numbad);  */
	  reunit(u[mu]);

	}

#if 0
	/* If using Schroedinger functional, reset the boundaries */
	if ( SchrFun > 0 )
	{
	  copymask(u[mu][cb], lSFmask[mu][cb], SFBndFld[mu][cb], REPLACE);
	}
#endif
      }        /* closes mu loop */
  }

  
  END_CODE("mciter");
}
