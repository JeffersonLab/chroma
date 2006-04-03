// $Id: cool.cc,v 3.0 2006-04-03 04:58:57 edwards Exp $

#error "NOT FULLY CONVERTED"


/* Make one interation of "cooling" of a gauge field configuration: */
/*      minimize the gauge action locally for each link. */
/* In the case of SU(3), for each link we loop over the 3 SU(2) subgroups. */

/* Warning: this works only for Nc = 2 and 3 ! */

/* u        -- gauge field ( Modify ) */
include(types.mh)

SUBROUTINE(cool, u)

multi1d<LatticeColorMatrix> u(Nd);

{ /* Local Variables */
  include(COMMON_DECLARATIONS)
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_tot;
  LatticeBoolean bad;                       /* not used here */
  
  int mu;
  int nu;
  int cb;
  int su2_index;
  int numbad;                            /* not used here */
  
  START_CODE();
  
    
  for(cb = 0;cb  <= ( 1); ++cb )
    for(mu = 0;mu  < ( Nd); ++mu )
    {
      tmp_tot = 0;

      
      for(nu = 0;nu  < ( Nd); ++nu )
      {
	if(nu == mu)  continue;

	/* Forward staple */
	/* tmp_1(x) = u(x+mu,nu) * u_dag(x+nu,mu) */
	tmp_1 = shift(u[nu][1-cb], cb, FORWARD, mu) * shift(adj[u[mu][1-cb]], cb, FORWARD, nu);

	/* tmp_tot(x) += tmp_1(x) * u_dag(x,nu)
		      += u(x+mu,nu) * u_dag(x+nu,mu)*u_dag(x,nu) */
	tmp_tot += tmp_1 * adj(u[nu][cb]);

	/* Backward staple */
	/* tmp_1(x) = u(x,mu) * u(x+mu,nu) */
	tmp_1 = u[mu][1-cb] * shift(u[nu][cb], 1-cb, FORWARD, mu);

	/* tmp_tot(x) += tmp_1_dag(x-nu) * u(x-nu,nu)
		      += u_dag(x+mu-nu,nu) * u_dag(x-nu,mu) * u(x-nu,nu) */
	tmp_tot += shift(adj[tmp_1], cb, BACKWARD, nu) * shift(u[nu][1-cb], cb, BACKWARD, nu);
      }

      
      /* Loop over SU(2) subgroup su2_index to maximize tr(u*tmp_tot) */
      if (SchrFun > 0)
      {
	/* Make it easy, since these are overwritten anyway! */
	FILLMASK(u(cb,mu), lSFmask(cb,mu), ONE);
	FILLMASK(tmp_tot, lSFmask(cb,mu), ONE);
      }
      for(su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	su3proj (u[mu][cb], tmp_tot, su2_index);

            /* Reunitarize */
      /*#reunit (u[mu][cb], bad, OPTION[REUNITARIZE_ERROR], numbad); */
      reunit (u[mu][cb], bad, OPTION[REUNITARIZE], numbad);
      
      if (SchrFun > 0)
      {
	/* Now do the overwrite. */
	copymask(u[mu][cb], lSFmask[mu][cb], SFBndFld[mu][cb], REPLACE);
      }
    }
  
    
  END_CODE();
}
