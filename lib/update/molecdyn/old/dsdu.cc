// $Id: dsdu.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $

//! Computes the derivative of the Wilson (and Symanzik improved)
/*!
/* Computes the derivative of the Wilson (and Symanzik improved)
 * action with respect to the link field
 *
 * \param u    gauge field ( Read )
 *
 *           |  dS
 * ds_u -- U | ----  ( Write )
 *           |  dU 
 */
void dsdu(multi1d<LatticeColorMatrix>& ds_u,
	  multi1d<LatticeColorMatrix>& u)
{
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_4;
  LatticeColorMatrix tmp_5;
  LatticeColorMatrix tmp_6;
  LatticeColorMatrix   tmp_3;
  LatticeColorMatrix   tmp_tot;

  int mu;
  int nu;
  int rho;
  int eta;
  int cb;
  int mu0;
  int nu0;
  int rho0;
  int dir;
  int igluetmp;

  int k;
  int j;
  multi1d<int> fdir(Nd);  

  START_CODE("subroutine");;

  igluetmp = fabs(GlueImp);
  
  if (igluetmp > 2)
    QDP_error_exit("Illegal value of GlueImp", GlueImp);
                    
  ds_u = 0;
  
  for(mu = 0; mu < Nd; ++mu)
  {
    tmp_tot = 0;
    
    for(nu=mu+1;nu<Nd;nu++)
    {
      for(cb=0;cb<2;cb++)
      {
	tmp_0 = shift(u[mu][1-cb], cb, FORWARD, nu) * shift(adj[u[nu][1-cb]], cb, FORWARD, mu);
	tmp_1 = tmp_0 * adj(u[mu][cb]);
	tmp_2 = u[nu][cb] * tmp_1;
	ds_u[nu][cb] += tmp_2;
	ds_u[mu][cb] += adj(tmp_2);
	ds_u[mu][1-cb] += shift(tmp_1, (1-cb), BACKWARD, nu) * shift(u[nu][cb], (1-cb), BACKWARD, nu);
	tmp_1 = adj(u[nu][cb]) * u[mu][cb];
	ds_u[nu][1-cb] += shift(adj[tmp_0], (1-cb), BACKWARD, mu) * shift(tmp_1, (1-cb), BACKWARD, mu);
      }
    }      
    
    /*+ */
    /* Symanzik corrections */
    /*- */
    if (igluetmp >= 1)
    {
      /*  2*1 rectangles */

      /* calculate double_mu links */
      tmp_3[1] = u[mu][1] * shift(u[mu][0], 1, FORWARD, mu);
      tmp_3[0] = u[mu][0] * shift(u[mu][1], 0, FORWARD, mu);
    
      for(cb=0;cb<2;cb++)
      {
	tmp_6 = tmp_3[1-cb] * GlueCoeffRT;
	j = 0;
	for(nu = 0; nu < Nd; nu++)
	{
	  if(nu == mu)
	    continue;
	  /* forward plaquette */
	  tmp_1 = u[nu][cb] * shift(tmp_6, cb, FORWARD, nu);
	  tmp_5 = shift(adj[tmp_1], (1-cb), BACKWARD, mu) * shift(tmp_3[cb], (1-cb), BACKWARD, mu);
	  ds_u[nu][cb] += u[nu][cb] * shift(tmp_5, cb, BACKWARD, mu);
	  tmp_5(rb[(1-cb)]) = shift(u[nu][cb], FORWARD, mu);

	  /* at this point we add the nu contribution directly to ds_u */
	  /* we could make tmp_tot carry a direction index in order to avoid that */
	  if(j++ == 0)
	  {
	    tmp_4 = tmp_1 * shift(adj[tmp_5], cb, FORWARD, mu);
	    ds_u[nu][cb] += tmp_4 * adj(tmp_3[cb]);
	  }
	  else
	  {
	    tmp_2 = tmp_1 * shift(adj[tmp_5], cb, FORWARD, mu);
	    tmp_4 += tmp_2; /* sum to the staple */
	    ds_u[nu][cb] += tmp_2 * adj(tmp_3[cb]);
	  }
	  
	  /* backward plaquette */
	  tmp_1 = adj(u[nu][1-cb]) * tmp_6;
	  tmp_2(rb[cb]) = shift(u[nu][1-cb], FORWARD, mu);
	  tmp_5(rb[(1-cb)]) = shift(tmp_2, FORWARD, mu);
	  tmp_4 += shift(tmp_1, cb, BACKWARD, nu) * shift(tmp_5, cb, BACKWARD, nu);
	}
	tmp_tot[cb] += shift(u[mu][1-cb], cb, FORWARD, mu) * adj(tmp_4);
	tmp_tot[1-cb] += shift(adj[tmp_4], (1-cb), BACKWARD, mu) * shift(u[mu][cb], (1-cb), BACKWARD, mu);
      }
    }
    
    if(igluetmp == 2)
    {
      /* generalized parallelogram */
      for(nu=mu+1;nu<Nd;nu++)
      {
	j = 0;
	for(k=0;k<Nd;k++)
	  if((k!=mu)&&(k!=nu))
	    fdir[j++] = k;

	tmp_3 = 0;
	
	for(cb=0;cb<2;cb++)
	{
	  /* |_  part fitting piece is ~| */
	  tmp_0 = shift(u[mu][cb], (1-cb), FORWARD, nu) * shift(adj[u[nu][cb]], (1-cb), FORWARD, mu);
	  tmp_5 = tmp_0 * GlueCoeffPG;

	  for(k=0;k<Nd-2;k++)
	  {
	    eta = fdir[k];
	    tmp_2 = shift(u[eta][1-cb], cb, FORWARD, nu) * shift(tmp_5, cb, FORWARD, eta);

	    if(k==0)
	    {
	      tmp_4 = tmp_2 * shift(adj[u[eta][1-cb]], cb, FORWARD, mu);
	    }
	    else
	    {
              tmp_4 += tmp_2 * shift(adj[u[eta][1-cb]], cb, FORWARD, mu);
	    }
	  }
	  tmp_3[cb] += tmp_4 * adj(u[mu][cb]);
	  tmp_tot[cb] += adj(tmp_4) * adj(u[nu][cb]);

	  /* _| piece */
	  tmp_0 = u[nu][cb] * shift(u[mu][1-cb], cb, FORWARD, nu);
	  tmp_5 = tmp_0 * GlueCoeffPG;
	  for(k=0;k<Nd-2;k++)
	  {
	    eta = fdir[k];
	    tmp_2 = shift(adj[u[eta][cb]], (1-cb), BACKWARD, mu) * shift(tmp_5, (1-cb), BACKWARD, mu);
	    tmp_1(rb[(1-cb)]) = shift(u[eta][cb], FORWARD, nu);
	    if(k==0)
	    {
	      tmp_4 = shift(tmp_2, cb, BACKWARD, eta) * shift(tmp_1, cb, BACKWARD, eta);
	    }
	    else
	    {
	      tmp_4 += shift(tmp_2, cb, BACKWARD, eta) * shift(tmp_1, cb, BACKWARD, eta);
	    }
	  }
	  tmp_3[cb] += adj(tmp_4) * shift(u[mu][1-cb], cb, BACKWARD, mu);
	  tmp_tot[1-cb] += shift(u[nu][cb], (1-cb), FORWARD, mu) * shift(adj[tmp_4], (1-cb), FORWARD, mu);

	  /* ~| part */
	  tmp_0 = adj(u[nu][1-cb]) * u[mu][1-cb];
	  tmp_5 = shift(tmp_0, cb, BACKWARD, nu) * GlueCoeffPG;
	  for(k=0;k<Nd-2;k++)
	  {
	    eta = fdir[k];
	    tmp_2 = u[eta][1-cb] * shift(tmp_5, (1-cb), FORWARD, eta);
	    if(k==0)
	    {
	      tmp_4 = shift(tmp_2, cb, BACKWARD, mu) * shift(adj[u[eta][1-cb]], cb, BACKWARD, nu);
	    }
	    else
	    {
	      tmp_4 += shift(tmp_2, cb, BACKWARD, mu) * shift(adj[u[eta][1-cb]], cb, BACKWARD, nu);
	    }
	  }

	  tmp_1(rb[cb]) = shift(u[nu][1-cb], BACKWARD, nu);
	  tmp_tot[1-cb] += shift(adj[tmp_1], (1-cb), FORWARD, mu) * shift(adj[tmp_4], (1-cb), FORWARD, mu);
	  tmp_1(rb[cb]) = shift(u[mu][1-cb], BACKWARD, mu);
	  tmp_3[1-cb] += shift(adj[tmp_1], (1-cb), FORWARD, nu) * shift(tmp_4, (1-cb), FORWARD, nu);
	  
	  /* Finally the |~ part */
	  tmp_0 = u[mu][cb] * shift(u[nu][1-cb], cb, FORWARD, mu);
	  tmp_5 = tmp_0 * GlueCoeffPG;
	  for(k=0;k<Nd-2;k++)
	  {
	    eta = fdir[k];
	    tmp_2 = shift(adj[u[eta][cb]], (1-cb), BACKWARD, nu) * shift(tmp_5, (1-cb), BACKWARD, nu);
	    tmp_1(rb[(1-cb)]) = shift(u[eta][cb], FORWARD, mu);
	    if(k==0)
	    {
	      tmp_4 = shift(tmp_2, cb, BACKWARD, eta) * shift(tmp_1, cb, BACKWARD, eta);
	    }
	    else
	    {
	      tmp_4 += shift(tmp_2, cb, BACKWARD, eta) * shift(tmp_1, cb, BACKWARD, eta);
	    }
	  }
	  tmp_tot[cb] += adj(tmp_4) * shift(u[nu][1-cb], cb, BACKWARD, nu);
	  tmp_3[1-cb] += shift(u[mu][cb], (1-cb), FORWARD, nu) * shift(adj[tmp_4], (1-cb), FORWARD, nu);

	} /* end loop over cb */	  

	/* mult tmp_3 by u(x,nu) matrix) */
	ds_u[nu][0] += u[nu][0] * tmp_3[0];
	ds_u[nu][1] += u[nu][1] * tmp_3[1];
      } /* end loop over nu */
    } /* end( igluetmp == 2)  */
    
    
    /* ds_u =  u(x,mu) * tmp_tot */
    for(cb = 0; cb < 2; ++cb)
      ds_u[mu][cb] += u[mu][cb] * tmp_tot[cb];
  }
  
  
                    
  END_CODE("subroutine");;
}
