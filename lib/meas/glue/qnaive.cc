// $Id: qnaive.cc,v 2.0 2005-09-25 21:04:34 edwards Exp $

#error "NOT FULLY CONVERTED"

/* Computes the naive topological charge, with the "clover" definition */
/* of the naive topological charge density. It also computes the */
/* ratio of the action to the continuum instaton action. */

/* u        -- gauge field ( Read ) */
/* qtop     -- topological charge ( Write ) */
/* S_ratio  -- ratio of actio to continuum instanton action ( Write ) */

include(types.mh)

SUBROUTINE(qnaive, u, qtop, S_ratio)

multi1d<LatticeColorMatrix> u(Nd);
Double qtop;
Double S_ratio;
{
  include(COMMON_DECLARATIONS)
  
  /* Local Variables */
  LatticeColorMatrix u_clov1;
  LatticeColorMatrix u_clov2;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeReal plaq_tmp;
  LatticeReal qtop_tmp;
  
  Double plaq;
  Double tmp;
  int mu1;
  int nu1;
  int mu2;
  int nu2;
  int cb;
  int n_pl;
  
  START_CODE();
  
  if( Nd != 4 )
    QDP_error_exit("Nd for the topological charge has to be 4 but: ", Nd);

          
  plaq = 0;
  qtop = 0;
  
  /* Loop over checkerboards and triplet of perpendicular planes */
  mu1 = 0;
  for(cb=0; cb<2; ++cb)
    for(nu1=1; nu1<Nd; ++nu1)
    {
      /* First "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1(rb[cb]) = shift(u[nu1][1-cb], FORWARD, mu1);

      /* tmp_2(x) = u(x+nu1,mu1) */
      tmp_2(rb[cb]) = shift(u[mu1][1-cb], FORWARD, nu1);

      /* tmp_1(x) = u(x+mu1,nu1) * u_dag(x+nu1,mu1) */
      tmp_1 = shift(u[nu1][1-cb], cb, FORWARD, mu1) * shift(adj[u[mu1][1-cb]], cb, FORWARD, nu1);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu1)
		  = u(x+mu1,nu1) * u_dag(x+nu1,mu1) * u_dag(x,nu1) */
      tmp_2 = tmp_1 * adj(u[nu1][cb]);

      /* u_clov1(x) = u(x,mu1) * tmp_2
		    = u(x,mu1) * u(x+mu1,nu1) *
		      u_dag(x+nu1,mu1) * u_dag(x,nu1) */
      u_clov1 = u[mu1][cb] * tmp_2;

            plaq_tmp = real(trace(u_clov1));
      plaq += sum(plaq_tmp);
      
      /* First "plus-minus" plaquette */
      /* tmp_1(x) = u(x,mu1) * u(x+mu1,nu1) */
      tmp_1 = u[mu1][1-cb] * shift(u[nu1][cb], 1-cb, FORWARD, mu1);

      /* tmp_2(x) = tmp_1_dag(x-nu1) * u(x-nu1,nu1)
		  = u_dag(x+mu1-nu1,nu1) * u_dag(x-nu1,mu1) * u(x-nu1,nu1) */
      tmp_2 = shift(adj[tmp_1], cb, BACKWARD, nu1) * shift(u[nu1][1-cb], cb, BACKWARD, nu1);

      /* tmp_1(x) = u(x,mu1) * tmp_2 = u(x,mu1) * u_dag(x-nu1+mu1,nu1) *
				       u_dag(x-nu1,mu1) * u(x-nu1,nu1) */
      tmp_1 = u[mu1][cb] * tmp_2;

      u_clov1 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      /* First "minus-minus" plaquette */
      /* tmp_1(x) = u(x,mu1) * u(x+mu1,nu1) */
      tmp_1 = u[mu1][cb] * shift(u[nu1][1-cb], cb, FORWARD, mu1);

      /* tmp_2(x) = u_dag(x-nu1,nu1) * tmp_1(x-nu1)
		  = u_dag(x-nu1,nu1) * u(x-nu1,mu1) * u(x+mu1-nu1,nu1) */
      tmp_2 = shift(adj[u[nu1][cb]], 1-cb, BACKWARD, nu1) * shift(tmp_1, 1-cb, BACKWARD, nu1);

      /* tmp_1(x) = u_dag(x-mu1,mu1) * tmp_2(x-mu1)
		  = u_dag(x-mu1,mu1) * u_dag(x-nu1-mu1,nu1) *
		    u(x-nu1-mu1,mu1) * u(x-nu1,nu1) */
      tmp_1 = shift(adj[u[mu1][1-cb]], cb, BACKWARD, mu1) * shift(tmp_2, cb, BACKWARD, mu1);

      u_clov1 += tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      /* First "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu1) * u(x,nu1) */
      tmp_1 = adj(u[mu1][1-cb]) * u[nu1][1-cb];

      /* tmp_2(x) = tmp_1(x) * u(x+nu1,mu1)
		  = u_dag(x,mu1) * u(x,nu1) * u(x+nu1,mu1) */
      tmp_2 = tmp_1 * shift(u[mu1][cb], 1-cb, FORWARD, nu1);

      /* tmp_1(x) = tmp_2(x-mu1) * u_dag(x,nu1)
		  = u_dag(x-mu1,mu1) * u(x-mu1,nu1) *
		    u(x-mu1+nu1,mu1) * u_dag(x,nu1) */
      tmp_1 = shift(tmp_2, cb, BACKWARD, mu1) * adj(u[nu1][cb]);

      u_clov1 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      mu2 = INTEGER_MOD_FUNCTION(nu1,3) + 1;
      nu2 = INTEGER_MOD_FUNCTION(mu2,3) + 1;

      /* Second "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1(rb[cb]) = shift(u[nu2][1-cb], FORWARD, mu2);

      /* tmp_2(x) = u(x+nu2,mu2) */
      tmp_2(rb[cb]) = shift(u[mu2][1-cb], FORWARD, nu2);

      /* tmp_1(x) = u(x+mu2,nu2) * u_dag(x+nu2,mu2) */
      tmp_1 = shift(u[nu2][1-cb], cb, FORWARD, mu2) * shift(adj[u[mu2][1-cb]], cb, FORWARD, nu2);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu2)
		  = u(x+mu2,nu2) * u_dag(x+nu2,mu2) * u_dag(x,nu2) */
      tmp_2 = tmp_1 * adj(u[nu2][cb]);

      /* u_clov2(x) = u(x,mu2) * tmp_2
		    = u(x,mu2) * u(x+mu2,nu2) *
		      u_dag(x+nu2,mu2) * u_dag(x,nu2) */
      u_clov2 = u[mu2][cb] * tmp_2;

            plaq_tmp = real(trace(u_clov2));
      plaq += sum(plaq_tmp);
      
      /* Second "plus-minus" plaquette */
      /* tmp_1(x) = u(x,mu2) * u(x+mu2,nu2) */
      tmp_1 = u[mu2][1-cb] * shift(u[nu2][cb], 1-cb, FORWARD, mu2);

      /* tmp_2(x) = tmp_1_dag(x-nu2) * u(x-nu2,nu2)
		  = u_dag(x+mu2-nu2,nu2) * u_dag(x-nu2,mu2) * u(x-nu2,nu2) */
      tmp_2 = shift(adj[tmp_1], cb, BACKWARD, nu2) * shift(u[nu2][1-cb], cb, BACKWARD, nu2);

      /* tmp_1(x) = u(x,mu2) * tmp_2 = u(x,mu2) * u_dag(x-nu2+mu2,nu2) *
				       u_dag(x-nu2,mu2) * u(x-nu2,nu2) */
      tmp_1 = u[mu2][cb] * tmp_2;

      u_clov2 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      /* Second "minus-minus" plaquette */
      /* tmp_1(x) = u(x,mu2) * u(x+mu2,nu2) */
      tmp_1 = u[mu2][cb] * shift(u[nu2][1-cb], cb, FORWARD, mu2);

      /* tmp_2(x) = u_dag(x-nu2,nu2) * tmp_1(x-nu2)
		  = u_dag(x-nu2,nu2) * u(x-nu2,mu2) * u(x+mu2-nu2,nu2) */
      tmp_2 = shift(adj[u[nu2][cb]], 1-cb, BACKWARD, nu2) * shift(tmp_1, 1-cb, BACKWARD, nu2);

      /* tmp_1(x) = u_dag(x-mu2,mu2) * tmp_2(x-mu2)
		  = u_dag(x-mu2,mu2) * u_dag(x-nu2-mu2,nu2) *
		    u(x-nu2-mu2,mu2) * u(x-nu2,nu2) */
      tmp_1 = shift(adj[u[mu2][1-cb]], cb, BACKWARD, mu2) * shift(tmp_2, cb, BACKWARD, mu2);

      u_clov2 += tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      /* Second "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu2) * u(x,nu2) */
      tmp_1 = adj(u[mu2][1-cb]) * u[nu2][1-cb];

      /* tmp_2(x) = tmp_1(x) * u(x+nu2,mu2)
		  = u_dag(x,mu2) * u(x,nu2) * u(x+nu2,mu2) */
      tmp_2 = tmp_1 * shift(u[mu2][cb], 1-cb, FORWARD, nu2);

      /* tmp_1(x) = tmp_2(x-mu2) * u_dag(x,nu2)
		  = u_dag(x-mu2,mu2) * u(x-mu2,nu2) *
		    u(x-mu2+nu2,mu2) * u_dag(x,nu2) */
      tmp_1 = shift(tmp_2, cb, BACKWARD, mu2) * adj(u[nu2][cb]);

      u_clov2 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      plaq += sum(plaq_tmp);
      
      /* Now comes the contribution to the topological charge */
      tmp_1 = 1;
      tmp_2 = adj(u_clov1) * tmp_1;
      u_clov1 -= tmp_2;
      tmp_2 = adj(u_clov2) * tmp_1;
      u_clov2 -= tmp_2;
      tmp_2 = u_clov1 * u_clov2;

            qtop_tmp = real(trace(tmp_2));
      qtop -= sum(qtop_tmp);
      
    }
  
          
  if ( SchrFun == 10 )
  {
    mu1 = nrow[0] - 2;
    nu1 = nrow[1] - 2;
    mu2 = nrow[2] - 2;
    nu2 = nrow[3] - 2;
    n_pl = 6 * mu1 * nu1 * mu2 * nu2;
    n_pl -= 3 * (nu1*mu2*nu2 + mu1*mu2*nu2 + mu1*nu1*nu2 + mu1*nu1*mu2);
    n_pl += mu2*nu2 + nu1*nu2 + nu1*mu2 + mu1*nu2 + mu1*mu2 + mu1*nu1;
  }
  else
    n_pl = Nd * (Nd-1) * vol_cb;

  plaq = plaq / TO_DOUBLE(4*n_pl*Nc);
  
  /* Topological charge */
  FILL(tmp,PI);
  qtop = qtop / ( TO_DOUBLE(256) * tmp*tmp );
  
  /* Ratio of action to continuum instanton action */
  S_ratio = (WORD_VALUE(WORD_plaq,ONE) - plaq);
  S_ratio = S_ratio * TO_DOUBLE(n_pl*Nc) / ( TO_DOUBLE(4)* tmp*tmp );
  
  END_CODE();
}
