// $Id: qactden.cc,v 1.1 2005-02-17 02:52:11 edwards Exp $

#error "NOT FULLY CONVERTED"

/* QACTDEN - measure the lattice density of the lattice energy and the naive  */
/*           topological charge. */

/* u       -- gauge field (Read) */
/* lrqtop  -- topological charge density (Write) */
/* lract   -- action to continuum instanton action density (Write) */

include(types.mh)

SUBROUTINE(QActDen, u, lrqtop, lract)

multi1d<LatticeColorMatrix> u(Nd);
LatticeReal lract;
LatticeReal lrqtop;
{
  include(COMMON_DECLARATIONS)
  
  /* other local variables */
  LatticeColorMatrix u_clov1;
  LatticeColorMatrix u_clov2;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_3;
  LatticeReal qtop_tmp;
  LatticeReal plaq_tmp;
  
  Real rtmp;
  Double tmp;
  int t;
  int c;
  int cb;
  int cb_0;
  int mu;
  int nu;
  int mu1;
  int nu1;
  int mu2;
  int nu2;
  
  START_CODE();
  
            
  /* Lattice version of S_ratio */
  rtmp = TO_REAL(2*Nd*(Nd-1)*Nc);
  FILL(lract(0), rtmp);
  FILL(lract(1), rtmp);
  
  /* Lattice version of Qtop */
  lrqtop = 0;
  
  /* Loop over checkerboards and triplet of perpendicular planes */
  mu1 = 0;
  for(cb=0; cb < 2; ++cb)
    for(nu1=1; nu1 < Nd; ++nu1)
    {
      /* First "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1(rb[cb]) = shift(u[nu1][1-cb], FORWARD, mu1);

      /* tmp_2(x) = u(x+nu1,mu1) */
      tmp_2(rb[cb]) = shift(u[mu1][1-cb], FORWARD, nu1);

      /* tmp_3(x) = tmp_1 * tmp_2^dag = u(x+mu1,nu1)*u_dag(x+nu1,mu1) */
      tmp_3 = tmp_1 * adj(tmp_2);

      /* tmp_1(x) = tmp_3 * u_dag(x,nu1)= u(x+mu1,nu1)*u_dag(x+nu1,mu1)*u_dag(x,nu1) */
      tmp_1 = tmp_3 * adj(u[nu1][cb]);

      /* u_clov1(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1)* */
      /*                                 u_dag(x+nu1,mu1)*u_dag(x,nu1) */
      u_clov1 = u[mu1][cb] * tmp_1;
      
            plaq_tmp = real(trace(u_clov1));
      lract[cb] -= plaq_tmp;
      
      /* First "plus-minus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1(rb[1-cb]) = shift(u[nu1][cb], FORWARD, mu1);

      /* tmp_2(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1) */
      tmp_2 = u[mu1][1-cb] * tmp_1;

      /* tmp_1(x) = tmp_2_dag * u(x,nu1) = u_dag(x+mu1,nu1)*u_dag(x,mu1)*u(x,nu1) */
      tmp_1 = adj(tmp_2) * u[nu1][1-cb];

      /* tmp_2(x) = tmp_1(x-nu1) */
      tmp_2(rb[cb]) = shift(tmp_1, BACKWARD, nu1);

      /* tmp_1(x) = u(x,mu1) * tmp_2 = u(x,mu1)*u_dag(x-nu1+mu1,nu1)* */
      /*                               u_dag(x-nu1,mu1)*u(x-nu1,nu1) */
      tmp_1 = u[mu1][cb] * tmp_2;

      u_clov1 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      lract[cb] -= plaq_tmp;
      
      /* First "minus-minus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1(rb[cb]) = shift(u[nu1][1-cb], FORWARD, mu1);

      /* tmp_2(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1) */
      tmp_2 = u[mu1][cb] * tmp_1;

      /* tmp_1(x) = u_dag(x,nu1) * tmp_2 = u_dag(x,nu1)*u(x,mu1)*u(x+mu1,nu1) */
      tmp_1 = adj(u[nu1][cb]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-nu1) */
      tmp_2(rb[1-cb]) = shift(tmp_1, BACKWARD, nu1);

      /* tmp_1(x) = u_dag(x,mu1) * tmp_2 = u_dag(x,mu1)*u_dag(x-nu1,nu1)* */
      /*                                   u(x-nu1,mu1)*u(x-nu1+mu1,nu1) */
      tmp_1 = adj(u[mu1][1-cb]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-mu1) */
      tmp_2(rb[cb]) = shift(tmp_1, BACKWARD, mu1);

      u_clov1 += tmp_2;

            plaq_tmp = real(trace(tmp_2));
      lract[cb] -= plaq_tmp;
      
      /* First "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu1) * u(x,nu1) */
      tmp_1 = adj(u[mu1][1-cb]) * u[nu1][1-cb];

      /* tmp_2(x) = u(x+nu1,mu1) */
      tmp_2(rb[1-cb]) = shift(u[mu1][cb], FORWARD, nu1);

      /* tmp_3(x) = tmp_1 * tmp_2 = u_dag(x,mu1)*u(x,nu1)*u(x+nu1,mu1) */
      tmp_3 = tmp_1 * tmp_2;

      /* tmp_1(x) = tmp_3(x-mu1) */
      tmp_1(rb[cb]) = shift(tmp_3, BACKWARD, mu1);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu1) = u_dag(x-mu1,mu1)*u(x-mu1,nu1)* */
      /*                                   u(x-mu1+nu1,mu1)*u_dag(x,nu1) */
      tmp_2 = tmp_1 * adj(u[nu1][cb]);

      u_clov1 -= tmp_2;

            plaq_tmp = real(trace(tmp_2));
      lract[cb] -= plaq_tmp;
      
      mu2 = INTEGER_MOD_FUNCTION(nu1,3) + 1;
      nu2 = INTEGER_MOD_FUNCTION(mu2,3) + 1;

      /* Second "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1(rb[cb]) = shift(u[nu2][1-cb], FORWARD, mu2);

      /* tmp_2(x) = u(x+nu2,mu2) */
      tmp_2(rb[cb]) = shift(u[mu2][1-cb], FORWARD, nu2);

      /* tmp_3(x) = tmp_1 * tmp_2^dag = u(x+mu2,nu2)*u_dag(x+nu2,mu2) */
      tmp_3 = tmp_1 * adj(tmp_2);

      /* tmp_1(x) = tmp_3 * u_dag(x,nu2)= u(x+mu2,nu2)*u_dag(x+nu2,mu2)*u_dag(x,nu2) */
      tmp_1 = tmp_3 * adj(u[nu2][cb]);

      /* u_clov2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2)* */
      /*                                 u_dag(x+nu2,mu2)*u_dag(x,nu2) */
      u_clov2 = u[mu2][cb] * tmp_1;

            plaq_tmp = real(trace(u_clov2));
      lract[cb] -= plaq_tmp;
      
      /* Second "plus-minus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1(rb[1-cb]) = shift(u[nu2][cb], FORWARD, mu2);

      /* tmp_2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2) */
      tmp_2 = u[mu2][1-cb] * tmp_1;

      /* tmp_1(x) = tmp_2_dag * u(x,nu2) = u_dag(x+mu2,nu2)*u_dag(x,mu2)*u(x,nu2) */
      tmp_1 = adj(tmp_2) * u[nu2][1-cb];

      /* tmp_2(x) = tmp_1(x-nu2) */
      tmp_2(rb[cb]) = shift(tmp_1, BACKWARD, nu2);

      /* tmp_1(x) = u(x,mu2) * tmp_2 = u(x,mu2)*u_dag(x-nu2+mu2,nu2)* */
      /*                               u_dag(x-nu2,mu2)*u(x-nu2,nu2) */
      tmp_1 = u[mu2][cb] * tmp_2;

      u_clov2 -= tmp_1;

            plaq_tmp = real(trace(tmp_1));
      lract[cb] -= plaq_tmp;
      
      /* Second "minus-minus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1(rb[cb]) = shift(u[nu2][1-cb], FORWARD, mu2);

      /* tmp_2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2) */
      tmp_2 = u[mu2][cb] * tmp_1;

      /* tmp_1(x) = u_dag(x,nu2) * tmp_2 = u_dag(x,nu2)*u(x,mu2)*u(x+mu2,nu2) */
      tmp_1 = adj(u[nu2][cb]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-nu2) */
      tmp_2(rb[1-cb]) = shift(tmp_1, BACKWARD, nu2);

      /* tmp_1(x) = u_dag(x,mu2) * tmp_2 = u_dag(x,mu2)*u_dag(x-nu2,nu2)* */
      /*                                   u(x-nu2,mu2)*u(x-nu2+mu2,nu2) */
      tmp_1 = adj(u[mu2][1-cb]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-mu2) */
      tmp_2(rb[cb]) = shift(tmp_1, BACKWARD, mu2);

      u_clov2 += tmp_2;

            plaq_tmp = real(trace(tmp_2));
      lract[cb] -= plaq_tmp;
      
      /* Second "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu2) * u(x,nu2) */
      tmp_1 = adj(u[mu2][1-cb]) * u[nu2][1-cb];

      /* tmp_2(x) = u(x+nu2,mu2) */
      tmp_2(rb[1-cb]) = shift(u[mu2][cb], FORWARD, nu2);

      /* tmp_3(x) = tmp_1 * tmp_2 = u_dag(x,mu2)*u(x,nu2)*u(x+nu2,mu2) */
      tmp_3 = tmp_1 * tmp_2;

      /* tmp_1(x) = tmp_3(x-mu2) */
      tmp_1(rb[cb]) = shift(tmp_3, BACKWARD, mu2);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu2) = u_dag(x-mu2,mu2)*u(x-mu2,nu2)* */
      /*                                   u(x-mu2+nu2,mu2)*u_dag(x,nu2) */
      tmp_2 = tmp_1 * adj(u[nu2][cb]);

      u_clov2 -= tmp_2;

            plaq_tmp = real(trace(tmp_2));
      lract[cb] -= plaq_tmp;
      
      /* Now comes the contribution to the topological charge */
      tmp_1 = 1;
      tmp_2 = adj(u_clov1) * tmp_1;
      u_clov1 -= tmp_2;
      tmp_2 = adj(u_clov2) * tmp_1;
      u_clov2 -= tmp_2;
      tmp_2 = u_clov1 * u_clov2;

            qtop_tmp = real(trace(tmp_2));
      lrqtop[cb] -= qtop_tmp;
          }
  
            
  /* Lattice version of S_ratio */
  FILL(rtmp,PI);
  rtmp = TO_REAL(WORD_VALUE(WORD_rtmp,ONE) / (16*rtmp*rtmp));
  lract[0] = lract[0] * rtmp;
  lract[1] = lract[1] * rtmp;
  
  /* Lattice version of qtop */
  FILL(rtmp,PI);
  rtmp = TO_REAL(WORD_VALUE(WORD_rtmp,ONE) / (256*rtmp*rtmp));
  lrqtop[0] = lrqtop[0] * rtmp;
  lrqtop[1] = lrqtop[1] * rtmp;
  
  END_CODE();
}
