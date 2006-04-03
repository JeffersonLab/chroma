// $Id: meslate.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $

#error "NOT FULLY CONVERTED"


/* MESLATE - measure the lattice energy and the naive topological charge,  */
/*           with the "clover" definition */
/*           of the naive topological charge density. It also computes the */
/*           ratio of the action to the continuum instaton action. */

/* u -- gauge field (Read) */
/* i_cool  -- cooling sweep number (used for printing out headers) (Read) */
/* t_slice -- time slice number (Read) */
/* qtop    -- topological charge (Write) */
/* S_ratio -- action to continuum instanton action (Write) */
include(types.mh)

SUBROUTINE(meslate, u, i_cool, j_decay, t_slice, qtop, S_ratio)


multi1d<LatticeColorMatrix> u(Nd);
int i_cool;
int j_decay;
int t_slice;
Double qtop;
Double S_ratio;

{ /* Local Variables */
  include(COMMON_DECLARATIONS)
  /* local stuff for printing lattice action */
  int ix;
  int iy;
  int iz;
  int it;
  int x;
  int y;
  int z;
  int t;
  int lx;
  int ly;
  int lz;
  int lt;
  int i;
  multi1d<int> coord(Nd);
  LatticeReal lract;
  LatticeReal lrqtop;

  /* other local variables */
  LatticeColorMatrix u_clov1;
  LatticeColorMatrix u_clov2;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_3;
  LatticeReal qtop_tmp;
  LatticeReal plaq_tmp;

  multi1d<Double> action_slice(lt);
  multi1d<Double> qtop_slice(lt);
  Double plaq;
  Real rtmp;
  Double tmp;
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

  if ( Nd != 4 )
    QDP_error_exit("code requires 4 dimensions for printing[only], Nd);

  if ( j_decay != Nd-1 && j_decay > 0 )
    QDP_error_exit("code requires j_decay == Nd-1 and > 0", j_decay);

  lx = nrow[0];
  ly = nrow[1];
  lz = nrow[2];
  lt = nrow[3];

                  
  lract = 0;
  lrqtop = 0;
  qtop = 0;
  plaq = 0;

  /* Lattice version of S_ratio */
  rtmp = TO_REAL(2*Nd*(Nd-1)*Nc);
  FILL(lract(0), rtmp);
  FILL(lract(1), rtmp);

  /* Loop over checkerboards and triplet of perpendicular planes */
  mu1 = 0;
  for(cb = 0;cb  <= ( 1); ++cb )
    for(nu1 = 1;nu1  < ( Nd); ++nu1 )
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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
      plaq += sum(plaq_tmp);
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

          
  plaq = plaq / TO_DOUBLE(2*(vol)*Nd*(Nd-1)*Nc);

  /* Lattice version of S_ratio */
  FILL(rtmp,PI);
  rtmp = TO_REAL(WORD_VALUE(WORD_rtmp,ONE) / (16*rtmp*rtmp));
  lract[0] = lract[0] * rtmp;
  lract[1] = lract[1] * rtmp;

  /* Check */
  S_ratio = sum(lract[0]);
  S_ratio += sum(lract[1]);

  /* Lattice version of qtop */
  FILL(rtmp,PI);
  rtmp = TO_REAL(WORD_VALUE(WORD_rtmp,ONE) / (256*rtmp*rtmp));
  lrqtop[0] = lrqtop[0] * rtmp;
  lrqtop[1] = lrqtop[1] * rtmp;

  /* Check */
  qtop = sum(lrqtop[0]);
  qtop += sum(lrqtop[1]);

  /* sum the action density over each time slice */
  action_slice = sumMulti(lract[0], timeslice);
  action_slice += sumMulti(lract[1], timeslice);

  /* sum the topological charge density over each time slice */
  qtop_slice = sumMulti(lrqtop[0], timeslice);
  qtop_slice += sumMulti(lrqtop[1], timeslice);


  /*### Print out the action and topological charge density in SCIAN stf format */
  
  fprintf(unit_21," TIME %4d\n RANK 3\n DIMENSIONS %4d%4d%4d\n\
 BOUNDS 0 %4d 0 %4d 0 %4d\n SCALAR\n ORDER COLUMN\n INTERLACED\n DATA\n",
	  i_cool,lx+1,ly+1,lz+1,lx,ly,lz);
  fprintf(unit_22," TIME %4d\n RANK 3\n DIMENSIONS %4d%4d%4d\n\
 BOUNDS 0 %4d 0 %4d 0 %4d\n SCALAR\n ORDER COLUMN\n INTERLACED\n DATA\n",
	  i_cool,lx+1,ly+1,lz+1,lx,ly,lz);

  /* Extend the x,y,z axes for padding. Done for Scian and periodic BC. */
  for(z = 0;z  <= ( lz); ++z )
    for(y = 0;y  <= ( ly); ++y )
      for(x = 0;x  <= ( lx); ++x )
      {
	ix = INTEGER_MOD_FUNCTION(x,lx);
	iy = INTEGER_MOD_FUNCTION(y,ly);
	iz = INTEGER_MOD_FUNCTION(z,lz);
	it = INTEGER_MOD_FUNCTION(t_slice,lt);

	cb = INTEGER_MOD_FUNCTION(ix+iy+iz+it,2);
	i  = (int)(ix/2);

	coord[0] = i;
	coord[1] = iy;
	coord[2] = iz;
	coord[3] = it;

	INDEXING(lract(cb),coord,rtmp);
	fprintf(unit_21," %14.7g\n",rtmp);

	INDEXING(lrqtop(cb),coord,rtmp);
	fprintf(unit_22," %14.7g\n",rtmp);
      }

  fprintf(unit_21," END\n");
  fprintf(unit_22," END\n");

    /*### END ############### */


  /*### Print out the action and topological charge density in SCIAN stf format */
  
  lx = cb_nrow[0];


  fprintf(unit_23," TIME %4d\n RANK 4\n DIMENSIONS %4d%4d%4d%4d\n\
 BOUNDS 0 %4d 0 %4d 0 %4d 0 %4d\n SCALAR\n ORDER COLUMN\n INTERLACED\n DATA\n",
	  i_cool,lx+1,ly+1,lz+1,lt+1,lx,ly,lz,lt);


  /* Extend the x,y,z axes for padding. Done for Scian and periodic BC. */
  for(t = 0;t  <= ( lt); ++t )
    for(z = 0;z  <= ( lz); ++z )
      for(y = 0;y  <= ( ly); ++y )
	for(x = 0;x  <= ( lx); ++x )
	{
	  ix = INTEGER_MOD_FUNCTION(x,lx);
	  iy = INTEGER_MOD_FUNCTION(y,ly);
	  iz = INTEGER_MOD_FUNCTION(z,lz);
	  it = INTEGER_MOD_FUNCTION(t,lt);

	  coord[0] = ix;
	  coord[1] = iy;
	  coord[2] = iz;
	  coord[3] = it;

	  INDEXING(lract(0),coord,rtmp);
	  fprintf(unit_23," %14.7g\n",rtmp);
	}

  fprintf(unit_23," END\n");

    /*### END ############### */

  printf(" %d %d %g %g %g\n",i_cool,t_slice,plaq,qtop,S_ratio);
  push(xml_out,"MesLatE");
write(xml_out, "i_cool", i_cool);
write(xml_out, "t_slice", t_slice);
write(xml_out, "plaq", plaq);
write(xml_out, "qtop", qtop);
write(xml_out, "S_ratio", S_ratio);
write(xml_out, "qtop_slice", qtop_slice);
write(xml_out, "action_slice", action_slice);
pop(xml_out);

        
  END_CODE();
}
