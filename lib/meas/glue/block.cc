// $Id: block.cc,v 1.1 2005-02-17 02:52:11 edwards Exp $

#error "NOT FULLY CONVERTED"

/* Construct block links from: */

/*      x           x                       x------x */
/*      |           |           _______            | */
/*      |           |           \                  | */
/*      |           |            \                 | */
/*      |     =     x     +       \                x */
/*      |           |             /                | */
/*      |           |            /                 | */
/*      |           |           /______            | */
/*      x           x                       x------x */

/* projected back onto SU(Nc) */

/* Warning: this works only for Nc = 2 and 3 ! */

/* u -- gauge field ( Read ) */
/* u_block -- blocked gauge field ( Write ) */
/* cb -- checkerboard of blocked gauge field ( Read ) */
/* mu -- direction of blocked gauge field ( Read ) */
/* bl_level -- blocking level (of the u's) ( Read ) */
/* BlkAccu -- accuracy in fuzzy link projection ( Read ) */
/* BlkMax -- maximum number of iterations in fuzzy link projection ( Read ) */
/* j_decay -- no staple in direction j_decay ( Read ) */

include(types.mh)

SUBROUTINE(block, u, u_block, cb, mu, bl_level, BlkAccu, BlkMax, j_decay)

multi1d<LatticeColorMatrix> u(Nd);
LatticeColorMatrix u_block;
int cb;
int mu;
int bl_level;
Real BlkAccu;
int BlkMax;
int j_decay;
{				/* Local Variables */
  include(COMMON_DECLARATIONS)
  LatticeColorMatrix u_dble;
  LatticeColorMatrix u_unproj;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_3;
  LatticeReal trace_tmp;
  LatticeBoolean bad;

  int nu;
  int cbs;
  int iter;
  int su2_index;
  int numbad;		/* not used here */
  int n_blk;
  int wrswitch;

  Double old_tr;
  Double new_tr;
  Double ddummy;
  Double conver;

  START_CODE();

  /* Determine checkboard of nearest neighbors at blocking level bl_level */
  if( bl_level == 0 )
    cbs = 1 - cb;
  else
    cbs = cb;

        
  /* Construct straight line segment x------x------x */
  /* tmp_1(x) = u(x+mu*2^bl_level,mu) */
  NEIGHBOUR(u(cbs,mu), tmp_1, cb, FORWARD, mu, bl_level);

  /* u_block = u(x,mu) * tmp_1 */
  u_block = u[mu][cb] * tmp_1;

  /* For blocking level zero, do it also for the other checkerboard */
  if( bl_level == 0 )
  {
    /* tmp_1(x) = u(x+mu*2^bl_level,mu) */
    NEIGHBOUR(u(cb,mu), tmp_1, cbs, FORWARD, mu, bl_level);

    /* u_dble = u(x,mu) * tmp_1 */
    u_dble = u[mu][cbs] * tmp_1;
  }
  /* For higher blocking levels, just _copy over */
  else
  {
    u_dble = u_block;
  }

  /* Now construct and add the staples, except in direction j_decay */
  for(nu = 0; nu < Nd; ++nu)
  {
    if( nu != mu && nu != j_decay )
    {

      /* Forward staple */
      /* tmp_1(x) = u(x+mu*2**(bl_level+1),nu) */
      NEIGHBOUR(u(cb,nu), tmp_1, cb, FORWARD, mu, bl_level+1);

      /* tmp_2(x) = u_dble(x+nu*2^bl_level) */
      NEIGHBOUR(u_dble, tmp_2, cb, FORWARD, nu, bl_level);

      /* tmp_3 = tmp_2 * tmp_1_dagger */
      tmp_3 = tmp_2 * adj(tmp_1);

      /* u_block = u_block + u(x,nu) * tmp_3 */
      u_block += u[nu][cb] * tmp_3;

      /* Backward staple */
      /* tmp_1(x) = u(x+mu*2^(bl_level+1),nu) */
      NEIGHBOUR(u(cbs,nu), tmp_1, cbs, FORWARD, mu, bl_level+1);

      /* tmp_2 = u_dble(x,mu) * tmp_1 */
      tmp_2 = u_dble * tmp_1;

      /* tmp_1 = u_dagger(x,nu) * tmp_2 */
      tmp_1 = adj(u[nu][cbs]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-nu*2^bl_level) */
      NEIGHBOUR(tmp_1, tmp_2, cb, BACKWARD, nu, bl_level);

      /* u_block = u_block + tmp_2 */
      u_block += tmp_2;

    }
  }

        
  /* Now project back to SU(3) by maximizing tr(u_block*u_dble_dagger), */
  /* where u_dble is the unprojected block link. */
  /* This is done by looping proj_iter times over the 3 SU(2) subgroups */
    u_unproj = adj(u_block);

  
  /* Start with a unitarized version */
  reunit (u_block, bad, OPTION[REUNITARIZE], numbad);

  /* The initial trace */
    trace_tmp = real(trace(u_block * u_unproj));
  old_tr = sum(trace_tmp);

  old_tr = old_tr / TO_DOUBLE(vol_cb*Nc);
  n_blk = 0;
  wrswitch = NO;		/* Write out iterations? */
  conver = 1;

  while ( conver > BlkAccu  &&  n_blk < BlkMax )
  {
    n_blk = n_blk + 1;

    /* Loop over SU(2) subgroup su2_index */
    for(su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
      su3proj (u_block, u_unproj, su2_index);
    
    /* Reunitarize */
    reunit (u_block, bad, OPTION[REUNITARIZE_LABEL], numbad);
    if ( numbad > 0 )
    {
      FPRINTF(trm_out,"BLOCK: WARNING unitarity violation\n");
      FPRINTF(trm_out,"   n_blk= %d  bl_level= %d\n",n_blk,bl_level);
      FPRINTF(trm_out,"   cb= %d  mu= %d  numbad= %d\n",cb,mu,numbad);
    }
    
    /* Calculate the trace */
    trace_tmp = real(trace(u_block * u_unproj));
    new_tr = sum(trace_tmp);
    new_tr = new_tr / TO_DOUBLE(vol_cb*Nc);

    if( wrswitch == YES )
      FPRINTF(trm_out," BLOCK: %4d  %16.9g%16.9g\n",n_blk,old_tr,new_tr);

    /* Normalized convergence criterion: */
    ddummy = new_tr;
    ddummy -= old_tr;
    ddummy = ddummy / old_tr;
    conver = fabs(ddummy);
    old_tr = new_tr;
  }

  if ( wrswitch == YES )
    push(xml_out,"Final_blkb");
write(xml_out, "bl_level", bl_level);
write(xml_out, "cb", cb);
write(xml_out, "mu", mu);
write(xml_out, "n_blk", n_blk);
write(xml_out, "new_tr", new_tr);
pop(xml_out);

      
  END_CODE();
}
