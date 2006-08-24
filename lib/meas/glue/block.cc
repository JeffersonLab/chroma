// $Id: block.cc,v 3.3 2006-08-24 03:18:41 edwards Exp $
/*! \file
 *  \brief Construct "block" links
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/shift2.h"
#include "util/gauge/su3proj.h"
#include "meas/glue/block.h"

namespace Chroma 
{ 

  //! Construct block links
  /*!
   * \ingroup glue
   *
   * Construct block links from:
   *      x           x                       x------x
   *      |           |           _______            |
   *      |           |           \                  |
   *      |           |            \                 |
   *      |     =     x     +       \                x
   *      |           |             /                |
   *      |           |            /                 |
   *      |           |           /______            |
   *      x           x                       x------x

   * projected back onto SU(Nc)

   * Warning: this works only for Nc = 2 and 3 !

   * \param u_block    blocked gauge field ( Write )
   * \param u          gauge field ( Read )
   * \param mu         direction of blocked gauge field ( Read )
   * \param bl_level   blocking level (of the u's) ( Read )
   * \param BlkAccu    accuracy in fuzzy link projection ( Read )
   * \param BlkMax     maximum number of iterations in fuzzy link projection ( Read )
   * \param j_decay    no staple in direction j_decay ( Read ) 
   */

  void block(LatticeColorMatrix& u_block, 
	     const multi1d<LatticeColorMatrix>& u, 
	     int mu, int bl_level, 
	     const Real& BlkAccu, int BlkMax, int j_decay)
  {
    START_CODE();

    LatticeColorMatrix u_dble;
    LatticeColorMatrix u_unproj;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;

    /* Construct straight line segment x------x------x */
    /* tmp_1(x) = u(x+mu*2^bl_level,mu) */
    tmp_1 = shift2(u[mu], FORWARD, mu, bl_level);

    /* u_block = u(x,mu) * tmp_1 */
    u_block = u[mu] * tmp_1;

    /* A copy */
    u_dble = u_block;

    /* Now construct and add the staples, except in direction j_decay */
    for(int nu = 0; nu < Nd; ++nu)
    {
      if( nu != mu && nu != j_decay )
      {
	/* Forward staple */
	/* tmp_1(x) = u(x+mu*2**(bl_level+1),nu) */
	tmp_1 = shift2(u[nu], FORWARD, mu, bl_level+1);

	/* tmp_2(x) = u_dble(x+nu*2^bl_level) */
	tmp_1 = shift2(u_dble, FORWARD, nu, bl_level);

	/* tmp_3 = tmp_2 * tmp_1_dagger */
	tmp_3 = tmp_2 * adj(tmp_1);

	/* u_block = u_block + u(x,nu) * tmp_3 */
	u_block += u[nu] * tmp_3;

	/* Backward staple */
	/* tmp_1(x) = u(x+mu*2^(bl_level+1),nu) */
	tmp_1 = shift2(u[nu], FORWARD, mu, bl_level+1);

	/* tmp_2 = u_dble(x,mu) * tmp_1 */
	tmp_2 = u_dble * tmp_1;

	/* tmp_1 = u_dagger(x,nu) * tmp_2 */
	tmp_1 = adj(u[nu]) * tmp_2;

	/* tmp_2(x) = tmp_1(x-nu*2^bl_level) */
	tmp_2 = shift2(tmp_1, BACKWARD, nu, bl_level);

	/* u_block = u_block + tmp_2 */
	u_block += tmp_2;
      }
    }

        
    /* Now project back to SU(3) by maximizing tr(u_block*u_dble_dagger), */
    /* where u_dble is the unprojected block link. */
    /* This is done by looping proj_iter times over the 3 SU(2) subgroups */
    u_unproj = adj(u_block);

    /* Start with a unitarized version */
    reunit(u_block);

    /* The initial trace */
    Double old_tr = sum(real(trace(u_block * u_unproj)));
    old_tr /= Double(QDP::Layout::vol()*Nc);

    int n_blk = 0;
    bool wrswitch = true;		/* Write out iterations? */
    Double conver = 1;

    while ( toBool(conver > BlkAccu)  &&  n_blk < BlkMax )
    {
      n_blk = n_blk + 1;

      /* Loop over SU(2) subgroup su2_index */
      for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	su3proj(u_block, u_unproj, su2_index);
    
      /* Reunitarize */
      {
	LatticeBoolean bad;
	int numbad;
	reunit(u_block, bad, numbad, REUNITARIZE_LABEL);
	if ( numbad > 0 )
	{
	  QDPIO::cout << "BLOCK: WARNING unitarity violation\n";
	  QDPIO::cout << "   n_blk= " << n_blk << "  bl_level= " << bl_level << "\n";
	  QDPIO::cout << "   mu= " << mu << "  numbad= " << numbad << endl;
	}
      }
    
      /* Calculate the trace */
      Double new_tr = sum(real(trace(u_block * u_unproj)));
      new_tr /= Double(QDP::Layout::vol()*Nc);

      if( wrswitch )
      {
	QDPIO::cout << "BLOCK: n=" << n_blk 
		    << " old_tr = " << old_tr << " new_tr = " << new_tr << endl;
      }

      /* Normalized convergence criterion: */
      conver = fabs((new_tr - old_tr) / old_tr);
      old_tr = new_tr;
    }

    END_CODE();
  }

}
