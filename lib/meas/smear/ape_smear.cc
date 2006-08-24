//  $Id: ape_smear.cc,v 3.2 2006-08-24 02:34:18 edwards Exp $
/*! \file
 *  \brief APE-smearing of the gauge configuration
 */

#include "chromabase.h"
#include "meas/smear/ape_smear.h"
#include "util/gauge/reunit.h"
#include "util/gauge/su3proj.h"
#include "util/gauge/shift2.h"

namespace Chroma 
{ 
  //! Construct APE smeared links from:
  /*!
   * \ingroup smear
   *
   *
   *                                  _______ 
   *      x               x           \           x------x
   *      +               |            \                 |
   *      +     =     c * |     +       \                |
   *      +               |             /                |
   *      x               x            /          x------x
   *                                  /______ 
   *
   *
   * Arguments:
   *
   * where c is the smearing factor sm_fact, and projected back onto SU(Nc)

   * Warning: this works only for Nc = 2 and 3 !

   * \param u        gauge field ( Read )
   * \param u_smear  smeared gauge field ( Write )
   * \param mu       direction of smeared gauge field ( Read )
   * \param bl_level blocking level (of the u's) ( Read )
   * \param sm_fact  smearing factor ( Read )
   * \param BlkAccu  accuracy in fuzzy link projection ( Read )
   * \param BlkMax   maximum number of iterations in fuzzy link projection ( Read )
   * \param j_decay  no staple in direction j_decay ( Read )
   */

  void APE_Smear(const multi1d<LatticeColorMatrix>& u,
		 LatticeColorMatrix& u_smear,
		 int mu, int bl_level, 
		 const Real& sm_fact, const Real& BlkAccu, 
		 int BlkMax, int j_decay)
  {
    START_CODE();
  
    // Initialize smeared link: sm_fact * "old" link
    u_smear = u[mu] * sm_fact;

    // Now construct and add the staples, except in direction j_decay
    for(int nu = 0; nu < Nd; ++nu)
    {
      if( nu != mu && nu != j_decay )
      {
	/* Forward staple */
	/* u_smear = u_smear + u(x,nu) * u(x+nu*2^bl_level,mu) * adj(u(x+mu*2^bl_level,nu)) */
	u_smear += u[nu] * shift2(u[mu], FORWARD, nu, bl_level) 
	  * adj(shift2(u[nu], FORWARD, mu, bl_level));

	/* Backward staple */
	/* tmp_1 = u_dagger(x,nu) * u(x,mu) * u(x+mu*2^bl_level,nu) */
	LatticeColorMatrix tmp_1 = adj(u[nu]) * u[mu] * shift2(u[nu], FORWARD, mu, bl_level);

	/* u_smear = u_smear + tmp_1(x-nu*2^bl_level) */
	u_smear += shift2(tmp_1, BACKWARD, nu, bl_level);
      }
    }

      
    /* Now project back to SU(3) by maximizing tr(u_smear*u_unproj_dagger), */
    /* where u_unproj is the unprojected smear link. */
    /* This is done by looping proj_iter times over the 3 SU(2) subgroups */
    LatticeColorMatrix u_unproj = adj(u_smear);

#if 0
    /* Start with a unitarized version */
    reunit(u_smear);
#else
    /* Start with original link */
    u_smear = u[mu];
#endif


#if 0
    // Do not yet support schroedinger functional BC
    if (SchrFun > 0)
    {
      /* Make it easy, since these are overwritten anyway! */
      copymask(u_smear,  lSFmask[mu], LatticeReal(1));
      copymask(u_unproj, lSFmask[mu], LatticeReal(1));
    }
#endif

    /* The initial trace */
    Double old_tr = sum(real(trace(u_smear * u_unproj))) / toDouble(Layout::vol()*Nc);
    Double new_tr;

    int n_smr = 0;
    bool wrswitch = false;			/* Write out iterations? */
    Double conver = 1;
  
    while ( toBool(conver > BlkAccu)  &&  n_smr < BlkMax )
    {
      ++n_smr;

      // Loop over SU(2) subgroup index
      for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	su3proj(u_smear, u_unproj, su2_index);

      /* Reunitarize */
      reunit(u_smear);
    
      /* Calculate the trace */
      new_tr = sum(real(trace(u_smear * u_unproj))) / toDouble(Layout::vol()*Nc);

      if( wrswitch )
	QDPIO::cout << " BLOCK: " << n_smr << " old_tr= " << old_tr << " new_tr= " << new_tr;

      /* Normalized convergence criterion: */
      conver = fabs((new_tr - old_tr) / old_tr);
      old_tr = new_tr;
    }
  

    //  if ( wrswitch )
    //  {
    //    push(nml_out,"Final_smear");
    //    write(nml_out, "mu", mu);
    //    write(nml_out, "n_smr", n_smr);
    //    write(nml_out, "new_tr", new_tr);
    //    pop(nml_out);
    //  }

#if 0
    // Do not yet support schroedinger functional BC
    if (SchrFun > 0)
    {
      /* Now do the overwrite. */
      copymask(u_smear, lSFmask[mu], SFBndFld[mu]);
    }
#endif

    END_CODE();
  }

} // Namespace Chroma
