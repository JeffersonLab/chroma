// $Id: grelax.cc,v 1.4 2004-01-02 20:53:25 edwards Exp $
/*! \file
 *  \brief Perform a single gauge fixing iteration
 */

#include "chromabase.h"
#include "meas/gfix/axgauge.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

using namespace QDP;

/********************** HACK ******************************/
// Primitive way for now to indicate the time direction
static int tDir() {return Nd-1;}
static Real xi_0() {return 1.0;}
static bool anisoP() {return false;}
/******************** END HACK ***************************/

//! Perform a single gauge fixing iteration
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param ug         (gauge fixed) gauge field ( Modify )
 * \param u_neg      (gauge fixed) gauge field, negative links ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read ) 
 */

void grelax(multi1d<LatticeColorMatrix>& ug, 
	    const multi1d<LatticeColorMatrix>& u_neg,
	    int j_decay, int su2_index, bool ordo,
	    const Real& orpara)
{
  LatticeColorMatrix v;
  multi1d<LatticeReal> r(4);
  multi1d<LatticeReal> a(4);
  
  START_CODE("grelax");
      
  /* Sum links to be gauge transformed on site x not in the direction */
  /* j_decay into matrix V: */
  if (tDir() != j_decay)
  {
    v = ug[tDir()] + adj(u_neg[tDir()]);

    if (anisoP())
      v *= pow(xi_0(), 2);
  }
  else
  {
    v = 0;
  }

  for(int mu = 0; mu < Nd; ++mu)
  {
    if (mu != j_decay && mu != tDir())
      v += ug[mu] + adj(u_neg[mu]);
  }
  
  if (Nc > 1)
  {
                                    
    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the SU(Nc) matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    su2Extract(r, v, su2_index, all);
  
    /*
     * Now project onto SU(2)
     */
    LatticeReal r_l = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  
    // Normalize
    LatticeBoolean lbtmp = r_l > fuzz;
    LatticeReal lftmp = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
    a[2] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
    a[3] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      /* get angle */
      LatticeReal theta_old = acos(a[0]);

      /* old sin */
      LatticeReal oldsin = sin(theta_old);

      /* overrelax, i.e. multiply by the angle */
      LatticeReal theta_new = theta_old * orpara;

      /* compute sin(new)/sin(old) */
      /* set the ratio to 0, if sin(old) < FUZZ */
      lftmp = where(oldsin > fuzz, sin(theta_new) / oldsin, LatticeReal(0));

      /* get the new cos = a[0] */
      a[0] = cos(theta_new);

      /* get the new a_k, k = 1, 2, 3 */
      a[1] *= lftmp;
      a[2] *= lftmp;
      a[3] *= lftmp;
    }
  
          
    /* Now fill the SU(Nc) matrix V with the SU(2) submatrix 'su2_index' */
    /* paramtrized by a_k in the sigma matrix basis. */
    sunFill(v, a, su2_index, all);

  }
  else		/* Nc = 1 */
  {
    r[0] = real(peekColor(v, 0, 0));
    r[1] = imag(peekColor(v, 0, 0));
    LatticeReal r_l = sqrt(r[0] * r[0] + r[1] * r[1]);

    // Normalize
    LatticeBoolean lbtmp = r_l > fuzz;
    LatticeReal lftmp = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      Real pi = 0.5 * twopi;

      /* get angle */
      LatticeReal theta = acos(a[0]) + pi;     // Do I really want to add pi ??? This was in szin

      /* overrelax, i.e. multiply by the angle */
      theta *= orpara;

      /* get the new cos = a[0] */
      a[0] = cos(theta);

      /* get the new sin */
      a[1] = sin(theta);
    }

    pokeColor(v, cmplx(a[0],a[1]), 0, 0);
  }

    
  /* Now do the gauge transformation with matrix V:  */
  LatticeColorMatrix u_tmp;
  for(int mu = 0; mu < Nd; ++mu)
  {
    /* Forward link (ug(x,mu) = v(x)*ug(x,mu)) */
    u_tmp = v * ug[mu];
    ug[mu] = u_tmp;

    /* Backward link ug(x,mu) = u_neg(x+mu,mu)*v_dag(x+mu) */
    u_tmp = u_neg[mu] * adj(v);
    ug[mu] = shift(u_tmp, FORWARD, mu);
  }
  
  END_CODE("grelax");
}
