// $Id: grelax.cc,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Perform a single gauge fixing iteration
 */

#include "chromabase.h"
#include "meas/gfix/axgauge.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

namespace Chroma {


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
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param u          original gauge field ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read ) 
 */

void grelax(LatticeColorMatrix& g,
	    const multi1d<LatticeColorMatrix>& u, 
	    int j_decay, int su2_index, int cb, bool ordo,
	    const Real& orpara)
{
  LatticeColorMatrix v;
  multi1d<LatticeReal> r(4);
  multi1d<LatticeReal> a(4);
  
  START_CODE();
      
  // Rotate the matrices using the current gauge rotation
  multi1d<LatticeColorMatrix> ug(Nd);
  for(int mu = 0; mu < Nd; ++mu)
    ug[mu] = g * u[mu] * shift(adj(g), FORWARD, mu);
    
  /* Gather the Nd negative links attached to a site: */
  /* u_tmp(x,mu) = U(x-mu,mu) */
  multi1d<LatticeColorMatrix> u_neg(Nd);
  for(int mu=0; mu<Nd; ++mu)
    u_neg[mu][rb[cb]] = shift(ug[mu], BACKWARD, mu);
  
  /* Sum links to be gauge transformed on site x not in the direction */
  /* j_decay into matrix V: */
  v = 0;
  if (tDir() != j_decay)
  {
    v[rb[cb]] = ug[tDir()] + adj(u_neg[tDir()]);

    if (anisoP())
      v[rb[cb]] *= pow(xi_0(), 2);
  }

  for(int mu = 0; mu < Nd; ++mu)
  {
    if (mu != j_decay && mu != tDir())
      v[rb[cb]] += ug[mu] + adj(u_neg[mu]);
  }
  
  if (Nc > 1)
  {
                                    
    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the SU(Nc) matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    su2Extract(r, v, su2_index, rb[cb]);
  
    /*
     * Now project onto SU(2)
     */
    LatticeReal r_l;
    r_l[rb[cb]] = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  
    // Normalize
    LatticeBoolean lbtmp;
    lbtmp[rb[cb]] = r_l > fuzz;
    LatticeReal lftmp;
    lftmp[rb[cb]] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0][rb[cb]] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1][rb[cb]] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
    a[2][rb[cb]] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
    a[3][rb[cb]] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      /* get angle */
      LatticeReal theta_old;
      theta_old[rb[cb]] = acos(a[0]);

      /* old sin */
      LatticeReal oldsin;
      oldsin[rb[cb]] = sin(theta_old);

      /* overrelax, i.e. multiply by the angle */
      LatticeReal theta_new;
      theta_new[rb[cb]] = theta_old * orpara;

      /* compute sin(new)/sin(old) */
      /* set the ratio to 0, if sin(old) < FUZZ */
      lftmp[rb[cb]] = where(oldsin > fuzz, sin(theta_new) / oldsin, LatticeReal(0));

      /* get the new cos = a[0] */
      a[0][rb[cb]] = cos(theta_new);

      /* get the new a_k, k = 1, 2, 3 */
      a[1][rb[cb]] *= lftmp;
      a[2][rb[cb]] *= lftmp;
      a[3][rb[cb]] *= lftmp;
    }
  
          
    /* Now fill the SU(Nc) matrix V with the SU(2) submatrix 'su2_index' */
    /* paramtrized by a_k in the sigma matrix basis. */
    sunFill(v, a, su2_index, rb[cb]);

  }
  else		/* Nc = 1 */
  {
    r[0][rb[cb]] = real(peekColor(v, 0, 0));
    r[1][rb[cb]] = imag(peekColor(v, 0, 0));
    LatticeReal r_l;
    r_l[rb[cb]] = sqrt(r[0] * r[0] + r[1] * r[1]);

    // Normalize
    LatticeBoolean lbtmp;
    lbtmp[rb[cb]] = r_l > fuzz;
    LatticeReal lftmp;
    lftmp[rb[cb]] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0][rb[cb]] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1][rb[cb]] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      Real pi = 0.5 * twopi;

      /* get angle */
      LatticeReal theta;
      theta[rb[cb]] = acos(a[0]) + pi;     // Do I really want to add pi ??? This was in szin

      /* overrelax, i.e. multiply by the angle */
      theta[rb[cb]] *= orpara;

      /* get the new cos = a[0] */
      a[0][rb[cb]] = cos(theta);

      /* get the new sin */
      a[1][rb[cb]] = sin(theta);
    }

    pokeColor(v, cmplx(a[0],a[1]), 0, 0);
  }

    
  // Now do the gauge transformation with matrix V:
  // Multiply into the global "g" field
  LatticeColorMatrix u_tmp;
  u_tmp[rb[cb]] = v * g;
  g[rb[cb]] = u_tmp;
  
  END_CODE();
}

}  // end namespace Chroma
