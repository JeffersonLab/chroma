// $Id: grelax.cc,v 1.1 2003-12-06 20:56:56 edwards Exp $
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

void grelax(multi1d<LatticeColorMatrix> ug, 
	    const multi1d<LatticeColorMatrix> u_neg,
	    int j_decay, int su2_index, bool ordo,
	    const Real& orpara)
{
  LatticeColorMatrix u_tmp;
  LatticeColorMatrix v;
  multi2d<LatticeComplex> vv(Nc, Nc);
  multi1d<LatticeReal> r(4);
  multi1d<LatticeReal> a(4);
  LatticeReal r_l;
  LatticeReal lftmp1;
  LatticeReal lftmp2;
  
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
  
    /* Normalize */
             
    lftmp1 = 1;
    LatticeBoolean lbtmp = r_l > fuzz;
    copymask(lftmp1, lbtmp, r_l);
  
    lftmp2 = 1;
    lftmp1 = lftmp2 / lftmp1;
    a[0] = r[0] * lftmp1;
    a[1] = -(r[1] * lftmp1);
    a[2] = -(r[2] * lftmp1);
    a[3] = -(r[3] * lftmp1);
  
    /* Fill with (1,0,0,0) the sites with r_l < FUZZ */
    lbtmp = !lbtmp;
    lftmp1 = 0;
    copymask(a[0], lbtmp, lftmp2);
    copymask(a[1], lbtmp, lftmp1);
    copymask(a[2], lbtmp, lftmp1);
    copymask(a[3], lbtmp, lftmp1);
  
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      /* get angle */
      r[0] = acos(a[0]);

      /* get the old sin */
      r_l = sin(r[0]);

      /* overrelax, i.e. multiply by the angle */
      r[1] = r[0] * orpara;

      /* get the new cos = a[0] */
      a[0] = cos(r[1]);

      /* get the new sin */
      r[0] = sin(r[1]);

      /* compute sin(new)/sin(old) */
      lftmp1 = 1;
      lbtmp = r_l > fuzz;
      copymask(lftmp1, lbtmp, r_l);

      r[1] = r[0] / lftmp1;

      /* set the ratio to 0, if sin(old) < FUZZ */
      lbtmp = !lbtmp;
      lftmp1 = 0;
      copymask(r[1], lbtmp, lftmp1);

      /* get the new a_k, k = 1, 2, 3 */
      a[1] *= r[1];
      a[2] *= r[1];
      a[3] *= r[1];
    }
  
          
    /* Now fill the SU(Nc) matrix V with the SU(2) submatrix 'su2_index' */
    /* paramtrized by a_k in the sigma matrix basis. */
    sunFill(v, a, su2_index, all);

  }
  else		/* Nc = 1 */
  {
    for(int i=0; i < Nc; ++i)
      for(int j=0; j < Nc; ++j)
	vv[i][j] = peekColor(v, i, j);

    r[0] = real(vv[0][0]);
    r[1] = imag(vv[0][0]);
    r_l = r[0] * r[0];
    r_l += r[1] * r[1];
    r_l = sqrt(r_l);

    /* Normalize */
            
    lftmp1 = 1;
    LatticeBoolean lbtmp = r_l > fuzz;
    copymask(lftmp1, lbtmp, r_l);

    lftmp2 = 1;
    lftmp1 = lftmp2 / lftmp1;
    a[0] = r[0] * lftmp1;
    a[1] = -(r[1] * lftmp1);

    /* Fill with (1,0) the sites with r_l < FUZZ */
    lbtmp = !lbtmp;
    lftmp1 = 0;
    copymask(a[0], lbtmp, lftmp2);
    copymask(a[1], lbtmp, lftmp1);

    
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      /* get angle */
      r[0] = acos(a[0]);

//      lbtmp = a[1] < 0;    ?????? this was in SZIN, do not know why...
      lftmp1 = 0.5 * twopi;
      r[0] += lftmp1;

      /* overrelax, i.e. multiply by the angle */
      r[1] = r[0] * orpara;

      /* get the new cos = a[0] */
      a[0] = cos(r[1]);

      /* get the new sin */
      a[1] = sin(r[1]);
    }

    vv[0][0] = cmplx(a[0],a[1]);
    for(int i=0; i < Nc; ++i)
      for(int j=0; j < Nc; ++j)
	pokeColor(v, vv[i][j], i, j);
  }

    
  /* Now do the gauge transformation with matrix V:  */
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
