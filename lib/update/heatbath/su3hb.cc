// $Id: su3hb.cc,v 1.3 2004-07-28 02:38:05 edwards Exp $
/*! \file
 *  \brief Do one SU(2) subgroup heatbath update of SU(Nc) matrix U with action
 */

#error "Not tested (or even compiled). However, reasonably well converted."

#error "NEED TO ADD NARROWING TO SUBSETS IN ALL OPERANDS"

#include "chromabase.h"
#include "update/heatbath/su3hb.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

using namespace QDP;

//! Do one SU(2) subgroup heatbath update of SU(Nc) matrix U with action
/*!
 * \ingroup heatbath
 *
 * BetaMC * [tr(U*W) + hc] / (2Nc).
 * Do at most nheat update trials of the Kennedy-Pendleton or Creutz SU(2)
 * heatbath algorithm.
 *
 * \param u          field to be updated ( Modify )
 * \param w          "staple" field in the action ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param nheat      maximal number of update trials ( Read )
 * \param ntrials    total number of link trials ( Write )
 * \param nfails     total number of failed trials ( Write ) 
 * \param sub        Subset for updating ( Read )
 */

void su3hb(LatticeColorMatrix& u,
	   const LatticeColorMatrix& w,
	   int su2_index,
	   int nheat,
	   int& ntrials,
	   int& nfails,
	   const OrderedSubset& sub)
{
  START_CODE();
  
  LatticeColorMatrix v;
  LatticeReal lftmp1;
  LatticeReal lftmp2;
  LatticeBoolean lbtmp2;
  LatticeBoolean lupdate;
  int nhb;
  int n_done;
  int itrials;
  int itmp;
  
  /* V = U*W */
  v[sub] = u * w;
  
  /* Extract components r_k proportional to SU(2) submatrix su2_index */
  /* from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the */
  /* sigma matrix basis. */
  multi1d<LatticeReal> r(4);
  su2Extract(r, v, su2_index, sub);
  
  /*
   * Now project onto SU(2)
   */
  LatticeReal r_l;
  r_l[sub] = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);

  // Normalize
  LatticeBoolean lbtmp;
  lbtmp[sub] = r_l > fuzz;
  LatticeReal lftmp;
  lftmp[sub] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

  // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
  //  and   (1,0,0,0)  for sites with r_l < fuzz
  multi1d<LatticeReal> a(4);
  a[0][sub] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
  a[1][sub] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
  a[2][sub] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
  a[3][sub] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));

  /* Do up to nheat heatbath trials to create new r_0 */
  /* This is the only place BetaMC is used, so on AnisoP:  */
  Real dummy = TO_REAL(BetaMC) / TO_REAL(Nc);

#if 0
  if( AnisoP == YES ) dummy /= TO_REAL(xi_0);
#endif

  r_l[sub] *= dummy;
  itrials = vol_cb;
  ntrials = 0;
  nfails = 0;
  
  if ( Algorithm == KPHB )
  {
    /*
     * Kennedy-Pendleton heatbath
     */
    n_done = 0;
    nhb = 0;

    r[0][sub] = a_0;
    lupdate[sub] = true;

    while ( nhb < nheat && n_done < vol_cb )
    {
      ntrials += itrials;

      random(r[1], sub);
      r[1][sub] = log(r[1]);

      random(r[2], sub);
      r[2][sub] = log(r[2]);

      random(r[3]);
      dummy = 8*atan(1);
      r[3] = r[3] * dummy;
      r[3] = cos(r[3]);
      r[3] = r[3] * r[3];

      r[1] += r[2] * r[3];
      r[2] = r[1] / r_l;

      /* r[2] is now a trial for 1+r[0] */
      /* see if this is accepted */
      random(r[1]);
      r[1] = r[1] * r[1];
      lftmp1 = 1;
      dummy = 0.5;
      lftmp1 += r[2] * dummy;

      lbtmp = r[1] <= lftmp1;
      lbtmp2 = lbtmp & lupdate;

      lftmp1 = 1;
      lftmp1 += r[2];
      copymask(r[0], lbtmp2, lftmp1, REPLACE);

      itmp = sum(lbtmp2);
      n_done += itmp;
      itrials -= itmp;
      nfails += itrials;

      lbtmp = !lbtmp;
      lupdate = lupdate & lbtmp;

      nhb = nhb + 1;
    }
  }
  else
  {
    /*+ */
    /* Creutz heatbath */
    /*- */
    n_done = 0;
    nhb = 0;

    lftmp1 = 1;
    dummy = 2;
    lftmp2 = r_l * dummy;
    lftmp2 = exp(lftmp2);
    lftmp2 -= lftmp1;

    r[0] = a_0;
    FILL(lupdate, TRUE);

    while ( nhb < nheat && n_done < vol_cb )
    {
      ntrials += itrials;

      lftmp1 = 1;
      random(r[1]);
      lftmp1 += r[1] * lftmp2;
      lftmp1 = log(lftmp1);
      r[2] = lftmp1 / r_l;

      lftmp1 = 1;
      r[2] -= lftmp1;

      /* r[2] is now a trial for r[0] */
      /* see if this is accepted */
      random(r[1]);
      r[1] = r[1] * r[1];
      lftmp1 -= r[2] * r[2];

      lbtmp = r[1] <= lftmp1;
      lbtmp2 = lbtmp & lupdate;

      copymask(r[0], lbtmp2, r[2], REPLACE);

      itmp = sum(lbtmp2);
      n_done += itmp;
      itrials -= itmp;
      nfails += itrials;

      lbtmp = !lbtmp;
      lupdate = lupdate & lbtmp;

      nhb = nhb + 1;
    }
  }
  
      
  /* Now create r[1], r[2] and r[3] according to the spherical measure */
  /* Take absolute value to guard against round-off */
  r_l = fabs(1 - r[0]*r[0]);
  lftmp1[sub] = sqrt(r_l);
  random(lftmp2, sub);
  r[2] = 1;
  dummy = 2;
  r[2] -= lftmp2 * dummy;
  r[3] = -(lftmp1 * r[2]);
  
  r_l -= r[3] * r[3];
  /* Take absolute value to guard against round-off */
  r_l = fabs(r_l);
  r_l = sqrt(r_l);
  random(lftmp1);
  dummy = 8*atan(1);
  lftmp2 = lftmp1 * dummy;
  lftmp1 = cos(lftmp2);
  r[1] = r_l * lftmp1;
  lftmp1 = sin(lftmp2);
  r[2] = r_l * lftmp1;
  
        
  /* Update matrix is B = R * A, with B, R and A given by b_i, r_i and a_i */
  multi1d<LatticeReal> b(4);
  b[0] = r[0]*a[0] - r[1]*a[1] - r[2]*a[2] - r[3]*a[3];
  b[1] = r[0]*a[1] + r[1]*a[0] - r[2]*a[3] + r[3]*a[2];
  b[2] = r[0]*a[2] + r[2]*a[0] - r[3]*a[1] + r[1]*a[3];
  b[3] = r[0]*a[3] + r[3]*a[0] - r[1]*a[2] + r[2]*a[1];
  
  /*
   * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
   * paramtrized by a_k in the sigma matrix basis.
   */
  sunFill(v, b, su2_index, sub);

  // U = V*U
  LatticeColorMatrix tmp[sub] = v * u;
  u[sub] = tmp;

  END_CODE();
}
