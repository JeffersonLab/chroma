// $Id: su3hb.cc,v 3.1 2007-02-22 21:11:49 bjoo Exp $
/*! \file
 *  \brief Do one SU(2) subgroup heatbath update of SU(Nc) matrix U with action
 */

#include "chromabase.h"
#include "update/heatbath/su3hb.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

namespace Chroma 
{

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
	     HeatbathType algorithm,
	     const Subset& sub)
  {
    START_CODE();
  
    /* V = U*W */
    LatticeColorMatrix v;
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
#if 0
    Real dummy = TO_REAL(BetaMC) / TO_REAL(Nc);
    if( AnisoP == YES ) dummy /= TO_REAL(xi_0);
#else
    QDPIO::cerr << "FIX Beta" << endl;
    Real dummy = 1;   // **WARNING** FIX THIS
#endif
    r_l[sub] *= dummy;

    int num_sites = sub.numSiteTable();
    int itrials = num_sites;
    ntrials = 0;
    nfails = 0;
  
    LatticeBoolean lupdate;
    LatticeReal lftmp1;
    LatticeReal lftmp2;
    LatticeBoolean lbtmp2;
  
    switch (algorithm)
    {
    case HEATBATH_TYPE_KPHB:
    {
      /*
       * Kennedy-Pendleton heatbath
       */
      int n_done = 0;
      int nhb = 0;

      r[0][sub] = a[0];
      lupdate[sub] = true;

      while ( nhb < nheat && n_done < num_sites )
      {
	ntrials += itrials;

	random(r[1], sub);
	r[1][sub] = log(r[1]);

	random(r[2], sub);
	r[2][sub] = log(r[2]);

	random(lftmp, sub);
	r[3][sub] = pow(cos(Real(twopi)*lftmp),2);

	r[1][sub] += r[2] * r[3];
	r[2][sub] = r[1] / r_l;

	/* r[2] is now a trial for 1+r[0] */
	/* see if this is accepted */
	random(lftmp, sub);
	r[1][sub] = lftmp*lftmp;

	lbtmp[sub]  = r[1] <= (1 + 0.5*r[2]);
	lbtmp2[sub] = lbtmp & lupdate;
	r[0][sub] = where(lbtmp2, 1+r[2], r[0]);

	int itmp = toInt(sum(where(lbtmp2, LatticeInteger(1), LatticeInteger(0)),sub));
	n_done += itmp;
	itrials -= itmp;
	nfails += itrials;

	lupdate &= ! lbtmp;

	++nhb;
      }
    }
    break;

    case HEATBATH_TYPE_CrHB:
    {
      /*+ */
      /* Creutz heatbath */
      /*- */
      int n_done = 0;
      int nhb = 0;

      lftmp2 = exp(2 * r_l) - 1;

      r[0] = a[0];
      lupdate = true;

      while ( nhb < nheat && n_done < num_sites )
      {
	ntrials += itrials;

	random(r[1], sub);
	r[2][sub] = log(1 + r[1] * lftmp2) / r_l  - 1;

	/* r[2] is now a trial for r[0] */
	/* see if this is accepted */
	random(lftmp1, sub);
	r[1][sub] = lftmp1 * lftmp1;

	lbtmp = r[1] <= (1 - r[2] * r[2]);
	lbtmp2 = lbtmp & lupdate;

	r[0][sub] = where(lbtmp2, r[2], r[0]);

	int itmp = toInt(sum(where(lbtmp2, LatticeInteger(1), LatticeInteger(0)),sub));
	n_done += itmp;
	itrials -= itmp;
	nfails += itrials;

	lupdate &= ! lbtmp;

	++nhb;
      }
    }
    break;

    default:
      QDPIO::cerr << __func__ << ": unknown algorithm type" << endl;
      QDP_abort(1);
    }
  
      
    /* Now create r[1], r[2] and r[3] according to the spherical measure */
    /* Take absolute value to guard against round-off */
    random(lftmp1, sub);
    r[2][sub] = 1 - 2*lftmp1;

    lftmp1[sub] = fabs(1 - r[0]*r[0]);
    r[3][sub] = -(sqrt(lftmp1) * r[2]);
  
    /* Take absolute value to guard against round-off */
    r_l[sub] = sqrt(fabs(lftmp1 - r[3]*r[3]));

    random(lftmp1, sub);
    lftmp1[sub] *= twopi;
    r[1][sub] = r_l * cos(lftmp1);
    r[2][sub] = r_l * sin(lftmp1);
  

    /* Update matrix is B = R * A, with B, R and A given by b_i, r_i and a_i */
    multi1d<LatticeReal> b(4);
    b[0][sub] = r[0]*a[0] - r[1]*a[1] - r[2]*a[2] - r[3]*a[3];
    b[1][sub] = r[0]*a[1] + r[1]*a[0] - r[2]*a[3] + r[3]*a[2];
    b[2][sub] = r[0]*a[2] + r[2]*a[0] - r[3]*a[1] + r[1]*a[3];
    b[3][sub] = r[0]*a[3] + r[3]*a[0] - r[1]*a[2] + r[2]*a[1];
  
    /*
     * Now fill an SU(3) matrix V with the SU(2) submatrix su2_index
     * paramtrized by a_k in the sigma matrix basis.
     */
    sunFill(v, b, su2_index, sub);

    // U = V*U
    LatticeColorMatrix tmp;
    tmp[sub] = v * u;
    u[sub] = tmp;

    END_CODE();
  }

}  // end namespace Chroma
