//  $Id: hyp_smear.cc,v 3.1 2006-08-11 16:13:30 edwards Exp $
/*! \file
 *  \brief Hyp-smearing of the gauge configuration
 */

#include "chromabase.h"
#include "meas/smear/hyp_smear.h"
#include "util/gauge/sun_proj.h"

namespace Chroma 
{ 
  //! Construct the "hyp-smeared" links of Anna Hasenfratz
  /*!
   * \ingroup smear
   *
   * Construct the "hyp-smeared" links of Anna Hasenfratz, with
   * staple coefficients alpha1, alpha2 and alpha3
   *
   * Arguments:
   *
   *  \param u		gauge field (Read)
   *  \param u_hyp	"hyp-smeared" gauge field (Write)
   *  \param alpha1	staple coefficient "1" (Read)
   *  \param alpha2	staple coefficient "2" (Read)
   *  \param alpha3	staple coefficient "3" (Read)
   *  \param BlkAccu	accuracy in SU(Nc) projection (Read)
   *  \param BlkMax	max number of iterations in SU(Nc) projection (Read)
   */

  void Hyp_Smear(const multi1d<LatticeColorMatrix>& u,
		 multi1d<LatticeColorMatrix>& u_hyp,
		 const Real& alpha1, const Real& alpha2, const Real& alpha3,
		 const Real& BlkAccu, int BlkMax)
  {
    multi1d<LatticeColorMatrix> u_lv1(Nd*(Nd-1));
    multi1d<LatticeColorMatrix> u_lv2(Nd*(Nd-1));
    LatticeColorMatrix u_tmp;
    LatticeColorMatrix tmp_1;
    Real ftmp1;
    Real ftmp2;
    int rho;
    int sigma;
    int ii;
    int jj;
    int kk;

    START_CODE();
  
    if (Nd > 4)
      QDP_error_exit("Hyp-smearing only implemented for Nd<=4",Nd);


    /*
     * Construct "level 1" smeared links in mu-direction with
     * staples only in one orthogonal direction, nu
     */
    ftmp1 = 1.0 - alpha3;
    ftmp2 = alpha3 / 2;
    ii = -1;
    for(int mu = 0; mu < Nd; ++mu)
    {
      for(int nu = 0; nu < Nd; ++nu)
      {
	if(nu == mu) continue;

	ii++;
	/*
	 * Forward staple
	 *
	 * u_tmp(x) = u(x,nu)*u(x+nu,mu)*u_dag(x+mu,nu)
	 */
	u_tmp = u[nu] * shift(u[mu],FORWARD,nu) * adj(shift(u[nu],FORWARD,mu));

	/*
	 * Backward staple
	 *
	 * u_tmp(x) += u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)
	 */
	u_tmp += shift(adj(u[nu]) * u[mu] * shift(u[nu],FORWARD,mu),BACKWARD,nu);

	/*
	 * Unprojected level 1 link
	 */
	tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
	u_tmp = adj(tmp_1);

	/*
	 * Project onto SU(Nc)
	 */
	if (Nd == 2)
	{
	  u_hyp[mu] = u[mu];
	  sun_proj(u_tmp, u_hyp[mu], BlkAccu, BlkMax);
	}
	else
	{
	  u_lv1[ii] = u[mu];
	  sun_proj(u_tmp, u_lv1[ii], BlkAccu, BlkMax);
	}
      }
    }

    if (Nd == 3)
    {
      /*
       * Construct hyp-smeared links in mu-direction with
       * "level 1" staples in the orthogonal direction, nu,
       * and the "level 1" links decorated in the 3-th orthogonal direction
       */
      ftmp1 = 1.0 - alpha2;
      ftmp2 = alpha2 / 4;
      for(int mu = 0; mu < Nd; ++mu)
      {
	u_tmp = 0;
	for(int nu = 0; nu < Nd; ++nu)
	{
	  if(nu == mu) continue;

	  /* 3-th orthogonal direction: rho */
	  for(jj = 0; jj < Nd; ++jj)
	  {
	    if(jj != mu && jj != nu) rho = jj;
	  }
	  jj = (Nd-1)*mu + rho;
	  if(rho > mu ) jj--;
	  kk = (Nd-1)*nu + rho;
	  if(rho > nu ) kk--;

	  /*
	   * Forward staple
	   *
	   * u_tmp(x) += u_lv1(x,kk)*u_lv1(x+nu,jj)*u_lv1_dag(x+mu,kk)
	   */
	  u_tmp += u_lv1[kk] * shift(u_lv1[jj],FORWARD,nu) * adj(shift(u_lv1[kk],FORWARD,mu));

	  /*
	   * Backward staple
	   *
	   * u_tmp(x) += u_lv1_dag(x-nu,kk)*u_lv1(x-nu,jj)*u_lv1(x-nu+mu,kk)
	   */
	  u_tmp += shift(adj(u_lv1[kk]) * u_lv1[jj] * shift(u_lv1[kk],FORWARD,mu),BACKWARD,nu);
	}

	/*
	 * Unprojected hyp-smeared link
	 */
	tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
	u_tmp = adj(tmp_1);

	/*
	 * Project onto SU(Nc)
	 */
	u_hyp[mu] = u[mu];
	sun_proj(u_tmp, u_hyp[mu], BlkAccu, BlkMax);
      }
    }
    else if (Nd == 4)
    {
      /*
       * Construct "level 2" smeared links in mu-direction with
       * "level 1" staples not in the orthogonal direction, nu,
       * and the "level 1" links decorated in the 4-th orthogonal direction
       */
      ftmp1 = 1.0 - alpha2;
      ftmp2 = alpha2 / 4;
      ii = -1;
      for(int mu = 0; mu < Nd; ++mu)
      {
	for(int nu = 0; nu < Nd; ++nu)
	{
	  if(nu == mu) continue;

	  ii++;
	  u_tmp = 0;
	  for(rho = 0; rho < Nd; ++rho)
	  {
	    if(rho == mu || rho == nu) continue;

	    /* 4-th orthogonal direction: sigma */
	    for(jj = 0; jj < Nd; ++jj)
	    {
	      if(jj != mu && jj != nu && jj != rho) sigma = jj;
	    }
	    jj = (Nd-1)*mu + sigma;
	    if(sigma > mu ) jj--;
	    kk = (Nd-1)*rho + sigma;
	    if(sigma > rho ) kk--;

	    /*
	     * Forward staple
	     *
	     * u_tmp(x) += u_lv1(x,kk)*u_lv1(x+rho,jj)*u_lv1_dag(x+mu,kk)
	     */
	    u_tmp += u_lv1[kk] * shift(u_lv1[jj],FORWARD,rho) * adj(shift(u_lv1[kk],FORWARD,mu));

	    /*
	     * Backward staple
	     *
	     * u_tmp(x) += u_lv1_dag(x-rho,kk)*u_lv1(x-rho,jj)*u_lv1(x-rho+mu,kk)
	     */
	    u_tmp += shift(adj(u_lv1[kk]) * u_lv1[jj] * shift(u_lv1[kk],FORWARD,mu),BACKWARD,rho);
	  }

	  /*
	   * Unprojected level 2 link
	   */
	  tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
	  u_tmp = adj(tmp_1);

	  /*
	   * Project onto SU(Nc)
	   */
	  u_lv2[ii] = u[mu];
	  sun_proj(u_tmp, u_lv2[ii], BlkAccu, BlkMax);
	}
      }

      /*
       * Construct hyp-smeared links in mu-direction with
       * "level 2" staples in the orthogonal direction, nu,
       * and the "level 2" links not decorated in the mu and nu directions
       */
      ftmp1 = 1.0 - alpha1;
      ftmp2 = alpha1 / 6;
      for(int mu = 0; mu < Nd; ++mu)
      {
	u_tmp = 0;
	for(int nu = 0; nu < Nd; ++nu)
	{
	  if(nu == mu) continue;

	  jj = (Nd-1)*mu + nu;
	  if(nu > mu ) jj--;
	  kk = (Nd-1)*nu + mu;
	  if(mu > nu ) kk--;

	  /*
	   * Forward staple
	   *
	   * u_tmp(x) += u_lv2(x,kk)*u_lv2(x+nu,jj)*u_lv2_dag(x+mu,kk)
	   */
	  u_tmp += u_lv2[kk] * shift(u_lv2[jj],FORWARD,nu) * adj(shift(u_lv2[kk],FORWARD,mu));

	  /*
	   * Backward staple
	   *
	   * u_tmp(x) += u_lv2_dag(x-nu,kk)*u_lv2(x-nu,jj)*u_lv2(x-nu+mu,kk)
	   */
	  u_tmp += shift(adj(u_lv2[kk]) * u_lv2[jj] * shift(u_lv2[kk],FORWARD,mu),BACKWARD,nu);
	}

	/*
	 * Unprojected hyp-smeared link
	 */
	tmp_1 = ftmp1*u[mu] + ftmp2*u_tmp;
	u_tmp = adj(tmp_1);

	/*
	 * Project onto SU(Nc)
	 */
	u_hyp[mu] = u[mu];
	sun_proj(u_tmp, u_hyp[mu], BlkAccu, BlkMax);
      }
    }

    END_CODE();
  }

} // ENd namespace Chroma
