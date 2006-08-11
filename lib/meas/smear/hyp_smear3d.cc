//  $Id: hyp_smear3d.cc,v 3.1 2006-08-11 16:13:30 edwards Exp $
/*! \file
 *  \brief Hyp-smearing of the gauge configuration
 */

#include "chromabase.h"
#include "meas/smear/hyp_smear3d.h"
#include "util/gauge/sun_proj.h"

namespace Chroma 
{ 
  //! Construct the "hyp-smeared" links of Anna Hasenfratz involving only the spatial links
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
   *  \param alpha1	staple coefficient "1" (Read, not used)
   *  \param alpha2	staple coefficient "2" (Read)
   *  \param alpha3	staple coefficient "3" (Read)
   *  \param BlkAccu	accuracy in SU(Nc) projection (Read)
   *  \param BlkMax	max number of iterations in SU(Nc) projection (Read)
   *  \param j_decay	direction of no staple(Read)
   */

  void Hyp_Smear3d(const multi1d<LatticeColorMatrix>& u,
		   multi1d<LatticeColorMatrix>& u_hyp,
		   const Real& alpha1, const Real& alpha2, const Real& alpha3,
		   const Real& BlkAccu, int BlkMax,int j_decay)
  {
    multi1d<LatticeColorMatrix> u_lv1((Nd-1)*(Nd-2));
    LatticeColorMatrix u_tmp;
    LatticeColorMatrix tmp_1;
    Real ftmp1;
    Real ftmp2;
    int rho;
    int ii;
    int jj;
    int kk;

    START_CODE();
  
    if (Nd > 4)
      QDP_error_exit("Hyp-smearing-3D only implemented for Nd=4",Nd);

    /*
     * Construct "level 1" smeared links in mu-direction
     * with staples only in one orthogonal direction, nu
     */
    QDPIO::cout << "HYP-smearing-3D only involving spatial links!" << endl;

    u_hyp = u;   // only need to make sure j_decay direction is set

    ftmp1 = 1.0 - alpha3;
    ftmp2 = alpha3 / 2;
    ii = -1;
    for(int mu = 0; mu < Nd; ++mu) if ( mu != j_decay )
    {
      for(int nu = 0; nu < Nd; ++nu) if (nu != mu && nu != j_decay )
      {
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
	u_lv1[ii] = u[mu];
	sun_proj(u_tmp, u_lv1[ii], BlkAccu, BlkMax);
      }
    }

    /*
     * Construct hyp-smeared links in mu-direction with
     * "level 1" staples in the orthogonal direction, nu,
     * & "level 1" links decorated in 3-rd orthogonal direction
     */

    ftmp1 = 1.0 - alpha2;
    ftmp2 = alpha2 / 4;
    for(int mu = 0; mu < Nd; ++mu) if ( mu != j_decay )
    {
      u_tmp = 0;
      for(int nu = 0; nu < Nd; ++nu) if (nu != mu && nu != j_decay )
      {
	/* 3-th orthogonal direction: rho */
	for(jj = 0; jj < Nd; ++jj) 
	{
	  if(jj != mu && jj != nu && jj != j_decay) rho = jj;
	}
	jj = (Nd-2)*mu + rho;
	if(rho > mu ) jj--;
	kk = (Nd-2)*nu + rho;
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

    END_CODE();
  }

} // ENd namespace Chroma
