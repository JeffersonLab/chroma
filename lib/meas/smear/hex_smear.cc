//  $Id: hyp_smear.cc,v 3.1 2006-08-11 16:13:30 edwards Exp $
/*! \file
 *  \brief Hyp-smearing of the gauge configuration
 */

#include "chromabase.h"
#include "meas/smear/hex_smear.h"
#include "util/gauge/eesu3.h"
#include "util/gauge/taproj.h"

namespace Chroma 
{ 
  //! Construct one iteration of the "hex-smeared" links of Capitani et al
  /*! hep-lat/0607006
   * \ingroup smear
   *
   * Construct the "hyp-smeared" links of Capitani et al., with
   * staple coefficients alpha1, alpha2 and alpha3
   *
   * Arguments:
   *
   *  \param u		gauge field (Read)
   *  \param u_hyp	"hex-smeared" gauge field (Write)
   *
   */

  void Hex_Smear_onestepp(multi1d<LatticeColorMatrix>& u,
		 multi1d<LatticeColorMatrix>& u_hyp)
  {
    multi1d<LatticeColorMatrix> u_lv1(Nd*(Nd-1));
    multi1d<LatticeColorMatrix> u_lv2(Nd*(Nd-1));
    LatticeColorMatrix u_tmp;
    LatticeColorMatrix tmp_1;
    LatticeComplex tt;

    /** see equation 4 in hep-lat/0607006 **/
    /** the parameters below were empirically obtained
        by BMW personnel **/

    const Real hex_alpha1 = 0.95 ;
    const Real hex_alpha2 = 0.76 ; 
    const Real hex_alpha3 = 0.38 ;

    /**  to get agreement with BMW code **/
    const Real hex_alpha3_fact = hex_alpha3 / 2.0  ;
    const Real hex_alpha2_fact = hex_alpha2 / 4.0  ;
    const Real hex_alpha1_fact = hex_alpha1 / 6.0  ;

    int rho;
    int sigma;
    int ii;
    int jj;
    int kk;

    START_CODE();
  
    if (Nd != 4)
      QDP_error_exit("Hex-smearing only implemented for Nd=4",Nd);

    /*
     * Construct "level 1" smeared links in mu-direction with
     * staples only in one orthogonal direction, nu
     */

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
	tmp_1 = u_tmp * adj(u[mu]) ;

	/*
	 * Exponentiation of AntiHermitian traceless matrix
	 */
	  u_tmp = hex_alpha3_fact * tmp_1 ;
	  taproj(u_tmp) ;
	  eesu3(u_tmp) ;

	  u_lv1[ii] = u_tmp * u[mu];

      } // end loop over nu
    }   // end loop over mu

      /*
       * Construct "level 2" smeared links in mu-direction with
       * "level 1" staples not in the orthogonal direction, nu,
       * and the "level 1" links decorated in the 4-th orthogonal direction
       */

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
	  tmp_1 = u_tmp * adj(u[mu]) ;

	  /*
           * Exponentiation of AntiHermitian traceless matrix
	   */
	  u_tmp     = hex_alpha2_fact * tmp_1 ;
	  taproj(u_tmp) ;
	  eesu3( u_tmp  ) ;
	  u_lv2[ii] = u_tmp * u[mu] ;

	}
      }

      /*
       * Construct hyp-smeared links in mu-direction with
       * "level 2" staples in the orthogonal direction, nu,
       * and the "level 2" links not decorated in the mu and nu directions
       */

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
	tmp_1 =  u_tmp * adj(u[mu]) ;

	/*
         * Exponentiation of AntiHermitian traceless matrix
         */
	 u_tmp = hex_alpha1_fact * tmp_1 ;
	 taproj(u_tmp) ;
	 eesu3( u_tmp ) ;

	 u_hyp[mu] = u_tmp * u[mu] ;
      }


    END_CODE();
  }

  //! Construct nstep iterations of the "hex-smeared" links of Capitani et al
  /*! hep-lat/0607006
   * \ingroup smear
   *
   * Construct the "hyp-smeared" links of Capitani et al., 
   *
   * Arguments:
   *
   *  \param u		gauge field (Read)
   *  \param u_hyp	"hex-smeared" gauge field (Write)
   *  \param nstep	number of hex steps (Read)
   *
   */



  void Hex_Smear(const multi1d<LatticeColorMatrix>& u,
                 multi1d<LatticeColorMatrix>& u_hyp, const int nstep) 
  {
    multi1d<LatticeColorMatrix> u_last(Nd);
    int i ; 
    START_CODE();

    if( nstep == 0)
      {
	u_hyp = u ;
      }
    else
      {
	u_last = u  ;
	for(i=0 ; i < nstep ; ++i)
	  {
	    Hex_Smear_onestepp(u_last,u_hyp) ;
	    u_last = u_hyp  ;
	  }
      }

    END_CODE();
  }

} // ENd namespace Chroma
