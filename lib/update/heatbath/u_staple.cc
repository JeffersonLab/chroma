// $Id: u_staple.cc,v 2.0 2005-09-25 21:04:40 edwards Exp $
/*! \file
 *  \brief calculate the sum of all staples for a particular u_mu link
 */

#include "chromabase.h"
#include "update/heatbath/hb_params.h"
#include "update/heatbath/u_staple.h"

namespace Chroma 
{

  //! Compute staple from Wilson gauge action
  /*!
   * \ingroup heatbath
   *
   * \param u                    link field
   * \param mu                   direction of the link
   * \param u_mu_staple          staple attached to the link u
   * \param sub                  subset for updating
   * \param HBParams             container of HB parameters
   */
  void u_staple(LatticeColorMatrix& u_mu_staple, 
		const multi1d<LatticeColorMatrix>& u,
		int mu,
		const OrderedSubset& sub,
		const HBParams& hbp)
  {
    START_CODE();

    int t_dir=Nd-1; //for now it is here
    bool AnisoP=hbp.aniso();
    Double xi02=hbp.xi2();

    u_mu_staple = zero;
    for(int nu=0; nu<Nd; nu++) 
    {
      if( nu == mu ) continue;
      if( AnisoP && (mu == t_dir || nu == t_dir) )
      {
	// anisotropic lattice, time-like staple
	// +forward staple
	u_mu_staple[sub]+=
	  shift(u[nu],FORWARD,mu)*
	  shift(adj(u[mu]),FORWARD,nu)*
	  adj(u[nu])*
	  xi02;
	// +backward staple
	u_mu_staple[sub]+=
	  shift(shift(adj(u[nu]),FORWARD,mu),BACKWARD,nu)*
	  shift(adj(u[mu]),BACKWARD,nu)*
	  shift(u[nu],BACKWARD,nu)*
	  xi02;
      } else {
	// +forward staple
	u_mu_staple[sub]+=
	  shift(u[nu],FORWARD,mu)*
	  shift(adj(u[mu]),FORWARD,nu)*adj(u[nu]); 
	// +backward staple
	u_mu_staple[sub]+=
	  shift(shift(adj(u[nu]),FORWARD,mu),BACKWARD,nu)*
	  shift(adj(u[mu]),BACKWARD,nu)*
	  shift(u[nu],BACKWARD,nu);
      }
    }
    END_CODE();
  }

}  // end namespace Chroma
