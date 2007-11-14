//  $Id: stout_smear.cc,v 3.4 2007-11-14 02:06:15 bjoo Exp $
/*! \file
 *  \brief Stout smear a gauge field
 */

#include "chromabase.h"
#include "meas/smear/stout_smear.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

namespace Chroma 
{ 
  //! Stout-link smearing of the gauge configuration
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   * \param u_smear      smeared gauge field ( Write )
   * \param u            gauge field ( Read )
   * \param mu           direction of smeared gauge field ( Read )
   * \param sm_fact      smearing factor ( Read )
   * \param smear_dirs	 no staples in these directions (Read)
   */

  void stout_smear(LatticeColorMatrix& u_smear,
		   const multi1d<LatticeColorMatrix>& u,
		   int mu, 
		   const Real& sm_fact, 
		   multi1d<bool> smear_dirs)
  {
    START_CODE();
  
    // Construct and add the staples, only if direction smear_dirs is allowed
    LatticeColorMatrix u_staple = 0;

    for(int nu = 0; nu < Nd; ++nu)
    {
      if( nu != mu && smear_dirs[nu] )
      {
	// Forward staple
	u_staple += u[nu] * shift(u[mu], FORWARD, nu) * adj(shift(u[nu], FORWARD, mu));

	// Backward staple
	// tmp_1(x) = u_dag(x,nu)*u(x,mu)*u(x+mu,nu)
	LatticeColorMatrix tmp_1 = adj(u[nu]) * u[mu] * shift(u[nu], FORWARD, mu);

	// u_staple(x) += tmp_1_dag(x-nu)
	//             += u_dag(x+mu-nu,nu)*u_dag(x-nu,mu)*u(x-nu,nu)
	u_staple += shift(tmp_1, BACKWARD, nu);
      }
    }
    
    // The proto smeared link
    Real msmear_fact=Real(-1)*sm_fact;

    // - sign goes in here, because taproj() works with 
    // a minus sign relative to the stout paper.
    LatticeColorMatrix u_tmp = msmear_fact * u_staple * adj(u[mu]);

    // Take the trace-less anti-hermitian projection of the staple
    taproj(u_tmp);

    // Make it Hermitian, traceless
    // u_tmp = Q in Morningstar/Peardon's paper (hep-lat/0311018)
    // u_tmp = timesMinusI(u_tmp);
    // BJ: Changed convention.
    // eesu3 now takes iQ and multiplies by -i Internally so I can use it
    // for HMC too


    // Exactly exponentiate the Lie Algebra matrix
    // Now u_tmp = exp(iQ)
    expmat(u_tmp,EXP_EXACT);

    // Undo the back link
    u_smear = u_tmp * u[mu];

    END_CODE();
  }


} // End namespace Chroma
