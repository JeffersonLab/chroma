//  $Id: stout_smear.cc,v 1.1 2004-02-03 22:39:06 edwards Exp $
/*! \file
 *  \brief Stout-link smearing of the gauge configuration
 */

#error "Not finished"

#include "chromabase.h"
#include "meas/smear/stout_smear.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

using namespace QDP;

//! Stout-link smearing of the gauge configuration
/*!
 * \ingroup smear
 *
 * Arguments:
 *
 * \param u_smear  smeared gauge field ( Write )
 * \param u        gauge field ( Read )
 * \param mu       direction of smeared gauge field ( Read )
 * \param sm_fact  smearing factor ( Read )
 * \param j_decay  no staple in direction j_decay ( Read )
 */

void stout_smear(LatticeColorMatrix& u_smear,
		 const multi1d<LatticeColorMatrix>& u,
		 int mu, 
		 const Real& sm_fact, int j_decay)
{
  START_CODE("APE_Smear");
  
  // Construct and add the staples, except in direction j_decay
  LatticeColorMatrix u_staple = 0;

  for(int nu = 0; nu < Nd; ++nu)
  {
    if( nu != mu && nu != j_decay )
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
  LatticeColorMatrix u_tmp = sm_fact * u_staple * adj(u[mu]);

  // Take the trace-less anti-hermitian projection of the staple
  taproj(u_tmp);

  // Exactly exponentiate the Lie Algebra matrix
  expmat(u_tmp,EXP_EXACT);

  // Undo the back link
  u_smear = u_tmp * u[mu];

  END_CODE("stout_smear");
}

