// $Id: hotst.cc,v 1.1 2003-01-26 04:28:41 edwards Exp $
// HOTST

/*! \file
 *  \brief Set a gauge field from a sample of (almost) Haar measure
 */

/*!
 * u = exp(A), where A = random traceless antihermitian matrix.
 *
 * Arguments:
 *
 *  \param u          Gauge field               (Modify)
 *
 *
 * NOTE!!!!: we are using a hack. The twelth_order exponentiation followed
 *   by a reunitarization should be replace by an exact exponentiation.
 *   So, the gauge field distriubution is NOT from Haar measure, but
 *   if someone can think of a really good reason why it should be Haar
 *   measure, then something can be done about it. 
 */

#include "qdp.h"
#include "propagator.h"

using namespace QDP;

void HotSt(multi1d<LatticeColorMatrix>& u)
{
  START_CODE("HotSt");

  for(int mu=0; mu<u.size(); ++mu)
  {
    gaussian(u[mu]); // Gaussian fill
    taproj(u[mu]);   // Traceless anti-hermitian projection
    expm12(u[mu]);   // Exponentiate
    reunit(u[mu]);   // Reunitarize
  }

  END_CODE("HotSt");
}
