// $Id: hotst.cc,v 1.3 2003-02-15 05:54:26 edwards Exp $
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

#include "szinbase.h"
#include "util/gauge/hotst.h"
#include "util/gauge/taproj.h"
#include "util/gauge/expm12.h"
#include "util/gauge/reunit.h"

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
