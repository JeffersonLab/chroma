// $Id: unit_check.cc,v 1.2 2003-10-02 03:32:31 edwards Exp $

/*! \file
 *  \brief Test a gauge field is unitarized
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/unit_check.h"

using namespace QDP;


//! Check the unitarity of color matrix in SU(N)
/*!
 * \ingroup gauge
 *
 * \param  a  The LatticeColorMatrix to be tested
 */
void unitarityCheck(const multi1d<LatticeColorMatrix>& u)
{
  LatticeBoolean lbad;
  int numbad;

  for (int mu=0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp = u[mu];
    reunit(u_tmp, lbad, numbad, REUNITARIZE_ERROR);
  }
}


