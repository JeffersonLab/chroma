// $Id: unit_check.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $

/*! \file
 *  \brief Test a gauge field is unitarized
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/unit_check.h"


namespace Chroma {

//! Check the unitarity of color matrix in SU(N)
/*!
 * \ingroup gauge
 *
 * \param  u  The LatticeColorMatrix to be tested
 */
void unitarityCheck(const multi1d<LatticeColorMatrix>& u)
{
  int numbad;

  for (int mu=0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp = u[mu];
    reunit(u_tmp, numbad, REUNITARIZE_ERROR);
  }
}

};

