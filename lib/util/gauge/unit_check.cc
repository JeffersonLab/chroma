// $Id: unit_check.cc,v 1.5 2005-01-14 15:59:00 bjoo Exp $

/*! \file
 *  \brief Test a gauge field is unitarized
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/unit_check.h"

using namespace QDP;
using namespace Chroma; 

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

