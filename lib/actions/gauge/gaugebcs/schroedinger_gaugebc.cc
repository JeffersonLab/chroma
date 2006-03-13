// $Id: schroedinger_gaugebc.cc,v 2.1 2006-03-13 05:19:01 edwards Exp $
/*! \file
 *  \brief Schroedinger functional base class
 */

#include "gaugebc.h"

#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"

namespace Chroma 
{

  // Modify U fields in place
  void SchrGaugeBC::modify(multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    for(int mu=0; mu < u.size(); ++mu)
      copymask(u[mu], lSFmask()[mu], SFBndFld()[mu]);

    END_CODE();
  }

  // Zero some gauge-like field in place on the masked links
  void SchrGaugeBC::zero(multi1d<LatticeColorMatrix>& p) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for(int mu=0; mu < p.size(); ++mu)
      copymask(p[mu], lSFmask()[mu], z);

    END_CODE();
  }

}
