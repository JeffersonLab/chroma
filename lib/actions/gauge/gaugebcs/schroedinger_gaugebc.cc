// $Id: schroedinger_gaugebc.cc,v 3.0 2006-04-03 04:58:54 edwards Exp $
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
  void SchrGaugeBC::zero(multi1d<LatticeColorMatrix>& ds_u) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for(int mu=0; mu < ds_u.size(); ++mu)
      copymask(ds_u[mu], lSFmask()[mu], z);

    END_CODE();
  }

}
