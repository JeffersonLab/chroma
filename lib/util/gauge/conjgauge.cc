// $Id: conjgauge.cc,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Take the complex conjugate of a gauge field
 */

#include "chromabase.h"

namespace Chroma
{

  //! Take complex conjugate of gauge field u
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   *  \param g          Gauge transformation matrices (Write)
   */

  void conjgauge(multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();
  
    
    for(int mu = 0; mu < Nd; ++mu)
    {
      LatticeColorMatrix u_tmp =  conj(u[mu]);
      u[mu] = u_tmp;
    }
    
    END_CODE();
  }
}
