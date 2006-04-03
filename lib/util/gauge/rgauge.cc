// $Id: rgauge.cc,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Do a random gauge transformation on the u fields
 */

#include "chromabase.h"
#include "util/gauge/taproj.h"
#include "util/gauge/expm12.h"
#include "util/gauge/reunit.h"


namespace Chroma
{

  //! Do a random gauge transformation on the u fields
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   *  \param g          Gauge transformation matrices (Write)
   */

  void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g)
  {
    START_CODE();
  
    // g = exp(A), where A = random traceless antihermitian matrix.
  
    // NOTE!!!!: we are using a hack. The twelth_order exponentiation followed
    //   by a reunitarization should be replace by an exact exponentiation.
    //   So, the gauge field distriubution is NOT from Haar measure, but
    //   if someone can think of a really good reason why it should be Haar
    //   measure, then something can be done about it.
  
    gaussian(g);

    taproj(g);
    expm12(g);
    reunit(g);
    
    for(int mu = 0; mu < Nd; ++mu)
    {
      LatticeColorMatrix u_tmp = g * u[mu];
      u[mu] = u_tmp * shift(adj(g), FORWARD, mu);
    }
    
    END_CODE();
  }


  //! Do a random gauge transformation on the u fields
  /*!
   * \ingroup gauge
   *
   * Convenience function: does not return gauge transformation matrices
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   */

  void rgauge(multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();
  
    LatticeColorMatrix g;
    rgauge(u, g);

    END_CODE();
  }


}
