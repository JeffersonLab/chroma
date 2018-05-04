/*! \file
 *  \brief Gamma5 hermiticity
 */

#include "util/ferm/gamma5_herm_w.h"

namespace Chroma
{
  // Return gamma_5*adj(source)*gamma_f
  LatticePropagator gamma5Herm(const LatticePropagator& source_prop)
  {
    int G5 = Ns*Ns-1;
    return Gamma(G5) * adj(source_prop) * Gamma(G5);
  }

}
