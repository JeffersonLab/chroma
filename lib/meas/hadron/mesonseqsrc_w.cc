// $Id: mesonseqsrc_w.cc,v 1.6 2005-03-07 02:55:20 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesonseqsrc_w.h"

namespace Chroma 
{

  //! Construct pion sequential source
  /*!
   * \ingroup hadron
   *
   *  delta(tz-tx) exp(i p.z) \gamma_5 G \gamma_5
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return the sequential source before projection onto the sink
   */

  LatticePropagator mesPionSeqSrc(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3)
  {
    START_CODE();

    LatticePropagator src_prop_tmp;

    /*
     * Note, the seq. source is simply   src_prop_tmp = gamma_5 * U * gamma_5
     * where in this encoding gamma_5 = Gamma(15) .
     * However, to maintain compatibility with the calling  hadseqsrc_w.cc
     * code, the seqsource defined here is just   adj(U) .
     * That's because the hadseqsrc routine (for compatibility with baryons)
     * will construct finally   gamma_5 * adj(src_prop_tmp) * gamma_5 .
     * So, only the adj() is applied here to compensate for the adj()
     * that is subsequently done. The gamma_5 come along for free.
     */
    
//    int G5 = Ns*Ns-1;
//    src_prop_tmp = Gamma(G5) * quark_propagator_1 * Gamma(G5);

    src_prop_tmp = adj(quark_propagator_1);

    END_CODE();

    return src_prop_tmp;
  }

}  // end namespace Chroma
