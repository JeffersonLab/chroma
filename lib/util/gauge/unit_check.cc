
/*! \file
 *  \brief Test a gauge field is unitarized
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/unit_check.h"


namespace Chroma 
{

  //! Check the unitarity of color matrix in SU(N)
  /*!
   * \ingroup gauge
   *
   * \param  u  The LatticeColorMatrix to be tested
   */
  template<typename Q>
  void unitarityCheck(const multi1d<Q>& u)
  {
    START_CODE();

    int numbad;

    for (int mu=0; mu < Nd; ++mu)
    {
      Q u_tmp = u[mu];
      reunit(u_tmp, numbad, REUNITARIZE_ERROR);
    }
    
    END_CODE();
  }

  // Instantiate template for various known types.
#if QDP_NC!=3
  template void unitarityCheck(const multi1d<LatticeColorMatrix>& u);
#endif
  template void unitarityCheck(const multi1d<LatticeColorMatrixF3>& u);
  template void unitarityCheck(const multi1d<LatticeColorMatrixD3>& u);  
}

