
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
  void unitarityCheck(const multi1d<LatticeColorMatrixF3>& u)
  {
    START_CODE();

    int numbad;

    for (int mu=0; mu < Nd; ++mu)
    {
      LatticeColorMatrixF3 u_tmp = u[mu];
      reunit(u_tmp, numbad, REUNITARIZE_ERROR);
    }
    
    END_CODE();
  }

 void unitarityCheck(const multi1d<LatticeColorMatrixD3>& u)
  {
    START_CODE();

    int numbad;

    for (int mu=0; mu < Nd; ++mu)
    {
      LatticeColorMatrixD3 u_tmp = u[mu];
      reunit(u_tmp, numbad, REUNITARIZE_ERROR);
    }
    
    END_CODE();
  }

}

