// $Id: weak_field.cc,v 2.1 2005-12-18 03:51:10 edwards Exp $
/*! \file
 *  \brief Construct a weak field
 */

#include "chromabase.h"
#include "util/gauge/weak_field.h"
#include "util/gauge/sunfill.h"

namespace Chroma
{

  //! Construct a weak field
  /*!
   * \ingroup gauge
   *
   * Make a gauge field close to free field - e.g., small fluctations
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   */
  void weakField(multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();

    u.resize(Nd);

    Real amp = 0.02 * twopi;
    LatticeReal theta;
    LatticeColorMatrix s1, s2;
    multi1d<LatticeReal> a(4);

    for(int mu = 0; mu < Nd; mu++)
    {
      // Embed a SU(2) within SU(N)
      random(theta);
      theta *= amp;

      a[0] = cos(theta);
      a[1] = sin(theta);
      a[2] = zero;
      a[3] = zero;

      s1 = 1.0;
      sunFill(s1, a, 0, all);

      // Embed another SU(2) within SU(N)
      a[0] = cos(theta);
      a[1] = zero;
      a[2] = sin(theta);
      a[3] = zero;

      s2 = 1.0;
      sunFill(s2, a, 1, all);

      // Mix up the indices by multiplying - still unitary
      u[mu] = s1 * s2;
    }

    END_CODE();
  }
}
	




      
