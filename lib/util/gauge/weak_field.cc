// $Id: weak_field.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
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
    Real fact = 1.0 / sqrt(2.0);

    for(int mu = 0; mu < Nd; mu++)
    {
      u[mu] = 1.0;

      for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; su2_index++)
      {
	// Embed a SU(2) within SU(N)
	random(theta);
	theta *= amp;
	a[0] = fact*cos(theta);
	a[1] = fact*sin(theta);

	random(theta);
	theta *= amp;
	a[2] = fact*cos(theta);
	a[3] = fact*sin(theta);

	s1 = 1.0;
	sunFill(s1, a, su2_index, all);

	// Mix up the indices by multiplying - still unitary
	s2 = u[mu] * s1;
	u[mu] = s2;
      }
    }

    END_CODE();
  }
}
	




      
