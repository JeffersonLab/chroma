// $Id: constgauge.cc,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Constant diagonal gauge field
 */

#include "chromabase.h"
#include "util/gauge/taproj.h"
#include "util/gauge/expm12.h"
#include "util/gauge/reunit.h"

namespace Chroma
{

  //! Const diagonal gauge field
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   *  \param theta      Angles                        (Read)
   */

  void constgauge(multi1d<LatticeColorMatrix>& u, const multi2d<Real> theta)
  {
    START_CODE();

    u.resize(Nd);

    for(int mu = 0; mu < Nd; mu++)   // Loop over the directions
    {
      multi1d<Complex> phase(Nc); // the phases
      Complex tmp_u;

      tmp_u =
	cmplx(Real(0.0),Real(0.0)); // Set the matrix to zero

      u[mu] = tmp_u;

      phase[0] = cmplx(cos(theta(0, mu)), sin(theta(0,mu)));
      phase[1] = cmplx(cos(theta(1, mu)), sin(theta(1,mu)));
      phase[2] = cmplx(cos(theta(0, mu) + theta(1, mu)), 
		       -sin(theta(0,mu) + theta(1,mu)));

      // We now insert these elements into the appropriate places
      for(int i = 0; i < Nc; i++)
	pokeColor(u[mu], phase[i], i, i);

    }

    END_CODE();
  }
}
	




      
