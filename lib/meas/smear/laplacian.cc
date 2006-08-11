// $Id: laplacian.cc,v 3.1 2006-08-11 16:13:30 edwards Exp $
/*! \file
 *  \brief Laplacian smearing of a source
 */

#include "chromabase.h"
#include "meas/smear/laplacian.h"
#include "actions/boson/operator/klein_gord.h"

namespace Chroma 
{

  //! Do a covariant Gaussian smearing of a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param j_decay  direction of decay ( Read )
   *  \param power    number of times to apply laplacian ( Read )
   */

  template<typename T>
  void laplacian(const multi1d<LatticeColorMatrix>& u, 
		 T& chi, 
		 int j_decay,
		 int power)
  {
    T psi;

    for (int p=0; p<power; p++) {
      psi = -1 * chi;

      /* hit with laplacian (Klein-Gordon with m=0) */
      klein_gord(u, psi, chi, 0, j_decay);
    }
  }

  void laplacian(const multi1d<LatticeColorMatrix>& u, 
		 LatticeColorVector& chi, 
		 int j_decay,
		 int power)
  {
    laplacian<LatticeColorVector>(u, chi, j_decay, power);
  }

  void laplacian(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& chi, 
		 int j_decay,
		 int power)
  {
    laplacian<LatticePropagator>(u, chi, j_decay, power);
  }

}  // end namespace Chroma

