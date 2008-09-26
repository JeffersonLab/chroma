//  $Id: displacement.cc,v 3.1 2008-09-26 19:53:28 edwards Exp $
/*! \file
 *  \brief Parallel transport a lattice field
 */

#include "chromabase.h"
#include "meas/smear/displacement.h"

namespace Chroma 
{

  //! Apply a displacement operator to a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param length   length of displacement ( Read )
   *  \param dir      direction of displacement ( Read )
   *
   *
   * Description:
   *
   *  Suppose q(x) is a quark field.
   *  Displacement operator D_j^{(p)} moves quark field 
   *  for p lattice sites to the direction j in covariant
   *  fashion.
   *
   *  Namely, 
   *  D_j^{(p)} q(x) = U_j(x) U_j(x+j) U_j(x+2j)...U_j(x+(p-1)j) q(x+pj),
   *  where U is the gauge-link.
   *
   *  dir: x(0), y(1), z(2)
   *
   */

  template<typename T>
  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    T& chi, 
		    int length, int dir)
  {
    if (length > 0)
      for(int n = 0; n < length; ++n)
      {
	T tmp = shift(chi, FORWARD, dir);
	chi = u[dir] * tmp;
      }

    else // If length = or < 0.  If length == 0, does nothing.
      for(int n = 0; n > length; --n)
      {
	T tmp = shift(adj(u[dir])*chi, BACKWARD, dir);
	chi = tmp;
      }
  }


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeColorVector& chi, 
		    int length, int dir)
  {
    displacement<LatticeColorVector>(u, chi, length, dir);
  }


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticePropagator& chi, 
		    int length, int dir)
  {
    displacement<LatticePropagator>(u, chi, length, dir);
  }

  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeFermion& chi, 
		    int length, int dir)
  {
    displacement<LatticeFermion>(u, chi, length, dir);
  }


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeStaggeredFermion& chi, 
		    int length, int dir)
  {
    displacement<LatticeStaggeredFermion>(u, chi, length, dir);
  }



  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeStaggeredPropagator& chi, 
		    int length, int dir)
  {
    displacement<LatticeStaggeredPropagator>(u, chi, length, dir);
  }

}  // end namespace Chroma
