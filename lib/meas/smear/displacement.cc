

#include "chromabase.h"
#include "meas/smear/displacement.h"

using namespace QDP;

//! apply a displacement operator to a lattice field
/*!
 * Arguments:
 *
 *  \param u        gauge field ( Read )
 *  \param chi      color vector field ( Modify )
 *  \param length   length of displacement ( Read )
 *  \param dir      direction of displacement ( Read )
 *
 *
 * Discription:
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
 */

template<typename T>
void displacement(const multi1d<LatticeColorMatrix>& u, 
		  T& chi, 
		  int length, int dir)
{
  // Initial ferm field
  T psi = chi;

  if (length > 0)
    for(int n = 0; n < length; ++n)
      {
	T tmp = shift(psi, FORWARD, dir);
	psi = u[dir] * tmp;
      }

  else 
    for(int n = 0; n < -1*length; ++n)
      {
	T tmp = shift(psi, BACKWARD, dir);
	psi = u[dir] * tmp;
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


