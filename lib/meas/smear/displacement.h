
#ifndef __displacement_h__
#define __displacement_h__

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


void displacement(const multi1d<LatticeColorMatrix>& u, 
		  LatticeColorVector& chi, 
		  int length, int dir);

void displacement(const multi1d<LatticeColorMatrix>& u, 
		  LatticePropagator& chi, 
		  int length, int dir);

#endif
