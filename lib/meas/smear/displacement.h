// -*- C++ -*-
//  $Id: displacement.h,v 3.1 2008-09-26 19:53:28 edwards Exp $
/*! \file
 *  \brief Parallel transport a lattice field
 */

#ifndef __displacement_h__
#define __displacement_h__

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

  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeColorVector& chi, 
		    int length, int dir);


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

  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticePropagator& chi, 
		    int length, int dir);


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeFermion& chi, 
		    int length, int dir) ;


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeStaggeredFermion& chi, 
		    int length, int dir) ;


  void displacement(const multi1d<LatticeColorMatrix>& u, 
		    LatticeStaggeredPropagator& chi, 
		    int length, int dir) ;  

}  // end namespace Chroma

#endif
