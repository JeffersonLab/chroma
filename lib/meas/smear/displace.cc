//  $Id: displace.cc,v 3.1 2006-12-02 18:16:28 edwards Exp $
/*! \file
 *  \brief Parallel transport a lattice field
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

#include "meas/smear/displace.h"
#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"

namespace Chroma 
{

  //! Apply a displacement operator to a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param psi      lattice field ( Read )
   *  \param length   length of displacement ( Read )
   *  \param dir      direction of displacement ( Read )
   *
   *  \return  displaced field
   */
  template<typename T>
  T displace(const multi1d<LatticeColorMatrix>& u, 
	     const T& psi, 
	     int length, int dir)
  {
    if (dir < 0 || dir >= Nd)
    {
      QDPIO::cerr << __func__ << ": invalid direction: dir=" << dir << endl;
      QDP_abort(1);
    }

    T chi = psi;

    if (length > 0)
    {
      for(int n = 0; n < length; ++n)
      {
	T tmp = shift(chi, FORWARD, dir);
	chi = u[dir] * tmp;
      }
    }
    else // If length = or < 0.  If length == 0, does nothing.
    {
      for(int n = 0; n > length; --n)
      {
	T tmp = shift(adj(u[dir])*chi, BACKWARD, dir);
	chi = tmp;
      }
    }
    return chi;
  }


  // Apply a displacement operator to a lattice field
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, int dir)
  {
    return displace<LatticeColorVector>(u, chi, length, dir);
  }


  // Apply a displacement operator to a lattice field
  LatticePropagator displace(const multi1d<LatticeColorMatrix>& u, 
			     LatticePropagator& chi, 
			     int length, int dir)
  {
    return displace<LatticePropagator>(u, chi, length, dir);
  }


  // Apply a displacement operator to a lattice field
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi, 
			  int length, int dir)
  {
    return displace<LatticeFermion>(u, chi, length, dir);
  }


  // Apply a displacement operator to a lattice field
  LatticeStaggeredFermion displace(const multi1d<LatticeColorMatrix>& u, 
				   const LatticeStaggeredFermion& chi, 
				   int length, int dir)
  {
    return displace<LatticeStaggeredFermion>(u, chi, length, dir);
  }


  // Apply a displacement operator to a lattice field
  LatticeStaggeredPropagator displace(const multi1d<LatticeColorMatrix>& u, 
				      LatticeStaggeredPropagator& chi, 
				      int length, int dir)
  {
    displace<LatticeStaggeredPropagator>(u, chi, length, dir);
  }



  //! Apply first deriv to the right onto source
  /*!
   * \ingroup sources
   *
   * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   *
   * \return $\f \nabla_\mu F(x)\f$
   */
  LatticePropagator rightNabla(const LatticePropagator& F, 
			       const multi1d<LatticeColorMatrix>& u,
			       int mu, int length)
  {
    return displace(u, F, length, mu) - displace(u, F, -length, mu);
  }


  //! Apply "D_i" operator to the right onto source
  /*!
   * \ingroup sources
   *
   * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
   *
   * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
   * 
   * \return $\f F(z,0) D_\mu\f$
   */
  LatticePropagator rightD(const LatticePropagator& F,
			   const multi1d<LatticeColorMatrix>& u,
			   int mu, int length)
  {
    LatticePropagator tmp = zero;

    // Slow implementation - to speed up could compute once the \nabla_j deriv
    for(int j=0; j < 3; ++j)
      for(int k=0; k < 3; ++k)
      {
	if (symTensor3d(mu,j,k) != 0)
	  tmp += rightNabla(rightNabla(F,u,j,length), u,k,length);
      }

    return tmp;
  }


  //! Apply "B_i" operator to the right onto source
  /*!
   * \ingroup sources
   *
   * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
   *
   * \return $\f F(z,0) B_\mu\f$
   */
  LatticePropagator rightB(const LatticePropagator& F,
			   const multi1d<LatticeColorMatrix>& u,
			   int mu, int length)
  {
    LatticePropagator tmp = zero;

    // Slow implementation - to speed up could compute once the \nabla_j deriv
    for(int j=0; j < 3; ++j)
      for(int k=0; k < 3; ++k)
      {
	if (antiSymTensor3d(mu,j,k) != 0)
	  tmp += Real(antiSymTensor3d(mu,j,k)) * rightNabla(rightNabla(F,u,j,length), u,k,length);
      }

    return tmp;
  }


}  // end namespace Chroma
