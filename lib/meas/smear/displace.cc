//  $Id: displace.cc,v 3.8 2009-09-17 14:48:21 colin Exp $
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
	     int length, int dir,
	     const Subset& sub)
  {
    if (dir < 0 || dir >= Nd)
    {
      QDPIO::cerr << __func__ << ": invalid direction: dir=" << dir << endl;
      QDP_abort(1);
    }

    T tmp;
    T chi;
    chi[sub] = psi;

    if (length > 0)
    {
      for(int n = 0; n < length; ++n)
      {
	tmp[sub] = shift(chi, FORWARD, dir);
	chi[sub] = u[dir] * tmp;
      }
    }
    else // If length = or < 0.  If length == 0, does nothing.
    {
      for(int n = 0; n > length; --n)
      {
	tmp[sub] = shift(adj(u[dir])*chi, BACKWARD, dir);
	chi[sub] = tmp;
      }
    }
    return chi;
  }


  // Apply a displacement operator to a lattice field
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, int dir)
  {
    return displace<LatticeColorVector>(u, chi, length, dir, QDP::all);
  }


  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, int dir, const Subset& sub)
  {
    return displace<LatticeColorVector>(u, chi, length, dir, sub);
  }

  // Apply a displacement operator to a lattice field
  LatticePropagator displace(const multi1d<LatticeColorMatrix>& u, 
			     const LatticePropagator& chi, 
			     int length, int dir)
  {
    return displace<LatticePropagator>(u, chi, length, dir, QDP::all);
  }


  // Apply a displacement operator to a lattice field
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi, 
			  int length, int dir)
  {
    return displace<LatticeFermion>(u, chi, length, dir, QDP::all);
  }


  // Apply a displacement operator to a lattice field
  LatticeStaggeredFermion displace(const multi1d<LatticeColorMatrix>& u, 
				   const LatticeStaggeredFermion& chi, 
				   int length, int dir)
  {
    return displace<LatticeStaggeredFermion>(u, chi, length, dir, QDP::all);
  }


  // Apply a displacement operator to a lattice field
  LatticeStaggeredPropagator displace(const multi1d<LatticeColorMatrix>& u, 
				      const LatticeStaggeredPropagator& chi, 
				      int length, int dir)
  {
    return displace<LatticeStaggeredPropagator>(u, chi, length, dir, QDP::all);
  }


  // Apply a displacement operator to a lattice field
  LatticeColorMatrix displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorMatrix& chi, 
			      int length, int dir)
  {
    return displace<LatticeColorMatrix>(u, chi, length, dir, QDP::all);
  }



  //! Apply a displacement path to a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Read )
   *  \param length   displacement length - must be greater than zero ( Read )
   *  \param path     array of direction of displacement paths - pos/neg, or zero ( Read )
   *  \param sub      Subset of sites to act ( Read )
   *
   *  \return  displaced field
   */
  template<typename T>
  T displace(const multi1d<LatticeColorMatrix>& u, 
	     const T& psi, 
	     int displacement_length, 
	     const multi1d<int>& path,
	     const Subset& sub)
  {
    if (displacement_length < 0)
    {
      QDPIO::cerr << __func__ << ": invalid length=" << displacement_length << endl;
      QDP_abort(1);
    }

    T chi;
    chi[sub] = psi;

    for(int i=0; i < path.size(); ++i)
    {
      if (path[i] > 0)
      {
	int disp_dir = path[i] - 1;
	int disp_len = displacement_length;
	chi[sub] = displace<T>(u, chi, disp_len, disp_dir, sub);
      }
      else if (path[i] < 0)
      {
	int disp_dir = -path[i] - 1;
	int disp_len = -displacement_length;
	chi[sub] = displace<T>(u, chi, disp_len, disp_dir, sub);
      }
    }

    return chi;
  }


  // Apply a displacement path to a lattice field
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, const multi1d<int>& path)
  {
    return displace<LatticeColorVector>(u, chi, length, path, QDP::all);
  }

  // Apply a displacement path to a lattice field
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, const multi1d<int>& path,
			      const Subset& sub)
  {
    return displace<LatticeColorVector>(u, chi, length, path, sub);
  }


  // Apply a displacement path to a lattice field
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi, 
			  int length, const multi1d<int>& path)
  {
    return displace<LatticeFermion>(u, chi, length, path, QDP::all);
  }

  // Apply a displacement path to a lattice field
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi, 
			  int length, const multi1d<int>& path,
			  const Subset& sub)
  {
    return displace<LatticeFermion>(u, chi, length, path, sub);
  }


  //! Apply first deriv to the right onto source
  /*!
   * \ingroup sources
   *
   * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   *
   * \return $\f \nabla_\mu F(x)\f$
   */
  LatticeColorVector rightNabla(const LatticeColorVector& F, 
				const multi1d<LatticeColorMatrix>& u,
				int mu, int length)
  {
    return displace(u, F, length, mu) - displace(u, F, -length, mu);
  }

  //! Apply first deriv to the right onto source
  LatticeFermion rightNabla(const LatticeFermion& F, 
			    const multi1d<LatticeColorMatrix>& u,
			    int mu, int length)
  {
    return displace(u, F, length, mu) - displace(u, F, -length, mu);
  }

  //! Apply first deriv to the right onto source
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


  //! Apply "E_i" operator to the right onto source
  /*!
   * \ingroup smear
   *
   * \f$E_0 = (1/sqrt{2})*\nabla_x\nabla_x - \nabla_y\nabla_y\f$
   * \f$E_1 = -(1/sqrt{6})*\nabla_x\nabla_x + \nabla_y\nabla_y - 2*\nabla_z\nabla_z\f$
   *
   * \return $\f E_\alpha F(z,0) \f$
   */
  LatticePropagator rightE(const LatticePropagator& F,
			   const multi1d<LatticeColorMatrix>& u,
			   int mu, int length)
  {
    LatticePropagator tmp;

    switch (mu)
    {
    case 0:
      tmp  = rightNabla(rightNabla(F,u,0,length), u,0,length);
      tmp -= rightNabla(rightNabla(F,u,1,length), u,1,length);
      tmp *= Real(1)/Real(sqrt(Real(2)));
      break;

    case 1:
      tmp  = rightNabla(rightNabla(F,u,0,length), u,0,length);
      tmp += rightNabla(rightNabla(F,u,1,length), u,1,length);
      tmp -= Real(2)*rightNabla(rightNabla(F,u,2,length), u,2,length);
      tmp *= Real(-1)/Real(sqrt(Real(6)));
      break;

    default:
      QDPIO::cerr << __func__ << ": invalid direction for E: mu=" << mu << endl;
      QDP_abort(1);
    }

    return tmp;
  }


  //! Apply "Laplacian" operator to the right onto source
  /*!
   * \ingroup smear
   *
   * \f$Laplacian = \sum_{i=1}^3\nabla_i\nabla_i\f$
   *
   * \return $\f \nabla^2 F(z,0) \f$
   */
  LatticePropagator rightLap(const LatticePropagator& F,
			     const multi1d<LatticeColorMatrix>& u,
			     int length)
  {
    LatticePropagator tmp = zero;

    for(int i=0; i < 3; ++i)
    {
      tmp += rightNabla(rightNabla(F,u,i,length), u,i,length);
    }

    return tmp;
  }


}  // end namespace Chroma
