// -*- C++ -*-
// $Id: displace.h,v 3.8 2009-09-17 14:48:21 colin Exp $
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

#ifndef __displace_h__
#define __displace_h__

#include "chromabase.h"

namespace Chroma 
{

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
   *
   *  \return  displaced field
   */
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, const multi1d<int>& path);


  //! Apply a displacement path to a lattice field
  /*! \ingroup smear */
  LatticeColorMatrix displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorMatrix& chi, 
			      int length, const multi1d<int>& path,
			      const Subset& sub);


  //! Apply a displacement path to a lattice field
  /*! \ingroup smear */
  LatticeColorMatrix displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorMatrix& chi, 
			      int length, const multi1d<int>& path);


  //! Apply a displacement path to a lattice field
  /*! \ingroup smear */
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, const multi1d<int>& path,
			      const Subset& sub);


  //! Apply a displacement path to a lattice field
  /*! \ingroup smear */
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi,
			  int length, const multi1d<int>& path);


  //! Apply a displacement path to a lattice field
  /*! \ingroup smear */
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi,
			  int length, const multi1d<int>& path,
			  const Subset& sub);


  //! Apply a displacement operator to a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Read )
   *  \param length   length of displacement - can be negative ( Read )
   *  \param dir      direction of displacement ( Read )
   *
   *  \return  displaced field
   */
  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, int dir);

  LatticeColorVector displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorVector& chi, 
			      int length, int dir,
			      const Subset& sub);

  //! Apply a displacement operator to a lattice field
  /*! \ingroup smear */
  LatticeFermion displace(const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& chi, 
			  int length, int dir);


  //! Apply a displacement operator to a lattice field
  /*! \ingroup smear */
  LatticeStaggeredFermion displace(const multi1d<LatticeColorMatrix>& u, 
				   const LatticeStaggeredFermion& chi, 
				   int length, int dir);


  //! Apply a displacement operator to a lattice field
  /*! \ingroup smear */
  LatticeStaggeredPropagator displace(const multi1d<LatticeColorMatrix>& u, 
				      const LatticeStaggeredPropagator& chi, 
				      int length, int dir);  

  //! Apply a displacement operator to a lattice field
  /*! \ingroup smear */
  LatticeColorMatrix displace(const multi1d<LatticeColorMatrix>& u, 
			      const LatticeColorMatrix& chi, 
			      int length, int dir);  

  //! Apply a displacement operator to a lattice field
  /*! \ingroup smear */
  LatticePropagator displacement(const multi1d<LatticeColorMatrix>& u, 
				 const LatticePropagator& chi, 
				 int length, int dir);
  

  //----------------------------------------------------------------------------------
  //! Apply a right nabla path to a lattice field
  /*!
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Read )
   *  \param length   displacement length - must be greater than zero ( Read )
   *  \param path     array of direction of derivative paths - pos or zero ( Read )
   *
   *  \return  deriv field
   */
  LatticeColorVector rightNabla(const multi1d<LatticeColorMatrix>& u, 
				const LatticeColorVector& chi, 
				int length, const multi1d<int>& path);


  //! Apply a right nabla path to a lattice field
  /*! \ingroup smear */
  LatticeColorMatrix rightNabla(const multi1d<LatticeColorMatrix>& u, 
				const LatticeColorMatrix& chi, 
				int length, const multi1d<int>& path);



  //! Apply first deriv to the right onto source
  /*!
   * \ingroup smear
   *
   * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   *
   * \return $\f \nabla_\mu F(x)\f$
   */
  LatticeColorVector rightNabla(const LatticeColorVector& F, 
				const multi1d<LatticeColorMatrix>& u,
				int mu, int length);

  //! Apply first deriv to the right onto source 
  LatticeFermion rightNabla(const LatticeFermion& F, 
			    const multi1d<LatticeColorMatrix>& u,
			    int mu, int length);

  //! Apply first deriv to the right onto source
  LatticePropagator rightNabla(const LatticePropagator& F, 
			       const multi1d<LatticeColorMatrix>& u,
			       int mu, int length);


  //----------------------------------------------------------------------------------
  //! Apply left-right deriv to the right onto source
  /*!
   * \ingroup sources
   *
   * \f$\right\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   * \f$\left\nabla_\mu f(x) = f(x+\mu)U_{-\mu}(x) - f(x-\mu)U_\mu(x-\mu)\f$
   *
   * \return $\f (\left\nabla_\mu - \right\nabla_\mu) F(x)\f$
   */
  LatticeColorVector leftRightNabla(const LatticeColorVector& F, 
				    const multi1d<LatticeColorMatrix>& u,
				    int mu, int length,
				    int mom);

  //! Apply left-right deriv to the right onto source
  LatticeFermion leftRightNabla(const LatticeFermion& F, 
				const multi1d<LatticeColorMatrix>& u,
				int mu, int length,
				int mom);

  //! Apply left-right deriv to the right onto source
  LatticePropagator leftRightNabla(const LatticePropagator& F, 
				   const multi1d<LatticeColorMatrix>& u,
				   int mu, int length,
				   int mom);


  //! Apply left-right deriv to the right onto source
  LatticeColorVector leftRightNabla(const multi1d<LatticeColorMatrix>& u, 
				    const LatticeColorVector& chi, 
				    int length, 
				    const multi1d<int>& path,
				    const multi1d<int>& mom);

  //! Apply left-right deriv to the right onto source
  LatticeColorMatrix leftRightNabla(const multi1d<LatticeColorMatrix>& u, 
				    const LatticeColorMatrix& chi, 
				    int length, 
				    const multi1d<int>& path,
				    const multi1d<int>& mom);


  //----------------------------------------------------------------------------------
  //! Apply "D_i" operator to the right onto source
  /*!
   * \ingroup smear
   *
   * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
   *
   * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
   * 
   * \return $\f D_\mu F(z,0)\f$
   */
  LatticePropagator rightD(const LatticePropagator& F,
			   const multi1d<LatticeColorMatrix>& u,
			   int mu, int length);


  //! Apply "B_i" operator to the right onto source
  /*!
   * \ingroup smear
   *
   * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
   *
   * \return $\f B_\mu F(z,0) \f$
   */
  LatticePropagator rightB(const LatticePropagator& F,
			   const multi1d<LatticeColorMatrix>& u,
			   int mu, int length);


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
			   int mu, int length);


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
			     int length);

}  // end namespace Chroma

#endif
