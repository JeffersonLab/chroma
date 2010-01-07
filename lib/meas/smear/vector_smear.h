// -*- C++ -*-
// $Id: vector_smear.h,v 3.2 2008-11-04 17:26:16 edwards Exp $
/*! \file
 *  \brief vector smearing of color vector
 */

#ifndef __vector_smear_h__
#define __vector_smear_h__

#include "util/ferm/subset_vectors.h"

namespace Chroma 
{

  //! Do a covariant vector smearing of a lattice color vector field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param chi      color vector field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    parameter for the exponential smearing weight ( Read ) 
   *  \param j_decay  direction of decay ( Read )
   */
  void vectorSmear(LatticeColorVector& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay);


  //! Do a vector smearing of a color matrix
  /*!
   * Arguments:
   *
   *  \param chi      color matrix field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticeColorMatrix& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay);


  //! Do a t vector smearing of a lattice Fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    parameter for the exponential smearing weight ( Read ) 
   *  \param j_decay  direction of decay ( Read )
   */
  void vectorSmear(LatticeFermion& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay);

  //! Do a t vector smearing of a lattice Fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    parameter for the exponential smearing weight ( Read ) 
   *  \param j_decay  direction of decay ( Read )
   */
  void vectorSmear(LatticeStaggeredPropagator& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay);


  //! Do a t vector smearing of a lattice Fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    parameter for the exponential smearing weight ( Read ) 
   *  \param j_decay  direction of decay ( Read )
   */
  void vectorSmear(LatticePropagator& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay);

}  // end namespace Chroma

#endif
