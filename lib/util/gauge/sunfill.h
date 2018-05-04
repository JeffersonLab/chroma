// -*- C++ -*-
/*! \file
 *  \brief  Fill an SU(Nc) matrix with an SU(2) submatrix
 */

#ifndef __sunfill_h__
#define __sunfill_h__

namespace Chroma 
{

  //! Fill a dest(su2_index) <- r_0,r_1,r_2,r_3  under a subset
  /*!
   * \ingroup gauge
   *
   * Fill an SU(Nc) matrix V with the SU(2) submatrix su2_index
   * paramtrized by b_k in the sigma matrix basis.
   *
   * Fill in B from B_SU(2) = b0 + i sum_k bk sigma_k
   *
   * There are Nc*(Nc-1)/2 unique SU(2) submatrices in an SU(Nc) matrix.
   * The user does not need to know exactly which one is which, just that
   * they are unique.
   *
   * Arguments:
   *
   * \param dest       su(n) matrix
   * \param r          su2 matrix represented in the O(4) rep. - an array of LatticeReal 
   * \param su2_index  int lying in [0, Nc*(Nc-1)/2)
   * \param s          subset for operations       (Read)
   */
  void
  sunFill(LatticeColorMatrix& dest,
	  const multi1d<LatticeReal>& r,
	  int su2_index,
	  const Subset& s);


}
  // end namespace Chroma

#endif
