// -*- C++ -*-
// $Id: su2extract.h,v 1.3 2003-12-06 17:35:57 edwards Exp $
/*! \file
 *  \brief  Extract an unnormalized SU(2) matrix from a GL(3,C) matrix
 */

#ifndef __su2extract__
#define __su2extract__

//! Su2_extract: r_0,r_1,r_2,r_3 <- source(su2_index)  [SU(N) field]  under a subset
/*! 
 * \ingroup gauge
 *
 * Extract components r_k proportional to SU(2) submatrix su2_index
 * from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the
 * sigma matrix basis.
 *
 * There are Nc*(Nc-1)/2 unique SU(2) submatrices in an SU(Nc) matrix.
 * The user does not need to know exactly which one is which, just that
 * they are unique.
 *
 * Arguments:
 *
 * \param source      extract the su2 matrix from this su(n) gauge field
 * \param su2_index   int lying in [0, Nc*(Nc-1)/2)
 * \param s           subset for operations       (Read)
 *
 * \return su2 matrix represented in the O(4) rep. - an array of LatticeReal 
 */

multi1d<LatticeReal> 
su2Extract(const LatticeColorMatrix& source, 
	   int su2_index, 
	   const UnorderedSubset& s);


//! Su2_extract: r_0,r_1,r_2,r_3 <- source(su2_index)  [SU(N) field]  under a subset
/*! 
 * \ingroup gauge
 *
 * Extract components r_k proportional to SU(2) submatrix su2_index
 * from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the
 * sigma matrix basis.
 *
 * There are Nc*(Nc-1)/2 unique SU(2) submatrices in an SU(Nc) matrix.
 * The user does not need to know exactly which one is which, just that
 * they are unique.
 *
 * Arguments:
 *
 * \param source      extract the su2 matrix from this su(n) gauge field
 * \param su2_index   int lying in [0, Nc*(Nc-1)/2)
 * \param s           subset for operations       (Read)
 *
 * \return su2 matrix represented in the O(4) rep. - an array of LatticeReal 
 */

multi1d<LatticeReal> 
su2Extract(const LatticeColorMatrix& source, 
	   int su2_index, 
	   const OrderedSubset& s);

#endif
