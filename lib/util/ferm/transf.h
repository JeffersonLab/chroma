// -*- C++ -*-
// $Id: transf.h,v 1.2 2004-04-14 19:20:18 edwards Exp $

#ifndef __transf_h__
#define __transf_h__

//! Convert (insert) a LatticeColorVector into a LatticeFermion
/*!
 * \ingroup ferm
 *
 * \param a      Source color vector
 * \param b      Destination fermion
 * \param spin_index   Spin index
 */
void CvToFerm(const LatticeColorVector& a, LatticeFermion& b, 
	      int spin_index);


//! Insert a LatticeFermion into a LatticePropagator
/*!
 * \ingroup ferm
 *
 * \param a      Source fermion
 * \param b      Destination propagator
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void FermToProp(const LatticeFermion& a, LatticePropagator& b, 
		int color_index, int spin_index);


//! Extract a LatticeFermion from a LatticePropagator
/*!
 * \ingroup ferm
 *
 * \param a      Source propagator
 * \param b      Destination fermion
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void PropToFerm(const LatticePropagator& b, LatticeFermion& a, 
		int color_index, int spin_index);

#endif


