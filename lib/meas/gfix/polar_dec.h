// -*- C++ -*-
// $Id: polar_dec.h,v 1.1 2003-12-06 20:56:56 edwards Exp $
/*! \file
 *  \brief Decompose a complex matrix as C = exp(i\alpha) V P
 */

#ifndef __polar_dec_h__
#define __polar_dec_h__

//! Decompose a complex matrix as C = exp(i\alpha) V P
/*!
 * \ingroup gfix
 *
 * Decompose a complex matrix as C = exp(i\alpha) V P
 * with V SU(Nc) and P = (C^\dagger C)^{1/2} positive hermitian
 *
 * \param c        complex Nc x Nc matrix ( Modify )
 *                 on exit it contains the hermitian matrix P
 * \param v        the projected SU(Nc) Matrix ( Write )
 * \param alpha    the phase ( Write )
 * \param JacAccu  accuracy in the Jacobi iteration ( Read )
 * \param JacMax   maximum number of Jacobi iterations ( Read ) 
 */

void polar_dec(LatticeColorMatrix& c, LatticeColorMatrix& v,
	       LatticeReal& alpha, const Real& JacAccu, int JacMax);

#endif
