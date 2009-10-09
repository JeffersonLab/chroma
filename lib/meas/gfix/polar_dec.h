// -*- C++ -*-
// $Id: polar_dec.h,v 3.2 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Decompose a complex matrix as \f$C = exp(i\alpha) V P\f$
 */

#ifndef __polar_dec_h__
#define __polar_dec_h__

namespace Chroma {
//! Decompose a complex matrix as \f$C = exp(i\alpha) V P\f#
/*!
 * \ingroup gfix
 *
 * Decompose a complex matrix as \f$C = exp(i\alpha) V P\f#
 * with V SU(Nc) and \f#P = (C^\dagger C)^{1/2}\f# positive hermitian
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

}
#endif
