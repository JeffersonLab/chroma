// -*- C++ -*-
// $Id: expmat.h,v 1.2 2004-01-05 00:47:20 edwards Exp $
/*! \file
 *  \brief Exponentiate a SU(n) lie algebra element by some method,
 */

#ifndef __expmat_h__
#define __expmat_h__

enum ExpMat_t {EXP_EXACT, EXP_TWELTH_ORDER};

//! Exponentiate a SU(n) lie algebra element by some method.
/*!
 * \ingroup gauge
 *
 * This routine is a driver for other routines. For example, expsu3
 * will exponentiate EFFICIENTLY an SU(3) matrix and, if desired,
 * check for errors in the truncation of the Taylor series.
 * The routine eesu3 EXACTLY exponentiates a SU(3) matrix. 
 *
 *  \param a        LatticeColorMatrix          (Modify)
 *  \param opt      Method of exponentiation    (Read)
 */

void expmat(LatticeColorMatrix& a,
	    ExpMat_t opt);

#endif
