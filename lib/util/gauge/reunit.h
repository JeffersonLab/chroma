// -*- C++ -*-
// $Id: reunit.h,v 1.3 2003-03-31 19:46:08 edwards Exp $

#ifndef REUNIT_INCLUDE
#define REUNIT_INCLUDE

enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};

//! Reunitarize in place a color matrix to SU(N)
/*!
 * \ingroup gauge
 *
 * \param  a  The descriptor of matrices to be reunitarized.
 *            Must be of type LatticeColorMatrix
 */
void reunit(LatticeColorMatrix& xa);

//! Reunitarize in place a color matrix to SU(N)
/*!
 * \ingroup gauge
 *
 * \param  a  The descriptor of matrices to be reunitarized.
 *            Must be of type LatticeColorMatrix
 * \param bad Descriptor of flags indicating sites violating unitarity.
 *            Only used if ruflag = REUNITARIZE_LABEL or
 *            REUNITARIZE_ERROR.
 * \param ruflag Can also be REUNITARIZE in which case the
 *            matrices are reunitarized but no complaints are made.
 * \param numbad Total number of matrices violating unitarity.
 *            ONLY USED IF ruflag is testing for ERROR or LABEL. 
 */
void reunit(LatticeColorMatrix& xa, LatticeBoolean& bad, int& numbad, enum Reunitarize ruflag);

#endif
