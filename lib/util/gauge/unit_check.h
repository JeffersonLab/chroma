// -*- C++ -*-
// $Id: unit_check.h,v 1.1 2003-10-02 03:26:21 edwards Exp $

#ifndef UNIT_CHECK_INCLUDE
#define UNIT_CHECK_INCLUDE

//! Check the unitarity of color matrix in SU(N)
/*!
 * \ingroup gauge
 *
 * \param  u  The multi1d<LatticeColorMatrix> to be tested
 */
void unitarityCheck(const multi1d<LatticeColorMatrix>& u);

#endif
