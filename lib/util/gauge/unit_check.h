// -*- C++ -*-
// $Id: unit_check.h,v 2.0 2005-09-25 21:04:45 edwards Exp $

#ifndef UNIT_CHECK_INCLUDE
#define UNIT_CHECK_INCLUDE

namespace Chroma {
//! Check the unitarity of color matrix in SU(N)
/*!
 * \ingroup gauge
 *
 * \param  u  The multi1d<LatticeColorMatrix> to be tested
 */
void unitarityCheck(const multi1d<LatticeColorMatrix>& u);

};

#endif
