// -*- C++ -*-
// $Id: unit_check.h,v 3.0 2006-04-03 04:59:13 edwards Exp $

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
