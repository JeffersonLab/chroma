// -*- C++ -*-
// $Id: unit_check.h,v 1.2 2005-01-14 15:59:00 bjoo Exp $

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
