// -*- C++ -*-
// $Id: axgauge.h,v 1.2 2005-01-14 18:42:35 edwards Exp $
/*! \file
 *  \brief Axial gauge fixing 
 */

#ifndef __axgauge_h__
#define __axgauge_h__

namespace Chroma {

//! Axial gauge fixing
/*!
 * \ingroup gfix
 *
 * Transfroms a gauge configuration, in place, into axial
 *            gauge with special direction j_decay.
 *
 * Note: The non-unity time-like gauge fields from the last time slice
 *       will be copied to all other time slices.
 *
 * \param ug         gauge field and its axial gauge transform (Modify)
 * \param j_decay    time direction (Read) 
 */

void axGauge(multi1d<LatticeColorMatrix>& ug, int j_decay);

}  // end namespace Chroma

#endif
