// -*- C++ -*-
// $Id: temporal_gauge.h,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Axial gauge fixing 
 */

#ifndef __temporal_gauge_h__
#define __temporal_gauge_h__

namespace Chroma 
{

  //! Temporal gauge fixing
  /*!
   * \ingroup gfix
   *
   * Transfroms a gauge configuration, in place, into axial
   *            gauge with special direction decay_dir.
   *
   * Note: The non-unity time-like gauge fields from the last time slice
   *       will be copied to all other time slices.
   *
   * \param ug         gauge field and its axial gauge transform (Modify)
   * \param g          gauge rotation matrix (Write)
   * \param decay_dir  time direction (Read) 
   */
  void temporalGauge(multi1d<LatticeColorMatrix>& ug, LatticeColorMatrix& g, int decay_dir);

}  // end namespace Chroma

#endif
