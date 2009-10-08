// -*- C++ -*-
// $Id: axgauge.h,v 3.2 2009-10-08 17:45:37 edwards Exp $
/*! \file
 *  \brief Axial gauge fixing 
 */

#ifndef __axgauge_h__
#define __axgauge_h__

namespace Chroma 
{

  //! Axial gauge fixing
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
   * \param v          gauge rotation matrix (Write)
   * \param decay_dir  time direction (Read) 
   */
  void axGauge(multi1d<LatticeColorMatrix>& ug, LatticeColorMatrix& v, int decay_dir);


  //! Axial gauge fixing
  /*!
   * \ingroup gfix
   *
   * Transfroms a gauge configuration, in place, into axial
   *            gauge with special direction decay_dir.
   *
   * Note: The non-unity time-like gauge fields from the last time slice
   *       will be copied to all other time slices.
   *
   * \param ug           gauge field and its axial gauge transform (Modify)
   * \param decay_dir    time direction (Read) 
   */
  void axGauge(multi1d<LatticeColorMatrix>& ug, int decay_dir);

}  // end namespace Chroma

#endif
