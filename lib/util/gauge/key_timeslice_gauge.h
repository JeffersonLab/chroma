// -*- C++ -*-
/*! \file
 * \brief Key for time-sliced gauge fields
 */

#ifndef __key_timeslice_gauge_h__
#define __key_timeslice_gauge_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Prop operator
  /*! \ingroup ferm */
  struct KeyTimeSliceGauge_t
  {
    int        t_slice;       /*!< Source time slice */
    int        dir;           /*!< Direction */
  };


  //! Support for the keys of prop color vectors
  /*! \ingroup ferm */
  bool operator<(const KeyTimeSliceGauge_t& a, const KeyTimeSliceGauge_t& b);

  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! KeyTimeSliceGauge read
  void read(BinaryReader& bin, KeyTimeSliceGauge_t& param);

  //! KeyTimeSliceGauge write
  void write(BinaryWriter& bin, const KeyTimeSliceGauge_t& param);

  //! KeyTimeSliceGauge reader
  void read(XMLReader& xml, const std::string& path, KeyTimeSliceGauge_t& param);

  //! KeyTimeSliceGauge writer
  void write(XMLWriter& xml, const std::string& path, const KeyTimeSliceGauge_t& param);
  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
