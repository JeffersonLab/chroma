// -*- C++ -*-
/*! \file
 * \brief Key for time-sliced color eigenvectors
 */

#ifndef __key_timeslice_colorvec_h__
#define __key_timeslice_colorvec_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Prop operator
  /*! \ingroup ferm */
  struct KeyTimeSliceColorVec_t
  {
    int        t_slice;       /*!< Source time slice */
    int        colorvec;      /*!< Colorvector index */
  };


  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyTimeSliceColorVec_t& param);

  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! KeyTimeSliceColorVec read
  void read(BinaryReader& bin, KeyTimeSliceColorVec_t& param);

  //! KeyTimeSliceColorVec write
  void write(BinaryWriter& bin, const KeyTimeSliceColorVec_t& param);

  //! KeyTimeSliceColorVec reader
  void read(XMLReader& xml, const std::string& path, KeyTimeSliceColorVec_t& param);

  //! KeyTimeSliceColorVec writer
  void write(XMLWriter& xml, const std::string& path, const KeyTimeSliceColorVec_t& param);
  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
