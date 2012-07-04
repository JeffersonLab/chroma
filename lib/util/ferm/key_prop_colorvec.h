// -*- C++ -*-
/*! \file
 * \brief Key for propagator colorvector sources
 */

#ifndef __key_prop_colorvec_h__
#define __key_prop_colorvec_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Prop operator
  /*! \ingroup ferm */
  struct KeyPropColorVec_t
  {
    int        t_source;      /*!< Source time slice */
    int        colorvec_src;  /*!< Source colorvector index */
    int        spin_src;      /*!< Source spin index */
  };


  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! KeyPropColorVec read
  void read(BinaryReader& bin, KeyPropColorVec_t& param);

  //! KeyPropColorVec write
  void write(BinaryWriter& bin, const KeyPropColorVec_t& param);

  //! KeyPropColorVec reader
  void read(XMLReader& xml, const std::string& path, KeyPropColorVec_t& param);

  //! KeyPropColorVec writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropColorVec_t& param);
  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
