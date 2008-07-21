// -*- C++ -*-
// $Id: key_prop_colorvec.h,v 1.1 2008-07-21 02:32:24 edwards Exp $
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


  //! Support for the keys of prop color vectors
  /*! \ingroup ferm */
  bool operator<(const KeyPropColorVec_t& a, const KeyPropColorVec_t& b);

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
