// -*- C++ -*-
// $Id: key_block_prop.h,v 1.1 2009-01-30 05:33:24 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#ifndef __key_prop_block_h__
#define __key_prop_block_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Prop operator
  /*! \ingroup ferm */
  struct KeyBlockProp_t
  {
    int        t_source ;      /*!< Source time slice */
    int        spin     ;      /*!< Source spin index */
    int        color    ;      /*!< Source colorvector index */
    int        block     ;      /*!< Source block index */
  };


  //! Support for the keys of prop color vectors
  /*! \ingroup ferm */
  bool operator<(const KeyBlockProp_t& a, const KeyBlockProp_t& b);

  //---------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! KeyBlockProp read
  void read(BinaryReader& bin, KeyBlockProp_t& param);

  //! KeyBlockProp write
  void write(BinaryWriter& bin, const KeyBlockProp_t& param);

  //! KeyBlockProp reader
  void read(XMLReader& xml, const std::string& path, KeyBlockProp_t& param);

  //! KeyBlockProp writer
  void write(XMLWriter& xml, const std::string& path, const KeyBlockProp_t& param);
  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
