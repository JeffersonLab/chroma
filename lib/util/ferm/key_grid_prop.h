// -*- C++ -*-
// $Id: key_grid_prop.h,v 1.1 2008-11-27 03:22:17 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#ifndef __key_prop_grid_h__
#define __key_prop_grid_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! Prop operator
  /*! \ingroup ferm */
  struct KeyGridProp_t
  {
    int        t_source ;      /*!< Source time slice */
    int        spin     ;      /*!< Source spin index */
    int        color    ;      /*!< Source colorvector index */
    int        grid     ;      /*!< Source grid index */
  };


  //! Support for the keys of prop color vectors
  /*! \ingroup ferm */
  bool operator<(const KeyGridProp_t& a, const KeyGridProp_t& b);

  //---------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! KeyGridProp read
  void read(BinaryReader& bin, KeyGridProp_t& param);

  //! KeyGridProp write
  void write(BinaryWriter& bin, const KeyGridProp_t& param);

  //! KeyGridProp reader
  void read(XMLReader& xml, const std::string& path, KeyGridProp_t& param);

  //! KeyGridProp writer
  void write(XMLWriter& xml, const std::string& path, const KeyGridProp_t& param);
  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
