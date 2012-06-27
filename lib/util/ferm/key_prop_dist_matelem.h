// -*- C++ -*-
/*! \file
 * \brief Key for propagator distillution matrix elements
 */

#ifndef __key_prop_dist_matelem_h__
#define __key_prop_dist_matelem_h__

#include "chromabase.h"
#include "util/ferm/key_val_db.h"

namespace Chroma
{
  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */
  //! Prop operator
  struct KeyPropDistElemOp_t
  {
    int                quark_line;    /*!< The quark line */
    bool               annihP;        /*!< An annihilation line? */
    int                t_slice;       /*!< Propagator time slice */
    int                t_source;      /*!< Source time slice */
    int                spin_src;      /*!< Source spin index */
    int                spin_snk;      /*!< Sink spin index */
    std::string        mass;          /*!< A mass label */
  };


  //! Prop operator
  struct ValPropDistElemOp_t
  {
    multi2d<ComplexD>  mat;               /*!< Distillution source and colorvector sink */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValPropDistElemOp_t
  {
    SerialDBKey<KeyPropDistElemOp_t>  key;
    SerialDBData<ValPropDistElemOp_t> val;
  };

  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDistElemOp_t& a);

  //----------------------------------------------------------------------------
  //! PropDistElemOp reader
  void read(BinaryReader& bin, KeyPropDistElemOp_t& param);

  //! PropDistElemOp write
  void write(BinaryWriter& bin, const KeyPropDistElemOp_t& param);

  //! PropDistElemOp reader
  void read(XMLReader& xml, const std::string& path, KeyPropDistElemOp_t& param);

  //! PropDistElemOp writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDistElemOp_t& param);


  //----------------------------------------------------------------------------
  //! PropDistElemOp reader
  void read(BinaryReader& bin, ValPropDistElemOp_t& param);

  //! PropDistElemOp write
  void write(BinaryWriter& bin, const ValPropDistElemOp_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
