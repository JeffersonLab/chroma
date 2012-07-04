// -*- C++ -*-
/*! \file
 * \brief Key for propagator distillution matrix elements
 */

#ifndef __key_peram_dist_h__
#define __key_peram_dist_h__

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
  struct KeyPeramDist_t
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
  struct ValPeramDist_t
  {
    multi2d<ComplexD>  mat;               /*!< Distillution source and colorvector sink */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValPeramDist_t
  {
    SerialDBKey<KeyPeramDist_t>  key;
    SerialDBData<ValPeramDist_t> val;
  };

  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPeramDist_t& a);

  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, KeyPeramDist_t& param);

  //! PeramDist write
  void write(BinaryWriter& bin, const KeyPeramDist_t& param);

  //! PeramDist reader
  void read(XMLReader& xml, const std::string& path, KeyPeramDist_t& param);

  //! PeramDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPeramDist_t& param);


  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, ValPeramDist_t& param);

  //! PeramDist write
  void write(BinaryWriter& bin, const ValPeramDist_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
