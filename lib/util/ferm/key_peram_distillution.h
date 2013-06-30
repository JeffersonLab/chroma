// -*- C++ -*-
/*! \file
 * \brief Key for propagator distillution matrix elements
 */

#ifndef __key_peram_distillution_h__
#define __key_peram_distillution_h__

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
  struct KeyPeramDistillution_t
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
  struct ValPeramDistillution_t
  {
    multi2d<ComplexD>  mat;               /*!< Distillution source and colorvector sink */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValPeramDistillution_t
  {
    SerialDBKey<KeyPeramDistillution_t>  key;
    SerialDBData<ValPeramDistillution_t> val;
  };

  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPeramDistillution_t& a);

  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, KeyPeramDistillution_t& param);

  //! PeramDist write
  void write(BinaryWriter& bin, const KeyPeramDistillution_t& param);

  //! PeramDist reader
  void read(XMLReader& xml, const std::string& path, KeyPeramDistillution_t& param);

  //! PeramDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPeramDistillution_t& param);


  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, ValPeramDistillution_t& param);

  //! PeramDist write
  void write(BinaryWriter& bin, const ValPeramDistillution_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
