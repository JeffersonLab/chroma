// -*- C++ -*-
/*! \file
 * \brief Key for distillution propagator sources and solutions
 */

#ifndef __key_prop_distillution_h__
#define __key_prop_distillution_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */

  //! Distillution propagators
  struct KeyPropDist_t
  {
    std::string        prop_type;     /*!< Distillution source/solution type */
    bool               annihP;        /*!< An annihilation line? */
    int                t_source;      /*!< Propagator source time slice */
    int                t_slice;       /*!< Propagator sink time slice */
    int                dist_src;      /*!< Source dist index */
    int                spin_src;      /*!< Source spin index */
    int                spin_snk;      /*!< Sink spin index */
    int                quark_line;    /*!< The quark line */
    std::string        mass;          /*!< Quark mass label */
  };


  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDist_t& a);

  //----------------------------------------------------------------------------
  //! KeyPropDist read
  void read(BinaryReader& bin, KeyPropDist_t& param);

  //! KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDist_t& param);

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDist_t& param);

  //! KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDist_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
