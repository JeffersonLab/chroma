// -*- C++ -*-
/*! \file
 * \brief Key for vanilla distillation propagator sources and solutions
 */

#ifndef __key_prop_distillation_h__
#define __key_prop_distillation_h__

#include "chromabase.h"

namespace Chroma
{

  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */

  //! Distillation propagators
  struct KeyPropDistillation_t
  {
    int                t_source;      /*!< Propagator source time slice */
    int                t_slice;       /*!< Propagator sink time slice */
    int                colorvec_src;  /*!< Source colovec index */
    int                spin_src;      /*!< Source spin index */
    int                spin_snk;      /*!< Sink spin index */
    std::string        mass;          /*!< Quark mass label */
  };


  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDistillation_t& a);

  //----------------------------------------------------------------------------
  //! KeyPropDist read
  void read(BinaryReader& bin, KeyPropDistillation_t& param);

  //! KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDistillation_t& param);

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDistillation_t& param);

  //! KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDistillation_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
