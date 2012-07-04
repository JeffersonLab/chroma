// -*- C++ -*-
/*! \file
 * \brief Key for propagator colorvector matrix elements
 */

#ifndef __key_prop_matelem_h__
#define __key_prop_matelem_h__

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
  struct KeyPropElementalOperator_t
  {
    int                t_slice;       /*!< Propagator time slice */
    int                t_source;      /*!< Source time slice */
    int                spin_src;      /*!< Source spin index */
    int                spin_snk;      /*!< Sink spin index */
    std::string        mass_label;    /*!< A mass label */
  };


  //! Prop operator
  struct ValPropElementalOperator_t
  {
    multi2d<ComplexD>  mat;               /*!< Colorvector source and sink */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValPropElementalOperator_t
  {
    SerialDBKey<KeyPropElementalOperator_t>  key;
    SerialDBData<ValPropElementalOperator_t> val;
  };

  //----------------------------------------------------------------------------
  //! PropElementalOperator reader
  void read(BinaryReader& bin, KeyPropElementalOperator_t& param);

  //! PropElementalOperator write
  void write(BinaryWriter& bin, const KeyPropElementalOperator_t& param);

  //! PropElementalOperator reader
  void read(XMLReader& xml, const std::string& path, KeyPropElementalOperator_t& param);

  //! PropElementalOperator writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropElementalOperator_t& param);


  //----------------------------------------------------------------------------
  //! PropElementalOperator reader
  void read(BinaryReader& bin, ValPropElementalOperator_t& param);

  //! PropElementalOperator write
  void write(BinaryWriter& bin, const ValPropElementalOperator_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
