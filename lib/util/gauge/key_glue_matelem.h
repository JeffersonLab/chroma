// -*- C++ -*-
/*! \file
 * \brief Key for glueball colorvector matrix elements
 */

#ifndef __key_glue_matelem_h__
#define __key_glue_matelem_h__

#include "chromabase.h"
#include "util/ferm/key_val_db.h"

namespace Chroma
{
  //----------------------------------------------------------------------------
  /*!
   * \ingroup gauge
   * @{
   */
  //! Glue operator
  struct KeyGlueElementalOperator_t
  {
    int                t_slice;      /*!< Glue operator time slice */
    int                left;         /*!< Left B field direction */ 
    int                right;        /*!< Right B field direction */ 
    multi1d<int>       displacement; /*!< Displacement dirs of right colorvector */
    multi1d<int>       mom;          /*!< D-1 momentum of this operator */
  };


  //! Glue operator
  struct ValGlueElementalOperator_t
  {
    multi1d<ComplexD>  op;           /*!< Colorvector source/sink (diagonal) with momentum projection */
  };


  //----------------------------------------------------------------------------
  //! Holds key and value as temporaries
  struct KeyValGlueElementalOperator_t
  {
    SerialDBKey<KeyGlueElementalOperator_t>  key;
    SerialDBData<ValGlueElementalOperator_t> val;
  };

  //----------------------------------------------------------------------------
  //! GlueElementalOperator reader
  void read(BinaryReader& bin, KeyGlueElementalOperator_t& param);

  //! GlueElementalOperator write
  void write(BinaryWriter& bin, const KeyGlueElementalOperator_t& param);

  //! GlueElementalOperator reader
  void read(XMLReader& xml, const std::string& path, KeyGlueElementalOperator_t& param);

  //! GlueElementalOperator writer
  void write(XMLWriter& xml, const std::string& path, const KeyGlueElementalOperator_t& param);


  //----------------------------------------------------------------------------
  //! GlueElementalOperator reader
  void read(BinaryReader& bin, ValGlueElementalOperator_t& param);

  //! GlueElementalOperator write
  void write(BinaryWriter& bin, const ValGlueElementalOperator_t& param);

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
