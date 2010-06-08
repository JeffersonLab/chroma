// -*- C++ -*-
/*! \file
 * \brief Key for glueball colorvector matrix elements
 */

#include "util/gauge/key_glue_matelem.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! GlueElementalOperator reader
  void read(BinaryReader& bin, KeyGlueElementalOperator_t& param)
  {
    read(bin, param.t_slice);
    read(bin, param.left);
    read(bin, param.right);
    read(bin, param.displacement);
    read(bin, param.mom);
  }

  //! GlueElementalOperator write
  void write(BinaryWriter& bin, const KeyGlueElementalOperator_t& param)
  {
    write(bin, param.t_slice);
    write(bin, param.left);
    write(bin, param.right);
    write(bin, param.displacement);
    write(bin, param.mom);
  }

  //! GlueElementalOperator reader
  void read(XMLReader& xml, const std::string& path, KeyGlueElementalOperator_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "left", param.left);
    read(paramtop, "right", param.right);
    read(paramtop, "displacement", param.displacement);
    read(paramtop, "mom", param.mom);
  }

  //! GlueElementalOperator writer
  void write(XMLWriter& xml, const std::string& path, const KeyGlueElementalOperator_t& param)
  {
    push(xml, path);

    write(xml, "t_slice", param.t_slice);
    write(xml, "left", param.left);
    write(xml, "right", param.right);
    write(xml, "displacement", param.displacement);
    write(xml, "mom", param.mom);

    pop(xml);
  }


  //----------------------------------------------------------------------------
  //! GlueElementalOperator reader
  void read(BinaryReader& bin, ValGlueElementalOperator_t& param)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    param.op.resize(n);
  
    for(int i=0; i < param.op.size(); ++i)
    {
      read(bin, param.op[i]);
    }
  }

  //! GlueElementalOperator write
  void write(BinaryWriter& bin, const ValGlueElementalOperator_t& param)
  {
    write(bin, param.op.size());    // always write the size

    for(int i=0; i < param.op.size(); ++i)
    {
      write(bin, param.op[i]);
    }
  }

} // namespace Chroma
