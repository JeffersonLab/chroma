// -*- C++ -*-
/*! \file
 * \brief Key for propagator colorvector matrix elements
 */

#include "util/ferm/key_prop_matelem.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! PropElementalOperator reader
  void read(BinaryReader& bin, KeyPropElementalOperator_t& param)
  {
    read(bin, param.t_slice);
    read(bin, param.t_source);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    read(bin, param.mass_label, 32);
  }

  //! PropElementalOperator write
  void write(BinaryWriter& bin, const KeyPropElementalOperator_t& param)
  {
    write(bin, param.t_slice);
    write(bin, param.t_source);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    write(bin, param.mass_label);
  }

  //! PropElementalOperator reader
  void read(XMLReader& xml, const std::string& path, KeyPropElementalOperator_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "spin_src", param.spin_src);
    read(paramtop, "spin_snk", param.spin_snk);
    read(paramtop, "mass_label", param.mass_label);
  }

  //! PropElementalOperator writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropElementalOperator_t& param)
  {
    push(xml, path);

    write(xml, "t_slice", param.t_slice);
    write(xml, "t_source", param.t_source);
    write(xml, "spin_src", param.spin_src);
    write(xml, "spin_snk", param.spin_snk);
    write(xml, "mass_label", param.mass_label);

    pop(xml);
  }


  //----------------------------------------------------------------------------
  //! PropElementalOperator reader
  void read(BinaryReader& bin, ValPropElementalOperator_t& param)
  {
    int n1;
    int n2;
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n1);    // the size is always written, even if 0
    param.mat.resize(n2,n1);
  
    for(int i=0; i < param.mat.size1(); ++i)
    {
      for(int j=0; j < param.mat.size2(); ++j)
      {
	read(bin, param.mat[j][i]);
      }
    }
  }

  //! PropElementalOperator write
  void write(BinaryWriter& bin, const ValPropElementalOperator_t& param)
  {
    write(bin, param.mat.size2());    // always write the size
    write(bin, param.mat.size1());    // always write the size

    for(int i=0; i < param.mat.size1(); ++i)
    {
      for(int j=0; j < param.mat.size2(); ++j)
      {
	write(bin, param.mat[j][i]);
      }
    }
  }

} // namespace Chroma
