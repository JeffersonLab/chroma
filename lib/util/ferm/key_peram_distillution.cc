// -*- C++ -*-
/*! \file
 * \brief Key for propagator distillution matrix elements
 */

#include "util/ferm/key_peram_distillution.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPeramDistillution_t& param)
  {
    os << "KeyPeramDistillution_t:";
    os << " quark_line= " << param.quark_line;
    os << " annihP= " << param.annihP;
    os << " t_slice= " << param.t_slice;
    os << " t_source= " << param.t_source;
    os << " spin_snk= " << param.spin_snk;
    os << " spin_src= " << param.spin_src;
    os << " mass= " << param.mass;
    os << endl;

    return os;
  }

  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, KeyPeramDistillution_t& param)
  {
    read(bin, param.quark_line);
    read(bin, param.annihP);
    read(bin, param.t_slice);
    read(bin, param.t_source);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    readDesc(bin, param.mass);
  }

  //! PeramDist write
  void write(BinaryWriter& bin, const KeyPeramDistillution_t& param)
  { 
    write(bin, param.quark_line);
    write(bin, param.annihP);
    write(bin, param.t_slice);
    write(bin, param.t_source);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    writeDesc(bin, param.mass);
  }

  //! PeramDist reader
  void read(XMLReader& xml, const std::string& path, KeyPeramDistillution_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "quark_line", param.quark_line);
    read(paramtop, "annihP", param.annihP);
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "spin_src", param.spin_src);
    read(paramtop, "spin_snk", param.spin_snk);
    read(paramtop, "mass", param.mass);
  }

  //! PeramDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPeramDistillution_t& param)
  {
    push(xml, path);

    write(xml, "quark_line", param.quark_line);
    write(xml, "annihP", param.annihP);
    write(xml, "t_slice", param.t_slice);
    write(xml, "t_source", param.t_source);
    write(xml, "spin_src", param.spin_src);
    write(xml, "spin_snk", param.spin_snk);
    write(xml, "mass", param.mass);

    pop(xml);
  }


  //----------------------------------------------------------------------------
  //! PeramDist reader
  void read(BinaryReader& bin, ValPeramDistillution_t& param)
  {
    int nrows;
    int ncols;
    read(bin, nrows);    // the size is always written, even if 0
    read(bin, ncols);    // the size is always written, even if 0
    param.mat.resize(nrows,ncols);
  
    for(int j=0; j < param.mat.nrows(); ++j)
    {
      for(int i=0; i < param.mat.ncols(); ++i)
      {
	read(bin, param.mat(j,i));
      }
    }
  }

  //! PeramDist write
  void write(BinaryWriter& bin, const ValPeramDistillution_t& param)
  {
    write(bin, param.mat.nrows());    // always write the size
    write(bin, param.mat.ncols());    // always write the size

    for(int j=0; j < param.mat.nrows(); ++j)
    {
      for(int i=0; i < param.mat.ncols(); ++i)
      {
	write(bin, param.mat(j,i));
      }
    }
  }

} // namespace Chroma
