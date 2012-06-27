/*! \file
 * \brief Key for distillution propagator sources and solutions
 */

#include <vector>
#include "util/ferm/key_prop_distillution.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDist_t& param)
  {
    os << "KeyPropDist_t:";
    os << " prop_type= " << param.prop_type;
    os << " annihP= " << param.annihP;
    os << " t_source= " << param.t_source;
    os << " t_slice= " << param.t_slice;
    os << " spin_snk= " << param.spin_snk;
    os << " spin_src= " << param.spin_src;
    os << " dist_src= " << param.dist_src;
    os << " quark_line= " << param.quark_line;
    os << " mass= " << param.mass;
    os << endl;

    return os;
  }

  //----------------------------------------------------------------------------
  // KeyPropDist read
  void read(BinaryReader& bin, KeyPropDist_t& param)
  {
    readDesc(bin, param.prop_type);
    read(bin, param.annihP);
    read(bin, param.t_source);
    read(bin, param.t_slice);
    read(bin, param.dist_src);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    read(bin, param.quark_line);
    readDesc(bin, param.mass);
  }

  // KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDist_t& param)
  {
    writeDesc(bin, param.prop_type);
    write(bin, param.annihP);
    write(bin, param.t_source);
    write(bin, param.t_slice);
    write(bin, param.dist_src);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    write(bin, param.quark_line);
    writeDesc(bin, param.mass);
  }

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDist_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "prop_type", param.prop_type);
    read(paramtop, "annihP", param.annihP);
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "dist_src", param.dist_src);
    read(paramtop, "spin_src", param.spin_src);
    read(paramtop, "spin_snk", param.spin_snk);
    read(paramtop, "quark_line", param.quark_line);
    read(paramtop, "mass", param.mass);
  }

  // KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDist_t& param)
  {
    push(xml, path);

    write(xml, "prop_type", param.prop_type);
    write(xml, "annihP", param.annihP);
    write(xml, "t_source", param.t_source);
    write(xml, "t_slice", param.t_slice);
    write(xml, "dist_src", param.dist_src);
    write(xml, "spin_src", param.spin_src);
    write(xml, "spin_snk", param.spin_snk);
    write(xml, "quark_line", param.quark_line);
    write(xml, "mass", param.mass);

    pop(xml);
  }

} // namespace Chroma
