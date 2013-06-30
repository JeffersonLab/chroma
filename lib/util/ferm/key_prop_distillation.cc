/*! \file
 * \brief Key for vanilla distillation propagator sources and solutions
 */

#include <vector>
#include "util/ferm/key_prop_distillation.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDistillation_t& param)
  {
    os << "KeyPropDistillation_t:";
    os << " t_slice= " << param.t_slice;
    os << " t_source= " << param.t_source;
    os << " spin_snk= " << param.spin_snk;
    os << " spin_src= " << param.spin_src;
    os << " colorvec_src= " << param.colorvec_src;
    os << " mass= " << param.mass;
    os << endl;

    return os;
  }

  //----------------------------------------------------------------------------
  // KeyPropDist read
  void read(BinaryReader& bin, KeyPropDistillation_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.t_slice);
    read(bin, param.colorvec_src);
    read(bin, param.spin_src);
    read(bin, param.spin_snk);
    readDesc(bin, param.mass);
  }

  // KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDistillation_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.t_slice);
    write(bin, param.colorvec_src);
    write(bin, param.spin_src);
    write(bin, param.spin_snk);
    writeDesc(bin, param.mass);
  }

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDistillation_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "colorvec_src", param.colorvec_src);
    read(paramtop, "spin_src", param.spin_src);
    read(paramtop, "spin_snk", param.spin_snk);
    read(paramtop, "mass", param.mass);
  }

  // KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDistillation_t& param)
  {
    push(xml, path);

    write(xml, "t_source", param.t_source);
    write(xml, "t_slice", param.t_slice);
    write(xml, "colorvec_src", param.colorvec_src);
    write(xml, "spin_src", param.spin_src);
    write(xml, "spin_snk", param.spin_snk);
    write(xml, "mass", param.mass);

    pop(xml);
  }

} // namespace Chroma
