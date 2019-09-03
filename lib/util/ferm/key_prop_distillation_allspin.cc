/*! \file
 * \brief Key for vanilla distillationAllSpin propagator sources and solutions
 */

#include "util/ferm/key_prop_distillation_allspin.h"
#include <vector>

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //! Diagnostics
  StandardOutputStream& operator<<(StandardOutputStream& os, const KeyPropDistillationAllSpin_t& param)
  {
    os << "KeyPropDistillationAllSpin_t:";
    os << " t_slice= " << param.t_slice;
    os << " t_source= " << param.t_source;
    os << " colorvec_src= " << param.colorvec_src;
    os << " mass= " << param.mass;
    os << std::endl;

    return os;
  }

  //----------------------------------------------------------------------------
  // KeyPropDist read
  void read(BinaryReader& bin, KeyPropDistillationAllSpin_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.t_slice);
    read(bin, param.colorvec_src);
    readDesc(bin, param.mass);
  }

  // KeyPropDist write
  void write(BinaryWriter& bin, const KeyPropDistillationAllSpin_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.t_slice);
    write(bin, param.colorvec_src);
    writeDesc(bin, param.mass);
  }

  //! KeyPropDist reader
  void read(XMLReader& xml, const std::string& path, KeyPropDistillationAllSpin_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "t_slice", param.t_slice);
    read(paramtop, "colorvec_src", param.colorvec_src);
    read(paramtop, "mass", param.mass);
  }

  // KeyPropDist writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropDistillationAllSpin_t& param)
  {
    push(xml, path);

    write(xml, "t_source", param.t_source);
    write(xml, "t_slice", param.t_slice);
    write(xml, "colorvec_src", param.colorvec_src);
    write(xml, "mass", param.mass);

    pop(xml);
  }

} // namespace Chroma
