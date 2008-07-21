// $Id: key_prop_colorvec.cc,v 1.1 2008-07-21 02:32:24 edwards Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#include "util/ferm/key_prop_colorvec.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Support for the keys of prop color vectors
  bool operator<(const KeyPropColorVec_t& a, const KeyPropColorVec_t& b)
  {
    multi1d<int> lgaa(3);
    lgaa[0] = a.t_source;
    lgaa[1] = a.colorvec_src;
    lgaa[2] = a.spin_src;

    multi1d<int> lgbb(3);
    lgbb[0] = b.t_source;
    lgbb[1] = b.colorvec_src;
    lgbb[2] = b.spin_src;

    return (lgaa < lgbb);
  }

  //----------------------------------------------------------------------------
  // KeyPropColorVec read
  void read(BinaryReader& bin, KeyPropColorVec_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.colorvec_src);
    read(bin, param.spin_src);
  }

  // KeyPropColorVec write
  void write(BinaryWriter& bin, const KeyPropColorVec_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.colorvec_src);
    write(bin, param.spin_src);
  }

  //! KeyPropColorVec reader
  void read(XMLReader& xml, const std::string& path, KeyPropColorVec_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "colorvec_src", param.colorvec_src);
    read(paramtop, "spin_src", param.spin_src);
  }

  // KeyPropColorVec writer
  void write(XMLWriter& xml, const std::string& path, const KeyPropColorVec_t& param)
  {
    push(xml, path);

    write(xml, "t_source", param.t_source);
    write(xml, "colorvec_src", param.colorvec_src);
    write(xml, "spin_src", param.spin_src);

    pop(xml);
  }

} // namespace Chroma
