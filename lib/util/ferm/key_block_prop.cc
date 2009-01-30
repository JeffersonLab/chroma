// $Id: key_block_prop.cc,v 1.1 2009-01-30 05:33:24 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#include "util/ferm/key_block_prop.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Support for the keys of prop color vectors
  bool operator<(const KeyBlockProp_t& a, const KeyBlockProp_t& b)
  {
    multi1d<int> lgaa(4);
    lgaa[0] = a.t_source;
    lgaa[1] = a.spin ;
    lgaa[2] = a.color ;
    lgaa[3] = a.block ;

    multi1d<int> lgbb(4);
    lgbb[0] = b.t_source;
    lgbb[1] = b.spin ;
    lgbb[2] = b.color ;
    lgbb[3] = b.block ;

    return (lgaa < lgbb);
  }

  //--------------------------------------------------------------------------
  // KeyBlockProp read
  void read(BinaryReader& bin, KeyBlockProp_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.spin);
    read(bin, param.color);
    read(bin, param.block);
  }

  // KeyBlockProp write
  void write(BinaryWriter& bin, const KeyBlockProp_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.spin    );
    write(bin, param.color   );
    write(bin, param.block    );
  }

  //! KeyBlockProp reader
  void read(XMLReader& xml, const std::string& path, KeyBlockProp_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "color"   , param.color);
    read(paramtop, "spin"    , param.spin);
    read(paramtop, "block"    , param.block);
  }

  // KeyBlockProp writer
  void write(XMLWriter& xml,const std::string& path,const KeyBlockProp_t& param)
  {
    push(xml, path);

    write(xml, "t_source", param.t_source);
    write(xml, "color"   , param.color   );
    write(xml, "spin"    , param.spin    );
    write(xml, "block"    , param.block    );

    pop(xml);
  }

} // namespace Chroma
