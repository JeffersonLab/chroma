// $Id: key_grid_prop.cc,v 1.1 2008-11-27 03:22:17 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#include "util/ferm/key_grid_prop.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  // Support for the keys of prop color vectors
  bool operator<(const KeyGridProp_t& a, const KeyGridProp_t& b)
  {
    multi1d<int> lgaa(4);
    lgaa[0] = a.t_source;
    lgaa[1] = a.spin ;
    lgaa[2] = a.color ;
    lgaa[3] = a.grid ;

    multi1d<int> lgbb(4);
    lgbb[0] = b.t_source;
    lgbb[1] = b.spin ;
    lgbb[2] = b.color ;
    lgbb[3] = b.grid ;

    return (lgaa < lgbb);
  }

  //--------------------------------------------------------------------------
  // KeyGridProp read
  void read(BinaryReader& bin, KeyGridProp_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.spin);
    read(bin, param.color);
    read(bin, param.grid);
  }

  // KeyGridProp write
  void write(BinaryWriter& bin, const KeyGridProp_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.spin    );
    write(bin, param.color   );
    write(bin, param.grid    );
  }

  //! KeyGridProp reader
  void read(XMLReader& xml, const std::string& path, KeyGridProp_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "t_source", param.t_source);
    read(paramtop, "color"   , param.color);
    read(paramtop, "spin"    , param.spin);
    read(paramtop, "grid"    , param.grid);
  }

  // KeyGridProp writer
  void write(XMLWriter& xml,const std::string& path,const KeyGridProp_t& param)
  {
    push(xml, path);

    write(xml, "t_source", param.t_source);
    write(xml, "color"   , param.color   );
    write(xml, "spin"    , param.spin    );
    write(xml, "grid"    , param.grid    );

    pop(xml);
  }

} // namespace Chroma
