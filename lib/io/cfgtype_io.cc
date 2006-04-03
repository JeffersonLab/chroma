// $Id: cfgtype_io.cc,v 3.0 2006-04-03 04:58:55 edwards Exp $
/*! \file
 *  \brief Configuration structure IO
 */

#include "chromabase.h"
#include "io/cfgtype_io.h"


namespace Chroma 
{

  // Configuration input
  void read(XMLReader& xml, const string& path, Cfg_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "cfg_type", input.cfg_type);
    read(inputtop, "cfg_file", input.cfg_file);
  }

  // Write a config struct
  void write(XMLWriter& xml, const string& path, const Cfg_t& cfg)
  {
    push(xml, "Cfg");
    write(xml, "cfg_type", cfg.cfg_type);
    write(xml, "cfg_file", cfg.cfg_file);
    pop(xml);
  }

}  // end namespace Chroma
