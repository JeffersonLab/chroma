// $Id: ptsrc_params.cc,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Point source construction
 *
 * This code is common to different fermions and propagators
 */

#include "chromabase.h"

#include "meas/sources/ptsrc_params.h"

namespace Chroma
{
  //! Read parameters
  PointSourceConstParams::PointSourceConstParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "t_source", t_source);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, PointSourceConstParams& param)
  {
    PointSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointSourceConstParams& param)
  {
    push(xml, path);
    write(xml, "t_source", param.t_source);
    pop(xml);
  }

}
