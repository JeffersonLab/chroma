// $Id: smearing_io.cc,v 2.0 2005-09-25 21:04:32 edwards Exp $
/*! \file
 *  \brief Smearing parameters
 */

#include "chromabase.h"
#include "io/smearing_io.h"


namespace Chroma 
{

  //! Initialize a smearing param struct
  SmearingParam_t::SmearingParam_t()
  {
    wvf_kind    = WVF_KIND_GAUSSIAN;
    wvf_param   = 0;
    wvfIntPar   = 0;
  }

  //! Read a smearing param struct
  void read(XMLReader& xml, const string& path, SmearingParam_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "wvf_kind", param.wvf_kind);
    read(paramtop, "wvf_param", param.wvf_param);
    read(paramtop, "wvfIntPar", param.wvfIntPar);
  }

  //! Write a smearing param struct
  void write(XMLWriter& xml, const string& path, const SmearingParam_t& param)
  {
    push(xml, path);

    write(xml, "wvf_kind", param.wvf_kind);
    write(xml, "wvf_param", param.wvf_param);
    write(xml, "wvfIntPar", param.wvfIntPar);

    pop(xml);
  }

}  // end namespace Chroma
