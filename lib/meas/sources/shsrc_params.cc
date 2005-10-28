// $Id: shsrc_params.cc,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Shell source construction
 *
 * This code is common to different fermions and propagators
 */

#include "chromabase.h"

#include "meas/sources/shsrc_params.h"

namespace Chroma
{
  //! Read parameters
  ShellSourceConstParams::ShellSourceConstParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "SourceType",  source_type);
    read(paramtop, "wave_state", wave_state);
    if (wave_state != WAVE_TYPE_S_WAVE)
    {
      read(paramtop, "direction",  direction);
    }

    {
      XMLReader shelltop(paramtop, "ShellSource");

      read(shelltop, "SourceSmearingParam", sourceSmearParam);
      read(shelltop, "laplace_power", laplace_power);
      read(shelltop, "disp_length", disp_length);
      read(shelltop, "disp_dir", disp_dir);

      {
	XMLReader xml_tmp(shelltop, "LinkSmearing");
	std::ostringstream os;
	xml_tmp.print(os);
	link_smearing = os.str();
      }
    }		

    read(paramtop, "j_decay",  j_decay);
    read(paramtop, "t_source", t_source);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, ShellSourceConstParams& param)
  {
    ShellSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const ShellSourceConstParams& param)
  {
    push(xml, path);

//    int version = 6;
//    write(xml, "version", version);
    write(xml, "wave_state", param.wave_state);
    if (param.wave_state != WAVE_TYPE_S_WAVE)
    {
      write(xml, "direction",  param.direction);
    }

    write(xml, "source_type", param.source_type);
    {
      push(xml, "ShellSource");
      write(xml, "SourceSmearingParam", param.sourceSmearParam);
      write(xml, "laplace_power", param.laplace_power);
      write(xml, "disp_length", param.disp_length);
      write(xml, "disp_dir", param.disp_dir);
      xml << param.link_smearing;
      pop(xml);
    }

    write(xml, "j_decay",  param.j_decay);
    write(xml, "t_source",  param.t_source);

    pop(xml);
  }

}
