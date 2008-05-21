// $Id: aniso_sym_gaugeact_params.cc,v 3.2 2008-05-21 17:07:49 bjoo Exp $
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#include "chromabase.h"

#include "actions/gauge/gaugeacts/aniso_sym_gaugeact_params.h"

namespace Chroma
{
 
  
  AnisoSymGaugeActParams::AnisoSymGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try 
    {
      read(paramtop, "beta", beta);
      read(paramtop, "u_s", u_s);
      read(paramtop, "u_t", u_t);
      read(paramtop, "AnisoParam", aniso);
      if( paramtop.count("ZeroEnergy") == 1 ) {
	read(paramtop, "ZeroEnergy", sub_zero);
	use_subtraction = true;
      }
      else { 
	sub_zero = Real(0);
	use_subtraction = false;
      }
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, AnisoSymGaugeActParams& p) 
  {
    AnisoSymGaugeActParams tmp(xml, path);
    p=tmp;
  }

  void write(XMLWriter& xml, const string& path, const AnisoSymGaugeActParams& param) 
  {
    push(xml, path);

    write(xml, "beta", param.beta);
    write(xml, "u_s", param.u_s);
    write(xml, "u_t", param.u_t);
    if( param.use_subtraction ) { 
      write(xml, "ZeroEnergy", param.sub_zero);
    }
    write(xml, "AnisoParam", param.aniso);

    pop(xml);
  }

  AnisoSymSpatialGaugeActParams::AnisoSymSpatialGaugeActParams(XMLReader& xml_in, const std::string& path) 
  {
    XMLReader paramtop(xml_in, path);

    try 
    {
      read(paramtop, "beta", beta);
      read(paramtop, "u_s", u_s);
      if( paramtop.count("ZeroEnergy") == 1 ) {
	read(paramtop, "ZeroEnergy", sub_zero);
	use_subtraction = true;
      }
      else { 
	sub_zero = Real(0);
	use_subtraction = false;
      }
      read(paramtop, "AnisoParam", aniso);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " <<  e << endl;
      QDP_abort(1);
    }
  }


  void read(XMLReader& xml, const string& path, AnisoSymSpatialGaugeActParams& p) 
  {
    AnisoSymSpatialGaugeActParams tmp(xml, path);
    p=tmp;
  }

  void write(XMLWriter& xml, const string& path, const AnisoSymSpatialGaugeActParams& param) 
  {
    push(xml, path);

    write(xml, "beta", param.beta);
    write(xml, "u_s", param.u_s);
    if( param.use_subtraction ) { 
      write(xml, "ZeroEnergy", param.sub_zero);
    }
    write(xml, "AnisoParam", param.aniso);

    pop(xml);
  }

}

