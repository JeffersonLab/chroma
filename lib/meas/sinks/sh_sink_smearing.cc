// $Id: sh_sink_smearing.cc,v 1.5 2005-11-08 18:51:44 edwards Exp $
/*! \file
 *  \brief Shell sink smearing
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sh_sink_smearing.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/quark_smearing.h"
#include "meas/smear/displacement.h"

namespace Chroma
{
  //! Initialize
  ShellQuarkSinkSmearingParams::ShellQuarkSinkSmearingParams()
  {
    disp_length = disp_dir = 0;
  }


  //! Read parameters
  ShellQuarkSinkSmearingParams::ShellQuarkSinkSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default:
      QDPIO::cerr << __func__ << ": parameter version " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "SinkType",  sink_type);

    {
      XMLReader xml_tmp(paramtop, "SmearingParam");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "wvf_kind", quark_smearing_type);
      quark_smearing = os.str();
    }

    if (paramtop.count("disp_length") != 0)
      read(paramtop, "disp_length", disp_length);
    else
      disp_length = 0;

    if (paramtop.count("disp_dir") != 0)
      read(paramtop, "disp_dir", disp_dir);
    else
      disp_dir = 0;

    if (paramtop.count("LinkSmearing") != 0)
    {
      XMLReader xml_tmp(paramtop, "LinkSmearing");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "LinkSmearingType", link_smearing_type);
      link_smearing = os.str();
    }

  }

  // Read parameters
  void read(XMLReader& xml, const string& path, ShellQuarkSinkSmearingParams& param)
  {
    ShellQuarkSinkSmearingParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const ShellQuarkSinkSmearingParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "SinkType", param.sink_type);
    xml << param.quark_smearing;
    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);

    pop(xml);
  }



  //! Hooks to register the class
  namespace ShellQuarkSinkSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
    {
      return new ShellQuarkSinkSmearing<LatticePropagator>(ShellQuarkSinkSmearingParams(xml_in, path), u);
    }

    //! Callback function
    QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						const std::string& path,
						const multi1d<LatticeColorMatrix>& u)
    {
      return new ShellQuarkSinkSmearing<LatticeFermion>(ShellQuarkSinkSmearingParams(xml_in, path), u);
    }

    //! Name to be used
    const std::string name("SHELL_SINK");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= QuarkSmearingEnv::registered;
      foo &= Chroma::ThePropSinkSmearingFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSinkSmearingFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the sink smearing
    const bool registered = registerAll();
  }


  //! Smear the sink
  template<>
  void
  ShellQuarkSinkSmearing<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
  {
    QDPIO::cout << "Shell sink" << endl;

    try
    {
      //
      // Create the quark smearing object
      //
      std::istringstream  xml_s(params.quark_smearing);
      XMLReader  smeartop(xml_s);
      const string smear_path = "/SmearingParam";
	
      Handle< QuarkSmearing<LatticePropagator> >
	quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing_type,
								      smeartop,
								      smear_path));

      //
      // Sink smear quark
      //
      (*quarkSmearing)(quark_sink, u_smr);

      displacement(u_smr, quark_sink,
		   params.disp_length, params.disp_dir);

    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSinkSmearingEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }
  }



  //! Smear the sink
  template<>
  void
  ShellQuarkSinkSmearing<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
  {
    QDPIO::cout << "Shell sink" << endl;

    try
    {
      //
      // Create the quark smearing object
      //
      std::istringstream  xml_s(params.quark_smearing);
      XMLReader  smeartop(xml_s);
      const string smear_path = "/SmearingParam";
	
      Handle< QuarkSmearing<LatticeFermion> >
	quarkSmearing(TheFermSmearingFactory::Instance().createObject(params.quark_smearing_type,
								      smeartop,
								      smear_path));

      //
      // Sink smear quark
      //
      (*quarkSmearing)(quark_sink, u_smr);

      displacement(u_smr, quark_sink,
		   params.disp_length, params.disp_dir);

    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSinkSmearingEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }
  }

}
