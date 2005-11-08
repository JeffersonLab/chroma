// $Id: sh_source_smearing.cc,v 2.3 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sh_source_smearing.h"
#include "meas/sources/source_smearing_factory.h"

#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/quark_smearing_factory.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/displacement.h"

namespace Chroma
{
  //! Read parameters
  ShellQuarkSourceSmearingParams::ShellQuarkSourceSmearingParams()
  {
    disp_length = disp_dir = 0;
  }

  //! Read parameters
  ShellQuarkSourceSmearingParams::ShellQuarkSourceSmearingParams(XMLReader& xml, const string& path)
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

    read(paramtop, "SourceType",  source_type);

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
  void read(XMLReader& xml, const string& path, ShellQuarkSourceSmearingParams& param)
  {
    ShellQuarkSourceSmearingParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceSmearingParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "SourceType", param.source_type);
    xml << param.quark_smearing;
    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);

    pop(xml);
  }



  //! Hooks to register the class
  namespace ShellQuarkSourceSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
    {
      return new ShellQuarkSourceSmearing<LatticePropagator>(ShellQuarkSourceSmearingParams(xml_in, path), u);
    }

    //! Callback function
    QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						const std::string& path,
						const multi1d<LatticeColorMatrix>& u)
    {
      return new ShellQuarkSourceSmearing<LatticeFermion>(ShellQuarkSourceSmearingParams(xml_in, path), u);
    }

    //! Name to be used
    const std::string name("SHELL_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= QuarkSmearingEnv::registered;
      foo &= Chroma::ThePropSourceSmearingFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }



  //! Construct the source smearing
  void
  ShellQuarkSourceSmearing<LatticePropagator>::create()
  {
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);
  }

  //! Smear the source
  void
  ShellQuarkSourceSmearing<LatticePropagator>::operator()(LatticePropagator& quark_source) const
  {
    QDPIO::cout << "Shell source smearing" << endl;
 
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
      // Smear quark source
      //
      displacement(u_smr, quark_source,
		   (-1)*params.disp_length, params.disp_dir);

      (*quarkSmearing)(quark_source, u_smr);

    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSourceSmearingEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }
  }


  //! Construct the source smearing
  void
  ShellQuarkSourceSmearing<LatticeFermion>::create()
  {
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);
  }

  //! Construct the source smearing
  void
  ShellQuarkSourceSmearing<LatticeFermion>::operator()(LatticeFermion& quark_source) const
  {
    QDPIO::cout << "Shell source smearing" << endl;

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
      // Smear quark source
      //
      displacement(u_smr, quark_source,
		   (-1)*params.disp_length, params.disp_dir);

      (*quarkSmearing)(quark_source, u_smr);

    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << ShellQuarkSourceSmearingEnv::name << ": Caught Exception smearing: " << e << endl;
      QDP_abort(1);
    }
  }


}
