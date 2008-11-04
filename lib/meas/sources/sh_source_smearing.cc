// $Id: sh_source_smearing.cc,v 3.9 2008-11-04 18:43:59 edwards Exp $
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

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, ShellQuarkSourceSmearingEnv::Params& param)
  {
    ShellQuarkSourceSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace ShellQuarkSourceSmearingEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						     const std::string& path,
						     const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmearing<LatticePropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								  const std::string& path,
								  const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmearing<LatticeStaggeredPropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmearing<LatticeFermion>(Params(xml_in, path), u);
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("SHELL_SOURCE");
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LinkSmearingEnv::registerAll();
	success &= QuarkSmearingEnv::registerAll();
	success &= QuarkDisplacementEnv::registerAll();
	success &= Chroma::ThePropSourceSmearingFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSourceSmearingFactory::Instance().registerObject(name, createStagProp);
	success &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Read parameters
    Params::Params()
    {
      j_decay = -1;
      quark_smear_lastP = true;
    }

    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      quark_smear_lastP = true;

      switch (version) 
      {
      case 1:
      {
	quark_displacement.id = SimpleQuarkDisplacementEnv::getName();
	int disp_length = 0;
	int disp_dir = 0;

	XMLBufferWriter xml_tmp;
	push(xml_tmp, "Displacement");
	write(xml_tmp, "DisplacementType", quark_displacement.id);

	if (paramtop.count("disp_length") != 0)
	  read(paramtop, "disp_length", disp_length);

	if (paramtop.count("disp_dir") != 0)
	  read(paramtop, "disp_dir", disp_dir);

	write(xml_tmp, "disp_length", disp_length);
	write(xml_tmp, "disp_dir",  disp_dir);

	pop(xml_tmp);  // Displacement

	quark_displacement.xml = xml_tmp.printCurrentContext();
      }
      break;

      case 2:
      {
	if (paramtop.count("Displacement") != 0)
	  quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
	else
	  quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
      }
      break;

      case 3:
      {
	read(paramtop, "quark_smear_lastP", quark_smear_lastP);

	if (paramtop.count("Displacement") != 0)
	  quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
	else
	  quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      quark_smearing = readXMLGroup(paramtop, "SmearingParam", "wvf_kind");

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();

      read(paramtop, "j_decay",  j_decay);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 3;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", ShellQuarkSourceSmearingEnv::name);
      xml << quark_smearing.xml;
      xml << quark_displacement.xml;
      xml << link_smearing.xml;

      write(xml, "j_decay",  j_decay);
      write(xml, "quark_smear_lastP",  quark_smear_lastP);

      pop(xml);
    }


    //! Smear the source
    template<>
    void
    SourceSmearing<LatticePropagator>::operator()(LatticePropagator& quark_source) const
    {
//      QDPIO::cout << "Shell source smearing" << endl;
 
      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
//	const string smear_path = "/SmearingParam";
	
	Handle< QuarkSmearing<LatticePropagator> >
	  quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
//	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	// Smear and displace
	if (params.quark_smear_lastP)
	{
	  // Smear the colour source
	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);

	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);
	}
	else
	{
	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);

	  // Smear the colour source
	  // smear the point source first, then displace
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }




    //! Smear the source
    template<>
    void
    SourceSmearing<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark_source) const
    {
//      QDPIO::cout << "Shell source smearing" << endl;
 
      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
//	const string smear_path = "/SmearingParam";
	
	Handle< QuarkSmearing<LatticeStaggeredPropagator> >
	  quarkSmearing(TheStagPropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									    smeartop,
									    params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
//	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticeStaggeredPropagator> >
	  quarkDisplacement(TheStagPropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										    displacetop,
										    params.quark_displacement.path));

	// Smear and displace
	if (params.quark_smear_lastP)
	{
	  // Smear the colour source
	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);

	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);
	}
	else
	{
	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);

	  // Smear the colour source
	  // smear the point source first, then displace
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }



    //! Smear the source
    template<>
    void
    SourceSmearing<LatticeFermion>::operator()(LatticeFermion& quark_source) const
    {
//      QDPIO::cout << "Shell source smearing" << endl;
 
      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
//	const string smear_path = "/SmearingParam";
	
	Handle< QuarkSmearing<LatticeFermion> >
	  quarkSmearing(TheFermSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
//	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticeFermion> >
	  quarkDisplacement(TheFermDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	// Smear and displace
	if (params.quark_smear_lastP)
	{
	  // Smear the colour source
	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);

	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);
	}
	else
	{
	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);

	  // Smear the colour source
	  // smear the point source first, then displace
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }

  }
}
