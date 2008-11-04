// $Id: pt_sink_smearing.cc,v 3.5 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Point sink construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/pt_sink_smearing.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, PointQuarkSinkSmearingEnv::Params& param)
  {
    PointQuarkSinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace PointQuarkSinkSmearingEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						     const std::string& path,
						     const multi1d<LatticeColorMatrix>& u)
      {
      return new SinkSmear<LatticePropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								  const std::string& path,
								  const multi1d<LatticeColorMatrix>& u)
      {
	return new SinkSmear<LatticeStaggeredPropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  const multi1d<LatticeColorMatrix>& u)
      {
	return new SinkSmear<LatticeFermion>(Params(xml_in, path), u);
      }

      //! Name to be used
      const std::string name("POINT_SINK");

      //! Local registration flag
      bool registered = false;
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
	success &= QuarkDisplacementEnv::registerAll();
	success &= Chroma::ThePropSinkSmearingFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSinkSmearingFactory::Instance().registerObject(name, createStagProp);
	success &= Chroma::TheFermSinkSmearingFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
      {
	quark_displacement = QuarkDisplacementEnv::nullXMLGroup();  // initialize
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

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
      int version = 1;
      write(xml, "version", version);

      write(xml, "SinkType", PointQuarkSinkSmearingEnv::getName());
      xml << link_smearing.xml;
      xml << quark_displacement.xml;
      pop(xml);
    }


    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
    {
      QDPIO::cout << "Point sink" << endl;
 
      try
      {
	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	//
	// Displace quark source
	//
	(*quarkDisplacement)(quark_sink, u_smr, PLUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in displacement: " << e << endl;
	QDP_abort(1);
      }
    }



    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark_sink) const
    {
      QDPIO::cout << "Point sink" << endl;
 
      try
      {
	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
	
	Handle< QuarkDisplacement<LatticeStaggeredPropagator> >
	  quarkDisplacement(TheStagPropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										    displacetop,
										    params.quark_displacement.path));
	
	//
	// Displace quark source
	//
	(*quarkDisplacement)(quark_sink, u_smr, PLUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in displacement: " << e << endl;
	QDP_abort(1);
      }
    }



    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
    {
//      QDPIO::cout << "Point sink" << endl;
 
      try
      {
	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
	
	Handle< QuarkDisplacement<LatticeFermion> >
	  quarkDisplacement(TheFermDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	//
	// Displace quark source
	//
	(*quarkDisplacement)(quark_sink, u_smr, PLUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in displacement: " << e << endl;
	QDP_abort(1);
      }
    }

  }

}
