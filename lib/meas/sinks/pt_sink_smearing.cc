// $Id: pt_sink_smearing.cc,v 1.7 2006-02-22 04:34:05 edwards Exp $
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
    //! Callback function
    QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
    {
      return new SinkSmear<LatticePropagator>(Params(xml_in, path), u);
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

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= QuarkDisplacementEnv::registered;
      foo &= Chroma::ThePropSinkSmearingFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSinkSmearingFactory::Instance().registerObject(name, createFerm);
      return true;
    }

    //! Register the sink smearing
    const bool registered = registerAll();


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
	quark_displacement_type = SimpleQuarkDisplacementEnv::name;
	int disp_length = 0;
	int disp_dir = 0;

	XMLBufferWriter xml_tmp;
	push(xml_tmp, "Displacement");
	write(xml_tmp, "DisplacementType", quark_displacement_type);

	if (paramtop.count("disp_length") != 0)
	  read(paramtop, "disp_length", disp_length);

	if (paramtop.count("disp_dir") != 0)
	  read(paramtop, "disp_dir", disp_dir);

	write(xml_tmp, "disp_length", disp_length);
	write(xml_tmp, "disp_dir",  disp_dir);

	pop(xml_tmp);  // Displacement

	quark_displacement = xml_tmp.printCurrentContext();
      }
      break;

      case 2:
      {
	if (paramtop.count("Displacement") != 0)
	{
	  XMLReader xml_tmp(paramtop, "Displacement");
	  std::ostringstream os;
	  xml_tmp.print(os);
	  read(xml_tmp, "DisplacementType", quark_displacement_type);
	  quark_displacement = os.str();
	}
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      if (paramtop.count("LinkSmearing") != 0)
      {
	XMLReader xml_tmp(paramtop, "LinkSmearing");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "LinkSmearingType", link_smearing_type);
	link_smearing = os.str();
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
      int version = 1;
      write(xml, "version", version);

      xml << link_smearing;
      xml << quark_displacement;
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
	std::istringstream  xml_d(params.quark_displacement);
	XMLReader  displacetop(xml_d);
	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement_type,
										displacetop,
										displace_path));

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
	std::istringstream  xml_d(params.quark_displacement);
	XMLReader  displacetop(xml_d);
	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticeFermion> >
	  quarkDisplacement(TheFermDisplacementFactory::Instance().createObject(params.quark_displacement_type,
										displacetop,
										displace_path));

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
