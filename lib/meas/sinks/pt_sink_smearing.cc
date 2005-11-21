// $Id: pt_sink_smearing.cc,v 1.6 2005-11-21 21:07:38 edwards Exp $
/*! \file
 *  \brief Point sink construction
 */

#include "chromabase.h"

#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/pt_sink_smearing.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/displacement.h"

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
      foo &= Chroma::ThePropSinkSmearingFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSinkSmearingFactory::Instance().registerObject(name, createFerm);
      return true;
    }

    //! Register the sink smearing
    const bool registered = registerAll();


    //! Initialize
    Params::Params()
    {
      disp_length = disp_dir = 0;
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
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
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


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
      int version = 1;
      write(xml, "version", version);

      write(xml, "disp_length", disp_length);
      write(xml, "disp_dir", disp_dir);
      xml << link_smearing;
      pop(xml);
    }


    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
    {
      QDPIO::cout << "Point sink" << endl;

      displacement(u_smr, quark_sink,
		   params.disp_length, params.disp_dir);
    }


    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
    {
//      QDPIO::cout << "Point sink" << endl;

      displacement(u_smr, quark_sink,
		   params.disp_length, params.disp_dir);
    }

  }
}
