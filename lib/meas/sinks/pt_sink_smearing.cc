// $Id: pt_sink_smearing.cc,v 1.4 2005-11-08 18:32:29 edwards Exp $
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
  //! Initialize
  PointQuarkSinkSmearingParams::PointQuarkSinkSmearingParams()
  {
    disp_length = disp_dir = 0;
  }


  //! Read parameters
  PointQuarkSinkSmearingParams::PointQuarkSinkSmearingParams(XMLReader& xml, const string& path)
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

  // Read parameters
  void read(XMLReader& xml, const string& path, PointQuarkSinkSmearingParams& param)
  {
    PointQuarkSinkSmearingParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSinkSmearingParams& param)
  {
    push(xml, path);
    int version = 1;
    write(xml, "version", version);

    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);
  }


  //! Hooks to register the class
  namespace PointQuarkSinkSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
    {
      return new PointQuarkSinkSmearing<LatticePropagator>(PointQuarkSinkSmearingParams(xml_in, path), u);
    }

    //! Callback function
    QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						const std::string& path,
						const multi1d<LatticeColorMatrix>& u)
    {
      return new PointQuarkSinkSmearing<LatticeFermion>(PointQuarkSinkSmearingParams(xml_in, path), u);
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
  }


  //! Construct the sink smearing
  template<>
  void
  PointQuarkSinkSmearing<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
  {
    QDPIO::cout << "Point sink" << endl;

    displacement(u_smr, quark_sink,
		 params.disp_length, params.disp_dir);
  }



  //! Construct the sink smearing
  template<>
  void
  PointQuarkSinkSmearing<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
  {
    QDPIO::cout << "Point sink" << endl;

    displacement(u_smr, quark_sink,
		 params.disp_length, params.disp_dir);
  }

}
