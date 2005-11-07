// $Id: pt_source_smearing.cc,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#include "chromabase.h"

#include "meas/sources/source_smearing_factory.h"
#include "meas/sources/pt_source_smearing.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/displacement.h"

namespace Chroma
{
  //! Read parameters
  PointQuarkSourceSmearingParams::PointQuarkSourceSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "disp_length", disp_length);
    read(paramtop, "disp_dir", disp_dir);

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
  void read(XMLReader& xml, const string& path, PointQuarkSourceSmearingParams& param)
  {
    PointQuarkSourceSmearingParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceSmearingParams& param)
  {
    push(xml, path);
    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
    pop(xml);
  }


  //! Hooks to register the class with the fermact factory
  namespace PointPropSourceSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
    {
      return new PointQuarkSourceSmearing<LatticePropagator>(PointQuarkSourceSmearingParams(xml_in, path), u);
    }

    //! Name to be used
    const std::string name("POINT_SINK");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= Chroma::ThePropSourceSmearingFactory::Instance().registerObject(name, createProp);
      return true;
    }

    //! Register the sink smearing
    const bool registered = registerAll();
  }


  //! Hooks to register the class with the fermact factory
  namespace PointFermSourceSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticeFermion>* createProp(XMLReader& xml_in,
						const std::string& path,
						const multi1d<LatticeColorMatrix>& u)
    {
      return new PointQuarkSourceSmearing<LatticeFermion>(PointQuarkSourceSmearingParams(xml_in, path), u);
    }
    
    //! Name to be used
    const std::string name("POINT_SINK");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createProp);
      return true;
    }

    //! Register the sink smearing
    const bool registered = registerAll();
  }


  //! Construct the source smearing
  void
  PointQuarkSourceSmearing<LatticePropagator>::create()
  {
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);
  }


  //! Construct the source smearing
  void
  PointQuarkSourceSmearing<LatticePropagator>::operator()(LatticePropagator& quark_source) const
  {
    QDPIO::cout << "Point source" << endl;

    // displace the point source first, then smear
    // displacement has to be taken along negative direction.
    displacement(u_smr, quark_source,
		 (-1)*params.disp_length, params.disp_dir);
  }



  //! Construct the source smearing
  void
  PointQuarkSourceSmearing<LatticeFermion>::create()
  {
    linkSmear(u_smr, params.link_smearing, params.link_smearing_type);
  }


  //! Construct the source smearing
  void
  PointQuarkSourceSmearing<LatticeFermion>::operator()(LatticeFermion& quark_source) const
  {
    QDPIO::cout << "Point source" << endl;

    // displace the point source first, then smear
    // displacement has to be taken along negative direction.
    displacement(u_smr, quark_source,
		 (-1)*params.disp_length, params.disp_dir);
  }

}
