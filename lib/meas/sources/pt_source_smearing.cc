// $Id: pt_source_smearing.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_smearing_factory.h"
#include "meas/sources/pt_source_smearing.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, PointQuarkSourceSmearingEnv::Params& param)
  {
    PointQuarkSourceSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceSmearingEnv::Params& param)
  {
    push(xml, path);
    int version = 1;
    write(xml, "version", version);

    write(xml, "SourceType", PointQuarkSourceSmearingEnv::name);
    xml << param.link_smearing;
    xml << param.quark_displacement;
    pop(xml);
  }



  //! Hooks to register the class with the fermact factory
  namespace PointQuarkSourceSmearingEnv
  {
    //! Callback function
    QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						const std::string& path,
						const multi1d<LatticeColorMatrix>& u)
    {
      return new SourceSmear<LatticeFermion>(Params(xml_in, path), u);
    }

    //! Name to be used
    const std::string name("POINT_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= LinkSmearingEnv::registered;
      foo &= QuarkDisplacementEnv::registered;
      foo &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createFerm);
      return true;
    }

    //! Register the source smearing
    const bool registered = registerAll();


    //! Read parameters
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
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      if (paramtop.count("Displacement") != 0)
      {
	XMLReader xml_tmp(paramtop, "Displacement");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "DisplacementType", quark_displacement_type);
	quark_displacement = os.str();
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
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);
      int version = 1;
      write(xml, "version", version);

      xml << link_smearing;
      xml << quark_displacement;
      pop(xml);
    }



    //! Construct the source smearing
    template<>
    void
    SourceSmear<LatticeFermion>::operator()(LatticeFermion& quark_source) const
    {
//      QDPIO::cout << "Point source" << endl;

      // displace the point source first, then smear
      // displacement has to be taken along negative direction.
 
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
	(*quarkDisplacement)(quark_source, u_smr, MINUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }

  }  // end namespace

} // end namespace Chroma
