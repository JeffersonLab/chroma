// $Id: pt_source_smearing.cc,v 2.8 2005-11-21 21:07:38 edwards Exp $
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

    write(xml, "disp_length", param.disp_length);
    write(xml, "disp_dir", param.disp_dir);
    xml << param.link_smearing;
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
      foo &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createFerm);
      return true;
    }

    //! Register the source smearing
    const bool registered = registerAll();


    //! Read parameters
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
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);
      int version = 1;
      write(xml, "version", version);

      write(xml, "disp_length", disp_length);
      write(xml, "disp_dir", disp_dir);
      xml << link_smearing;
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
      displacement(u_smr, quark_source,
		   (-1)*params.disp_length, params.disp_dir);
    }

  }  // end namespace

} // end namespace Chroma
