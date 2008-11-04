// $Id: pt_source_smearing.cc,v 3.6 2008-11-04 18:43:58 edwards Exp $
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

#include "meas/smear/simple_quark_displacement.h"

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
    param.writeXML(xml, path);
  }



  //! Hooks to register the class with the fermact factory
  namespace PointQuarkSourceSmearingEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						     const std::string& path,
						     const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmear<LatticePropagator>(Params(xml_in, path), u);
      }
      
      //! Callback function
      QuarkSourceSink<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								  const std::string& path,
								  const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmear<LatticeStaggeredPropagator>(Params(xml_in, path), u);
      }
      
      //! Callback function
      QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  const multi1d<LatticeColorMatrix>& u)
      {
	return new SourceSmear<LatticeFermion>(Params(xml_in, path), u);
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("POINT_SOURCE");
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

      read(paramtop, "j_decay", j_decay);

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 2;
      write(xml, "version", version);

      write(xml, "SourceType", PointQuarkSourceSmearingEnv::name);
      write(xml, "j_decay", j_decay);
      xml << link_smearing.xml;
      xml << quark_displacement.xml;
      pop(xml);
    }



    //! Construct the source smearing
    template<>
    void
    SourceSmear<LatticePropagator>::operator()(LatticePropagator& quark_source) const
    {
//      QDPIO::cout << "Point source" << endl;

      // displace the point source first, then smear
      // displacement has to be taken along negative direction.
 
      try
      {
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



    //! Construct the source smearing
    template<>
    void
    SourceSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark_source) const
    {
//      QDPIO::cout << "Point source" << endl;

      // displace the point source first, then smear
      // displacement has to be taken along negative direction.
 
      try
      {
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
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
//	const string displace_path = "/Displacement";
	
	Handle< QuarkDisplacement<LatticeFermion> >
	  quarkDisplacement(TheFermDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

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
