// $Id: norm_sh_sink_smearing.cc,v 3.1 2009-04-12 03:45:00 kostas Exp $
/*! \file
 *  \brief NormShell sink smearing
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/norm_sh_sink_smearing.h"

#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/quark_smearing.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, NormShellQuarkSinkSmearingEnv::Params& param)
  {
    NormShellQuarkSinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const NormShellQuarkSinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace NormShellQuarkSinkSmearingEnv
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
      const std::string name("NORM_SHELL_SINK");

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
	success &= QuarkSmearingEnv::registerAll();
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
      quark_smear_firstP = true;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      quark_smear_firstP = true;

      switch (version) 
      {
      case 1:
      {
	quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
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
	read(paramtop, "quark_smear_firstP", quark_smear_firstP);

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

      read(paramtop, "Nhits",  Nhits);
      
      site_normalized = false ;
      if (paramtop.count("SiteNormalize") != 0){
        read(paramtop, "SiteNormalize", site_normalized );
      }
      
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 3;
      write(xml, "version", version);

      write(xml, "SinkType", NormShellQuarkSinkSmearingEnv::getName());
      xml << quark_smearing.xml;
      xml << quark_displacement.xml;
      xml << link_smearing.xml;
      write(xml, "quark_smear_firstP", quark_smear_firstP);
      write(xml, "Nhits",  Nhits);
      write(xml, "SiteNormalize",  site_normalized);
      pop(xml);

      pop(xml);
    }



    //! Smear the sink
    template<>
    void
    SinkSmear<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
    {
      QDPIO::cout << "NormShell sink" << endl;

      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
        QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << endl;
	
	Handle< QuarkSmearing<LatticePropagator> >
	  quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
        QDPIO::cout << "Displacement type = " << params.quark_displacement.id << endl;
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	if (params.quark_smear_firstP)
	{
	  //
	  // Sink smear quark
	  //
	  for(int hit(0);hit<params.Nhits;hit++){
	    (*quarkSmearing)(quark_sink, u_smr);

	    for(int src_s(0);src_s<Ns;src_s++){
	      for(int snk_s(0);snk_s<Ns;snk_s++){
		LatticeColorMatrix cm = peekSpin(quark_sink,snk_s,src_s);
		multi1d<LatticeColorVector> v(Nc) ;
		for(int c(0);c<Nc;c++){
		  for(int i = 0; i < Nc; ++i)
		    pokeColor(v[c], peekColor(cm, i, c), i) ;

		  // now normalize the vectors
		  for(int k(0);k<c;k++){
		    DComplex cc = innerProduct(v[k],v[c]);
		    v[c] -= cc*v[k] ;
		  }
		  Double nn = sqrt(norm2(v[c])) ;
		  v[c] /= nn ;
		  for(int i = 0; i < Nc; ++i)
		    pokeColor(cm,peekColor(v[c],i),i,c);
		}//c
		pokeSpin(quark_sink,cm,snk_s,src_s);
	      }//snk_s
	    }//src_s
	  }
	  //
	  // Displace the sink
	  //
	  (*quarkDisplacement)(quark_sink, u_smr, PLUS);
	}
	else
	{
	  QDPIO::cout<<__func__<<": WARNING! smearing last not implemented...\n" ;
	  //
	  // Displace the sink
	  //
	  (*quarkDisplacement)(quark_sink, u_smr, PLUS);

	  //
	  // Sink smear quark
	  //
	  (*quarkSmearing)(quark_sink, u_smr);
	}

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }




    //! Smear the sink
    template<>
    void
    SinkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark_sink) const
    {
      QDPIO::cout << "NormShell sink" << endl;

      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
        QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << endl;
	
	Handle< QuarkSmearing<LatticeStaggeredPropagator> >
	  quarkSmearing(TheStagPropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									    smeartop,
									    params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
        QDPIO::cout << "Displacement type = " << params.quark_displacement.id << endl;
	
	Handle< QuarkDisplacement<LatticeStaggeredPropagator> >
	  quarkDisplacement(TheStagPropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										    displacetop,
										    params.quark_displacement.path));

	if (params.quark_smear_firstP)
	{
	  //
	  // Sink smear quark
	  //
	  (*quarkSmearing)(quark_sink, u_smr);

	  //
	  // Displace the sink
	  //
	  (*quarkDisplacement)(quark_sink, u_smr, PLUS);
	}
	else
	{
	  //
	  // Displace the sink
	  //
	  (*quarkDisplacement)(quark_sink, u_smr, PLUS);

	  //
	  // Sink smear quark
	  //
	  (*quarkSmearing)(quark_sink, u_smr);
	}

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }



    //! Smear the sink
    template<>
    void
    SinkSmear<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
    {
//      QDPIO::cout << "NormShell sink" << endl;

      try
      {
	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
	
	Handle< QuarkSmearing<LatticeFermion> >
	  quarkSmearing(TheFermSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

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
	// Sink smear quark
	//
	(*quarkSmearing)(quark_sink, u_smr);


	//
	// Displace the sink
	//
	(*quarkDisplacement)(quark_sink, u_smr, PLUS);

      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }

  }

}
