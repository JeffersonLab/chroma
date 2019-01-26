/*! \file
 *  \brief Shell source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sh_zN_grid_source_const.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"
#include "meas/sources/zN_src.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const std::string& path, ShellZnGridQuarkSourceConstEnv::Params& param)
  {
    ShellZnGridQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const ShellZnGridQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace ShellZnGridQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }

      //! Callback function
      QuarkSourceConstruction<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
									  const std::string& path)
      {
	return new SourceConst<LatticeStaggeredPropagator>(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("ZNC_GRID_SHELL_SOURCE");
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
	success &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSourceConstructionFactory::Instance().registerObject(name, createStagProp);

	registered = true;
      }
      return success;
    }


    //! Read parameters
    Params::Params()
    {
      j_decay = -1;
      t_srce.resize(Nd);
      t_srce = 0;
      quark_smear_lastP = true;
    }

    //! Read parameters
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      quark_smear_lastP = true;

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
	read(paramtop, "quark_smear_lastP", quark_smear_lastP);

	if (paramtop.count("Displacement") != 0)
	  quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
	else
	  quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }

      quark_smearing = readXMLGroup(paramtop, "SmearingParam", "wvf_kind");

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();

      read(paramtop, "t_srce", t_srce);
      read(paramtop, "grd", grd);
      read(paramtop, "j_decay",  j_decay);
      read(paramtop, "ran_seed",  ran_seed);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 3;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", ShellZnGridQuarkSourceConstEnv::name);
      xml << quark_smearing.xml;
      xml << quark_displacement.xml;
      xml << link_smearing.xml;

      write(xml, "t_srce",  t_srce);
      write(xml, "j_decay",  j_decay);
      write(xml, "quark_smear_lastP",  quark_smear_lastP);
      write(xml, "grd",  grd);
      write(xml, "ran_seed",  ran_seed);

      pop(xml);
    }
    


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "ZNc random Shell source" << std::endl;

      Seed curr_seed;
      QDP::RNG::savern(curr_seed);

      //Seed the random number generator
      QDP::RNG::setrn(params.ran_seed);

      LatticePropagator quark_source;

      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	{
	  std::istringstream  xml_l(params.link_smearing.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smearing.id << std::endl;
	
	  Handle< LinkSmearing >
	    linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
									 linktop,
									 params.link_smearing.path));
	  (*linkSmearing)(u_smr);
	}

	//
	// Create the quark smearing object
	//
	LatticeComplex noise ;
	
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
        QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << std::endl;
	
	Handle< QuarkSmearing<LatticePropagator> >
	  quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
        QDPIO::cout << "Displacement type = " << params.quark_displacement.id << std::endl;
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));

	//
	// Create quark source
	//

	//setup mask
	LatticeBoolean mask = false ;
	LatticeBoolean btmp = true;
	for(int j=0; j < params.grd.size(); ++j){
	  LatticeInteger X= Layout::latticeCoordinate(j) - params.t_srce[j];
	  if(params.grd[j] > 0)
	    X = X % params.grd[j];
	  btmp &= (X==0) ;
	}
	mask |= btmp;
	
	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << std::endl; 

	  

	  // Make a point source at coordinates t_srce
	  //srcfil(src_color_vec, params.t_srce, color_source);

	  for(int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin = " << spin_source << std::endl; 

	    LatticeColorVector src_color_vec = zero;
	    LatticeColorVector cv=zero ;
	    LatticeComplex ZNc = zN_rng(Nc);
	    pokeColor(cv,ZNc,color_source) ;
	    src_color_vec = where(mask,cv,LatticeColorVector(zero));
	    // Insert a ColorVector into spin index spin_source
	    // This only overwrites sections, so need to initialize first
	    LatticeFermion chi = zero;

	    CvToFerm(src_color_vec, chi, spin_source);
      
	    /*
	     *  Move the source to the appropriate components
	     *  of quark source.
	     */
	    FermToProp(chi, quark_source, color_source, spin_source);
	  }
	}


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
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << std::endl;
	QDP_abort(1);
      }

      //Return the seed to its previous value
      QDP::RNG::setrn(curr_seed);
      
      return quark_source;
    }



    //! Construct the source
    template<>
    LatticeStaggeredPropagator
    SourceConst<LatticeStaggeredPropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Shell source" << std::endl;

      LatticeStaggeredPropagator quark_source;
      QDPIO::cout << "IT ONLY does one point here. I will not bother with staggered"<<std::endl;
      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	{
	  std::istringstream  xml_l(params.link_smearing.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smearing.id << std::endl;
	
	  Handle< LinkSmearing >
	    linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
									 linktop,
									 params.link_smearing.path));
	  (*linkSmearing)(u_smr);
	}

	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
        QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << std::endl;
	
	Handle< QuarkSmearing<LatticeStaggeredPropagator> >
	  quarkSmearing(TheStagPropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									    smeartop,
									    params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
        QDPIO::cout << "Displacement type = " << params.quark_displacement.id << std::endl;
	
	Handle< QuarkDisplacement<LatticeStaggeredPropagator> >
	  quarkDisplacement(TheStagPropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										    displacetop,
										    params.quark_displacement.path));

	//
	// Create quark source
	//
	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << std::endl; 

	  LatticeColorVector src_color_vec = zero;

	  // Make a point source at coordinates t_srce
	  srcfil(src_color_vec, params.t_srce, color_source);

	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  LatticeStaggeredFermion chi = zero;

	  CvToFerm(src_color_vec, chi);
      
	  /*
	   *  Move the source to the appropriate components
	   *  of quark source.
	   */
	  FermToProp(chi, quark_source, color_source);
	}


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
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << std::endl;
	QDP_abort(1);
      }

      return quark_source;
    }

  }
}
