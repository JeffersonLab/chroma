// $Id: rndzNwall_source_const.cc,v 1.3 2009-05-01 22:52:02 kostas Exp $
/*! \file
 *  \brief Random ZN wall source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/rndzNwall_source_const.h"
#include "util/ferm/transf.h"


#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, RandZNWallQuarkSourceConstEnv::Params& param)
  {
    RandZNWallQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const RandZNWallQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace RandZNWallQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }
      
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("RAND_ZN_WALL_SOURCE");
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
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      t_source = -1;
      N=4 ;
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

      read(paramtop, "ran_seed", ran_seed);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_source", t_source);
      read(paramtop, "N", N);
      
      if (paramtop.count("Displacement") != 0)
	quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
      else
	quark_displacement = QuarkDisplacementEnv::nullXMLGroup();

      read(paramtop, "quark_smear_lastP", quark_smear_lastP);
      quark_smearing = readXMLGroup(paramtop, "SmearingParam", "wvf_kind");

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
      write(xml, "SourceType", RandZNWallQuarkSourceConstEnv::name);
      xml << quark_smearing.xml;
      xml << quark_displacement.xml;
      xml << link_smearing.xml;
      write(xml, "version", version);
      write(xml, "ran_seed", ran_seed);
      write(xml, "j_decay", j_decay);
      write(xml, "t_source", t_source);
      write(xml, "N", N);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Rand Z"<<params.N<<" Wall source" << endl;

      
      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      // Set the seed to desired value
      QDP::RNG::setrn(params.ran_seed);
      
      // Create the quark source
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
	    QDPIO::cout << "Link smearing type = " << params.link_smearing.id << endl;
	    
	    Handle< LinkSmearing >
	      linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,  linktop,   params.link_smearing.path));
	    (*linkSmearing)(u_smr);
	  }
	  
	  //
	  // Create the quark smearing object
	  //
	  std::istringstream  xml_s(params.quark_smearing.xml);
	  XMLReader  smeartop(xml_s);
	  QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << endl;
	  
	  Handle< QuarkSmearing<LatticePropagator> >
	    quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing.id, smeartop, params.quark_smearing.path));
	  
	  //
	  // Create the quark displacement object
	  //
	  std::istringstream  xml_d(params.quark_displacement.xml);
	  XMLReader  displacetop(xml_d);
	  QDPIO::cout << "Displacement type = " << params.quark_displacement.id << endl;
	  
	  Handle< QuarkDisplacement<LatticePropagator> >
	    quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,	  displacetop,	  params.quark_displacement.path));
	  
	  multi1d<LatticeColorVector> tmp_color_vec(Nc);
	  
	  LatticeReal rnd,theta;
	  LatticeComplex z;

	  random(rnd);
	  
	  Real twopiN = Chroma::twopi / params.N;
	  theta = twopiN * floor(params.N*rnd);
	  z = cmplx(cos(theta),sin(theta));
	  
	  
	  for(int i=0; i<Nc; i++) {
	    tmp_color_vec[i] = zero;
	    pokeColor(tmp_color_vec[i], z, i);
	  }

	  for(int color_source = 0; color_source < Nc; ++color_source)
	    {
	      QDPIO::cout << "color = " << color_source << endl; 
	      int mu = params.j_decay;
	      int slice = params.t_source;
	      LatticeColorVector src_color_vec ;
	      if((slice>=Layout::lattSize()[mu])||(slice<0)){
		QDPIO::cout<<"Doing the full Volume source"<<endl ;
		src_color_vec = tmp_color_vec[color_source] ;
	      }
	      else{
		//QDPIO::cout<<"Doing Wall source on timeslice "<<slice<<endl ;
		src_color_vec = where( Layout::latticeCoordinate(mu) == slice,
				       tmp_color_vec[color_source],
				       LatticeColorVector(zero));
	      }
	      for(int spin_source = 0; spin_source < Ns; ++spin_source)
		{
		  QDPIO::cout << "spin = " << spin_source << endl; 
		  
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
	  if (params.quark_smear_lastP){
	    // Smear the colour source
	    // displace the point source first, then smear
	    // displacement has to be taken along negative direction.
	    (*quarkDisplacement)(quark_source, u_smr, MINUS);
	    
	    // do the smearing
	    (*quarkSmearing)(quark_source, u_smr);
	  }
	  else{
	    // do the smearing
	    (*quarkSmearing)(quark_source, u_smr);
	    
	    // Smear the colour source
	    // smear the point source first, then displace
	    // displacement has to be taken along negative direction.
	    (*quarkDisplacement)(quark_source, u_smr, MINUS);
	  }
	}// try
      catch(const std::string& e){
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }

      // Reset the seed
      QDP::RNG::setrn(ran_seed);

      return quark_source;
    }

  }

}
