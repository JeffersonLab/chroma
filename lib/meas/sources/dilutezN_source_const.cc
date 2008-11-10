// $Id: dilutezN_source_const.cc,v 3.3 2008-11-10 00:04:58 kostas Exp $
/*! \file
 *  \brief Random ZN wall source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/quark_smearing_factory.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, DiluteZNQuarkSourceConstEnv::Params& param)
  {
    DiluteZNQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DiluteZNQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  // Hooks to register the class
  namespace DiluteZNQuarkSourceConstEnv
  {
    // Anonymous namespace
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SourceConst<LatticeFermion>(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("RAND_DILUTE_ZN_SOURCE");
    }  // end namespace

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
	//	success &= Chroma::TheFermSourceSmearingFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
      smear = false ;
      j_decay = -1;
      t_source = -1;
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

      smear = false ;
      if(paramtop.count("Smearing") !=0 ) {
	smr = readXMLGroup(paramtop, "Smearing", "wvf_kind");
	link_smear = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
	displace = readXMLGroup(paramtop, "Displacement","DisplacementType");
	smear = true ;
      }

      read(paramtop, "ran_seed", ran_seed);
      read(paramtop, "N", N);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_source", t_source);

      read(paramtop, "spatial_mask_size", spatial_mask_size);
      read(paramtop, "spatial_mask", spatial_mask);
      read(paramtop, "color_mask", color_mask);
      read(paramtop, "spin_mask", spin_mask);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;

      if(smear){
	 xml << smr.xml;
	 xml << displace.xml ;
	 xml << link_smear.xml ;
      }

      write(xml, "version", version);
      write(xml, "ran_seed", ran_seed);
      write(xml, "N", N);
      write(xml, "j_decay", j_decay);
      write(xml, "t_source", t_source);

      write(xml, "spatial_mask_size", spatial_mask_size);
      write(xml, "spatial_mask", spatial_mask);
      write(xml, "color_mask", color_mask);
      write(xml, "spin_mask", spin_mask);

      pop(xml);
    }



    //! Construct the source
    template<>
    LatticeFermion
    SourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Diluted random complex ZN source" << endl;

      //
      // Sanity checks
      //
      if (params.spatial_mask_size.size() != Nd-1)
      {
	QDPIO::cerr << name << ": spatial mask size incorrect 1" << endl;
	QDP_abort(1);
      }

      if (params.spatial_mask.size() == 0)
      {
	QDPIO::cerr << name << ": spatial mask incorrect 2" << endl;
	QDP_abort(1);
      }

      multi1d<int> lookup_dir(Nd-1);
      int mu = 0;
      for(int j=0; j < params.spatial_mask_size.size(); ++j, ++mu)
      {
	if (j == params.j_decay) ++mu;  // bump up to next dir

	lookup_dir[j] = mu;
      }

      for(int j=0; j < params.spatial_mask.size(); ++j)
      {
	if (params.spatial_mask[j].size() != Nd-1)
	{
	  QDPIO::cerr << name << ": spatial mask incorrect 3" << endl;
	  QDP_abort(1);
	}
      }

      // More sanity checks
      for(int c=0; c < params.color_mask.size(); ++c)
      {
	if (params.color_mask[c] < 0 || params.color_mask[c] >= Nc)
	{
	  QDPIO::cerr << name << ": color mask incorrect 6" << endl;
	  QDP_abort(1);
	}
      }

      for(int s=0; s < params.spin_mask.size(); ++s)
      {
	if (params.spin_mask[s] < 0 || params.spin_mask[s] >= Ns)
	{
	  QDPIO::cerr << name << ": spin mask incorrect 7" << endl;
	  QDP_abort(1);
	}
      }

      //
      // Finally, do something useful
      //

      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      // Set the seed to desired value
      QDP::RNG::setrn(params.ran_seed);

      // Create the noisy quark source on the entire lattice
      LatticeFermion quark_noise;
      zN_src(quark_noise, params.N);

      // This is the filtered noise source to return
      LatticeFermion quark_source = zero;

      // Filter over the color and spin indices first
      for(int s=0; s < params.spin_mask.size(); ++s)
      {
	int spin_source = params.spin_mask[s];
	LatticeColorVector colvec = peekSpin(quark_noise, spin_source);
	LatticeColorVector dest   = zero;

	for(int c=0; c < params.color_mask.size(); ++c)
	{ 
	  int color_source = params.color_mask[c];
	  LatticeComplex comp = peekColor(colvec, color_source);

	  pokeColor(dest, comp, color_source);
	}

	pokeSpin(quark_source, dest, spin_source);
      }

      quark_noise = quark_source;  // reset

      // Filter over the spatial sites
      LatticeBoolean    mask = false;  // this is the starting mask

      for(int n=0; n < params.spatial_mask.size(); ++n)
      {
	LatticeBoolean btmp = true;

	for(int j=0; j < params.spatial_mask[n].size(); ++j)
	  btmp &= (Layout::latticeCoordinate(lookup_dir[j]) % params.spatial_mask_size[j]) == params.spatial_mask[n][j];

	mask |= btmp;
      }

      // Filter over the time slices
      mask &= Layout::latticeCoordinate(params.j_decay) == params.t_source;

      // Zap the unused sites
      quark_source = where(mask, quark_noise, Fermion(zero));


      if(params.smear){// do the smearing
	Handle< QuarkSmearing<LatticeFermion> >  Smearing ;
	try{
	  std::istringstream  xml_l(params.smr.xml);
	  XMLReader  smrtop(xml_l);
	  QDPIO::cout << "Quark smearing type = " <<params.smr.id ; 
	  QDPIO::cout << endl;
	  
	  Smearing = 
	    TheFermSmearingFactory::Instance().createObject(params.smr.id,smrtop, 
							    params.smr.path);
	}
	catch(const std::string& e){
	  QDPIO::cerr <<name<< ": Caught Exception creating quark smearing object: " << e << endl;
	  QDP_abort(1);
	}
	catch(...){
	  QDPIO::cerr <<name<< ": Caught generic exception creating smearing object" << endl;
	  QDP_abort(1);
	}
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;

	try{
	  std::istringstream  xml_l(params.link_smear.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smear.id ; 
	  QDPIO::cout << endl;
	  
	  
	  Handle<LinkSmearing> linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smear.id, linktop, params.link_smear.path));
	  
	  (*linkSmearing)(u_smr);
	  //MesPlq(xml_out, "Smeared_Observables", u_smr);
	}
	catch(const std::string& e){
	  QDPIO::cerr<<name<<": Caught Exception link smearing: " << e << endl;
	  QDP_abort(1);
	}
	

	(*Smearing)(quark_source, u_smr);

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.displace.xml);
	XMLReader  displacetop(xml_d);
	
	try{
	  Handle< QuarkDisplacement<LatticeFermion> >
	    quarkDisplacement(TheFermDisplacementFactory::Instance().createObject(params.displace.id, displacetop, params.displace.path));
	  QDPIO::cout << "Quark displacement type = " << params.displace.id ; 
	  QDPIO::cout << endl;

	  // displacement has to be taken along negative direction.
	  // Not sure why MINUS....
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);
	}
	catch(const std::string& e){
	  QDPIO::cerr<<name<<": Caught Exception quark displacement: "<<e<< endl;
	  QDP_abort(1);
	}
      }// if(smear) ends here
      
      
      // Reset the seed
      QDP::RNG::setrn(ran_seed);

      return quark_source;
    }

  } // end namespace

} // end namespace Chroma
