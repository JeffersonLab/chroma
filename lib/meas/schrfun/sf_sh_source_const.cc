// $Id: sf_sh_source_const.cc,v 1.1 2007-08-25 04:07:40 edwards Exp $
/*! \file
 *  \brief Shell source construction within Schroedinger Functional
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/schrfun/sf_source_const_factory.h"
#include "meas/schrfun/sf_sh_source_const.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace ShellSFSourceConstEnv
  {
    namespace
    {
      //! Callback function
      SFSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    //! Name to be used
    const std::string name("SHELL_SOURCE");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= QuarkSmearingEnv::registerAll();
	success &= Chroma::TheSFSourceConstructionFactory::Instance().registerObject(name, createProp);

	registered = true;
      }
      return success;
    }


    //! Read parameters
    Params::Params()
    {
      t_srce.resize(Nd-1);
      t_srce = 0;
    }

    //! Read parameters
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "SF shell source: parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      quark_smearing = readXMLGroup(paramtop, "SmearingParam", "wvf_kind");
      read(paramtop, "t_srce", t_srce);

      if (t_srce.size() != (Nd-1))
      {
	QDPIO::cerr << "SF shell source params: source location size is not Nd-1" << endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", ShellSFSourceConstEnv::name);
      xml << quark_smearing.xml;

      write(xml, "t_srce",  t_srce);

      pop(xml);
    }
    


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u, 
					       int t_source, int decay_dir) const
    {
      QDPIO::cout << "Shell source" << endl;

      // Sanity check
      if (decay_dir < 0 || decay_dir >= Nd)
      {
	QDPIO::cerr << name << ": invalid decay_dir=" << decay_dir << endl;
	QDP_abort(1);
      }
      
      if (t_source < 0 || t_source >= QDP::Layout::lattSize()[decay_dir])
      {
	QDPIO::cerr << name << ": invalid time_slice=" << t_source << endl;
	QDP_abort(1);
      }
      
      //
      // Actual Nd source location
      //
      multi1d<int> t_srce(Nd);
      for(int mu=0,i=0; mu < Nd; ++mu)
      {
	if (mu == decay_dir)
	{
	  t_srce[mu] = t_source;
	}
	else
	{
	  t_srce[mu] = params.t_srce[i++];
	}
      }

      // Construct source
      LatticePropagator quark_source;

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
	// Create quark source
	//
	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << endl; 

	  LatticeColorVector src_color_vec = zero;

	  // Make a point source at coordinates t_srce
	  srcfil(src_color_vec, t_srce, color_source);

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

	// do the smearing
	(*quarkSmearing)(quark_source, u);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }

      return quark_source;
    }

  }
}
