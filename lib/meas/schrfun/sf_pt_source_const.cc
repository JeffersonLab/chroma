// $Id: sf_pt_source_const.cc,v 1.1 2007-08-25 04:07:40 edwards Exp $
/*! \file
 *  \brief Point source construction for Schroedinger Functional
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/schrfun/sf_source_const_factory.h"
#include "meas/schrfun/sf_pt_source_const.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

namespace Chroma
{

  //! Hooks to register the class
  namespace PointSFSourceConstEnv
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
    const std::string name("POINT_SOURCE");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheSFSourceConstructionFactory::Instance().registerObject(name, createProp);
	registered = true;
      }
      return success;
    }

    //! Initialize
    Params::Params()
    {
      t_srce.resize(Nd-1);
      t_srce = 0;
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

      read(paramtop, "t_srce", t_srce);
      if (t_srce.size() != (Nd-1))
      {
	QDPIO::cerr << "SF shell source params: source location size is not Nd-1" << endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "SourceType", PointSFSourceConstEnv::name);
      write(xml, "t_srce", t_srce);
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

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	// Make a point source at coordinates t_source
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

      return quark_source;
    }
  }

}
