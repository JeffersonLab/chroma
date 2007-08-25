// $Id: sf_wall_source_const.cc,v 1.1 2007-08-25 04:07:41 edwards Exp $
/*! \file
 *  \brief Wall source construction for Schroedinger Functional
 */

#include "chromabase.h"

#include "meas/schrfun/sf_source_const_factory.h"
#include "meas/schrfun/sf_wall_source_const.h"
#include "meas/sources/walfil_w.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Hooks to register the class
  namespace WallSFSourceConstEnv
  {
    //! Callback function
    SFSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							const std::string& path)
    {
      return new SourceConst<LatticePropagator>(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name("WALL_SOURCE");

    //! Local registration flag
    static bool registered = false;

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
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "SourceType", WallSFSourceConstEnv::name);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u, 
					       int t_source, int decay_dir) const
    {
      QDPIO::cout << "Wall source" << endl;

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
      
      // Create the quark source
      LatticePropagator quark_source;

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // Wall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  walfil(chi, t_source, decay_dir, color_source, spin_source);
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      return quark_source;
    }

  }

}
