// $Id: wall_source_const.cc,v 3.3 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/wall_source_const.h"
#include "meas/sources/walfil_w.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const string& path, WallQuarkSourceConstEnv::Params& param)
  {
    WallQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const WallQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace WallQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
      return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }
      
      //! Name to be used
      const std::string name("WALL_SOURCE");

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

      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_source", t_source);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "SourceType", WallQuarkSourceConstEnv::name);
      write(xml, "j_decay", j_decay);
      write(xml, "t_source", t_source);
      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Wall source" << endl;

      // Create the quark source
      LatticePropagator quark_source;

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // Wall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  walfil(chi, 
		 params.t_source,
		 params.j_decay, 
		 color_source, spin_source);
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      return quark_source;
    }

  }

}
