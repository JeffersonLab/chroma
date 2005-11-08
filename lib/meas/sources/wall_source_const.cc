// $Id: wall_source_const.cc,v 2.1 2005-11-08 05:29:02 edwards Exp $
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
  //! Initialize
  WallQuarkSourceConstParams::WallQuarkSourceConstParams()
  {
    j_decay = -1;
    t_source = -1;
  }


  //! Read parameters
  WallQuarkSourceConstParams::WallQuarkSourceConstParams(XMLReader& xml, const string& path)
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

  // Read parameters
  void read(XMLReader& xml, const string& path, WallQuarkSourceConstParams& param)
  {
    WallQuarkSourceConstParams tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const WallQuarkSourceConstParams& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "j_decay", param.j_decay);
    write(xml, "t_source", param.t_source);
    pop(xml);
  }



  //! Hooks to register the class
  namespace WallQuarkSourceConstEnv
  {
    //! Callback function
    QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							   const std::string& path)
    {
      return new WallQuarkSourceConst<LatticePropagator>(WallQuarkSourceConstParams(xml_in, path));
    }

    //! Callback function
    QuarkSourceConstruction<LatticeFermion>* createFerm(XMLReader& xml_in,
							const std::string& path)
    {
      return new WallQuarkSourceConst<LatticeFermion>(WallQuarkSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("WALL_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSourceConstructionFactory::Instance().registerObject(name, createFerm);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Construct the source
  LatticePropagator
  WallQuarkSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
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



  //! Construct the source
  LatticeFermion
  WallQuarkSourceConst<LatticeFermion>::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Wall source" << endl;

    // Create the quark source
    LatticeFermion quark_source = zero;

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	walfil(quark_source, 
	       params.t_source,
	       params.j_decay, 
	       color_source, spin_source);
      }
    }

    return quark_source;
  }

}
