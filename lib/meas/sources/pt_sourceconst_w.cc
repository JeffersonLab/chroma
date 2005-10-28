// $Id: pt_sourceconst_w.cc,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#include "chromabase.h"

#include "meas/sources/prop_source_factory_w.h"
#include "meas/sources/pt_sourceconst_w.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace PropPointSourceConstEnv
  {
    //! Callback function
    SourceConstruction<LatticePropagator>* createSource(XMLReader& xml_in,
							const std::string& path)
    {
      return new PropPointSourceConst(PointSourceConstParams(xml_in, path));
    }

    //! Name to be used
    const std::string name("POINT_SOURCE");

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createSource);
    }

    //! Register the source construction
    const bool registered = registerAll();
  }



  //! Construct the source
  LatticePropagator
  PropPointSourceConst::operator()(const multi1d<LatticeColorMatrix>& u) const
  {
    QDPIO::cout << "Point source" << endl;

    LatticePropagator quark_source;

    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      QDPIO::cout << "color = " << color_source << endl; 

      LatticeColorVector src_color_vec = zero;

      // Make a point source at coordinates t_source
      srcfil(src_color_vec, params.t_source, color_source);

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
