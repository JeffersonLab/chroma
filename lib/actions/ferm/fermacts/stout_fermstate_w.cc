// -*- C++ -*-
// $Id: stout_fermstate_w.cc,v 1.19 2006-09-06 18:59:18 bjoo Exp $
/*! @file 
 *  @brief Connection State for Stout state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/stout_fermstate_w.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

#include "actions/ferm/fermacts/ferm_createstate_factory_w.h"
#include "actions/ferm/fermacts/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"
#include "actions/ferm/fermacts/stout_utils.h"

namespace Chroma 
{ 

  /*! \ingroup fermacts */
  namespace CreateStoutFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreateStoutFermState< LatticeFermion,
		                       multi1d<LatticeColorMatrix>, 
		                       multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, path),
								     StoutFermStateParams(xml, path));
    }

    const std::string name = "STOUT_FERM_STATE";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, 
									  createFerm);
      return foo;
    }

    const bool registered = registerAll();
  }

  /*! \ingroup fermacts */
  namespace CreateSLICFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreateSLICFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, path),
				     StoutFermStateParams(xml, path));
    }

    const std::string name = "SLIC_FERM_STATE";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, 
									  createFerm);
      return foo;
    }

    const bool registered = registerAll();
  }





} // End namespace Chroma
