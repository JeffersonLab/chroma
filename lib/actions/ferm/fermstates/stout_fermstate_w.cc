// -*- C++ -*-
// $Id: stout_fermstate_w.cc,v 1.2 2006-09-20 20:31:41 edwards Exp $
/*! @file 
 *  @brief Connection State for Stout state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma 
{ 

  /*! \ingroup fermstates */
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }

  /*! \ingroup fermstates */
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

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }





} // End namespace Chroma
