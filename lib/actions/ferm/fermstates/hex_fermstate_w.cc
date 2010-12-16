// -*- C++ -*-
// $Id: hex_fermstate_w.cc,v 1.2 2006-09-20 20:31:41 edwards Exp $
/*! @file 
 *  @brief Connection State for Hex state (.cpp file)
 */



#include "chromabase.h"
#include "actions/ferm/fermstates/hex_fermstate_w.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"



namespace Chroma 
{ 

  /*! \ingroup fermstates */
  namespace CreateHexFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
  multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
	const std::string& path) 
    {
      return new CreateHexFermState< LatticeFermion,
		                       multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, path),
	  HexFermStateParams(xml, path));
    }

    const std::string name = "HEX_FERM_STATE" ;

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


