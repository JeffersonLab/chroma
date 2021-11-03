// -*- C++ -*-
/*! @file 
 *  @brief Connection State for Hyp state (.cpp file)
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/hyp_fermstate_w.h"
#include "util/gauge/expmat.h"
#include "util/gauge/taproj.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma 
{ 

  /*! \ingroup fermstates */
  namespace CreateHypFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
                    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
                                                              const std::string& path) 
    {
      return new CreateHypFermState< LatticeFermion,
                                     multi1d<LatticeColorMatrix>, 
                                     multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, path),
                                                                   HypFermStateParams(xml, path));
    }
    
    const std::string name = "HYP_FERM_STATE" ;

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


