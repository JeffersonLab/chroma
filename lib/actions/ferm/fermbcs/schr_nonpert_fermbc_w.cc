// $Id: schr_nonpert_fermbc_w.cc,v 3.1 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - use for non-pertubative tuning of clover action
 */

#include "actions/ferm/fermbcs/schr_nonpert_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrNonPertFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrNonPertFermBC(SchrNonPertGaugeBC(SchrGaugeBCParams(xml, path)),
				   SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_NONPERT_FERMBC";

    static bool registered = false;

    // Register all objects
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
	registered = true;
      }
      return success;
    }
  }

}
