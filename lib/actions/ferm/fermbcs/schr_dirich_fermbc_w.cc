/*! \file
 *  \brief Schroedinger BC - dirichlet BC
 */

#include "actions/ferm/fermbcs/schr_dirich_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrDirichletFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml, const std::string& path)
    {
      return new SchrDirichletFermBC(SchrDirichletGaugeBC(SchrGaugeBCParams(xml, path)),
				     SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_DIRICHLET_FERMBC";

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
