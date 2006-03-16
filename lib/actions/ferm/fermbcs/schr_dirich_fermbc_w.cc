// $Id: schr_dirich_fermbc_w.cc,v 2.1 2006-03-16 03:00:12 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - dirichlet ferm BC
 */

#include "actions/ferm/fermbcs/schr_dirich_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrDirichletFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrDirichletFermBC<LatticeFermion>(SchrDirichletGaugeBC(SchrGaugeBCParams(xml, path)),
						     SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_DIRICHLET_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
