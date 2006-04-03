// $Id: schr_dirich_fermbc_w.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
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
    FermBC<LatticeFermion,
	   multi1d<LatticeColorMatrix>, 
	   multi1d<LatticeColorMatrix> >* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrDirichletFermBC(SchrDirichletGaugeBC(SchrGaugeBCParams(xml, path)),
				     SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_DIRICHLET_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
