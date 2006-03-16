// $Id: schr_nonpert_fermbc_w.cc,v 2.1 2006-03-16 03:00:12 edwards Exp $
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
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrNonPertFermBC<LatticeFermion>(SchrNonPertGaugeBC(SchrGaugeBCParams(xml, path)),
						   SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_NONPERT_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
