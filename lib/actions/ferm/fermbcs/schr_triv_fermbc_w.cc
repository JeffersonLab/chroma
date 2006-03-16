// $Id: schr_triv_fermbc_w.cc,v 2.1 2006-03-16 03:00:12 edwards Exp $
/*! \file
 *  \brief Schroedinger functional trivial ferm BC
 */

#include "actions/ferm/fermbcs/schr_triv_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrTrivialFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrTrivialFermBC<LatticeFermion>(SchrTrivialGaugeBC(SchrGaugeBCParams(xml, path)),
						   SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_TRIVIAL_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
