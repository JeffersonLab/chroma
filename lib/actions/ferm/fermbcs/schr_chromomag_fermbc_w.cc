// $Id: schr_chromomag_fermbc_w.cc,v 2.1 2006-03-16 03:00:12 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - chromo-magnetic ferm BC
 */

#include "actions/ferm/fermbcs/schr_chromomag_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{

  namespace SchrChromoMagFermBCEnv 
  { 
    //! Callback function to register with the factory
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml, const string& path)
    {
      return new SchrChromoMagFermBC<LatticeFermion>(SchrChromoMagGaugeBC(SchrGaugeBCParams(xml, path)),
						     SchrFermBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_CHROMOMAG_FERMBC";
    const bool registered = TheWilsonTypeFermBCFactory::Instance().registerObject(name,
										  createFermBC);
  }

}
